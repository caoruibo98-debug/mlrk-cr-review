"""
module_router.py — Phase 0

把任意 SMILES 路由到四个食品组分模块之一：
  B = 碳水/糖苷 (carbohydrate / glycoside)
  C = 蛋白/肽/氨基酸 (protein / peptide / amino acid)
  D = 脂肪/胆汁/固醇 (lipid / fatty acyl / bile / steroid)
  A = 小分子（其余，多酚等）

判定全靠 RDKit 结构（离线，非打分路径）。优先级按"肠道菌首步反应"逻辑：
  氨基酸/肽(C) > 糖苷(B) > 脂质/固醇(D) > 小分子(A)。
（糖苷黄酮如 rutin 首步是脱糖→归 B；游离氨基酸→C；胆汁酸/固醇→D；多酚苷元→A）

CLI：在 benchmark(67, 有 module 真值)上报路由准确率 + 混淆矩阵。
  python scripts/ssrf/ml_ranking_kernel/modular/module_router.py
"""
from __future__ import annotations

import sys
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import kio  # noqa: E402

RDLogger.DisableLog("rdApp.*")

MODULES = ["A", "B", "C", "D"]
MODULE_NAMES = {"A": "small_molecule", "B": "carb_glycoside",
                "C": "protein_peptide_aa", "D": "lipid_bile_steroid"}

# --- SMARTS ---
_AA = Chem.MolFromSmarts("[NX3;H1,H2;!$(NC=O)][CX4H]([!$([CX3]=O)])[CX3](=O)[OX1H0-,OX2H1]")
_PEPTIDE = Chem.MolFromSmarts("[NX3;H1][CX4][CX3](=O)[NX3;H1][CX4][CX3](=O)")  # ≥2 肽键骨架
_GLYCOSIDIC = Chem.MolFromSmarts("[C;R]([O;R])[O;!R][#6]")  # 端基碳-环O + 外接O-C（糖苷键）
_ESTER = Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]")
_ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")


def _has(mol, patt) -> bool:
    return patt is not None and mol.HasSubstructMatch(patt)


def _has_steroid_core(mol) -> bool:
    """甾核 = 6-6-6-5 全碳稠环系。环分析：≥3 个全碳六元环 + ≥1 个全碳五元环，且互相稠合。"""
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    def all_carbon(r):
        return all(mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in r)
    six = [set(r) for r in rings if len(r) == 6 and all_carbon(r)]
    five = [set(r) for r in rings if len(r) == 5 and all_carbon(r)]
    if len(six) < 3 or len(five) < 1:
        return False
    # 稠合：构造环邻接（共享 ≥2 原子），看是否有含≥3六+≥1五的连通块
    allr = six + five
    n = len(allr)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if len(allr[i] & allr[j]) >= 2:
                adj[i].add(j); adj[j].add(i)
    seen = set()
    for s in range(n):
        if s in seen:
            continue
        comp = []; stack = [s]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x); comp.append(x)
            stack.extend(adj[x] - seen)
        n6 = sum(1 for k in comp if k < len(six))
        n5 = sum(1 for k in comp if k >= len(six))
        if n6 >= 3 and n5 >= 1:
            return True
    return False


def _sugar_ring(mol) -> bool:
    """5/6 元环：恰 1 个环氧 + 其余环碳，且 ≥2 个环碳带 -OH（吡喃/呋喃糖）。"""
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) not in (5, 6):
            continue
        o = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "O")
        c = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "C")
        if o != 1 or c != len(ring) - 1:
            continue
        oh = 0
        for i in ring:
            a = mol.GetAtomWithIdx(i)
            if a.GetSymbol() != "C":
                continue
            for nb in a.GetNeighbors():
                if nb.GetSymbol() == "O" and nb.GetIdx() not in ring and nb.GetTotalNumHs() >= 1:
                    oh += 1
                    break
        if oh >= 2:
            return True
    return False


def _longest_carbon_chain(mol) -> int:
    """最长全碳链长度（粗略：碳-碳子图最长路径，用 BFS 近似）。"""
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C" and not a.GetIsAromatic()]
    cset = set(carbons)
    best = 0
    for start in carbons:
        # BFS 最远（链近似；带支链会高估，作粗筛足够）
        seen = {start}
        stack = [(start, 1)]
        while stack:
            node, d = stack.pop()
            best = max(best, d)
            for nb in mol.GetAtomWithIdx(node).GetNeighbors():
                j = nb.GetIdx()
                if j in cset and j not in seen and nb.GetDegree() <= 3:
                    seen.add(j)
                    stack.append((j, d + 1))
    return best


def classify(smiles: str) -> tuple[str, str]:
    """返回 (module, reason)。"""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) and smiles else None
    if mol is None:
        return "A", "unparseable->default_A"
    # C 蛋白/肽/氨基酸
    if _has(mol, _PEPTIDE):
        return "C", "peptide_backbone"
    if _has(mol, _AA):
        return "C", "free_amino_acid"
    # B 碳水/糖苷
    if _has(mol, _GLYCOSIDIC) or _sugar_ring(mol):
        return "B", "glycoside_or_sugar_ring"
    # D 脂肪/胆汁/固醇
    if _has_steroid_core(mol):
        return "D", "steroid_core"
    if (_has(mol, _ESTER) or _has(mol, _ACID)) and _longest_carbon_chain(mol) >= 8:
        return "D", "long_chain_acyl_ester"
    # A 其余
    return "A", "small_molecule_default"


def _bench_module_to_letter(m: str) -> str:
    m = str(m).upper()
    if m.startswith("A"):
        return "A"
    if m.startswith("B"):
        return "B"
    if m.startswith("C"):
        return "C"
    if m.startswith("D"):
        return "D"
    return "?"


def main() -> None:
    import pandas as pd
    kio.setup_logging()
    b = pd.read_csv(kio.resolve_input("analysis/benchmark_ssrf_rclss_2026-05-22.csv"))
    b = b.dropna(subset=["substrate_smiles", "module"]).copy()
    b["true"] = b["module"].map(_bench_module_to_letter)
    b["pred"] = b["substrate_smiles"].map(lambda s: classify(s)[0])
    acc = (b["true"] == b["pred"]).mean()
    kio.log.info("router accuracy on benchmark (n=%d): %.1f%%", len(b), 100 * acc)
    print("\n混淆矩阵 (行=真, 列=预测):")
    print(pd.crosstab(b["true"], b["pred"], dropna=False).to_string())
    print("\n错分样例:")
    wrong = b[b["true"] != b["pred"]][["substrate_name", "true", "pred", "food_component_class"]]
    print(wrong.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
