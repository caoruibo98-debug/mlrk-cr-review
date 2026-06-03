"""
ltr_candidates.py — 步骤1+2：同底物 hard 负样本 + 反应类型特征

对每个底物 S，用模块规则 RunReactants 生成"同底物候选产物"，每个候选带其生成规则的 EC class
（=反应类型特征）。真产物注入候选集；其余 = 同底物 hard decoy（结构真实但错的变换）。

这把任务从"词表检索(相似度主导)"换成"同底物候选判别(变换主导)"——逼模型学'这个底物走哪个反应'，
而不是'哪个产物最像底物'。这才能跑赢 Tanimoto。

输出 modular/outputs/ltr/candidates.parquet
用法：python scripts/ssrf/ml_ranking_kernel/modular/ltr_candidates.py [--n-rules 600] [--max-decoy 100] [--weak-cap 1200] [--force]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio  # noqa: E402
from generate_candidates import run_reactants  # 复用 RunReactants + 标准化  # noqa: E402
from ltr_build import LTR_OUT  # noqa: E402

RDLogger.DisableLog("rdApp.*")
POOL = r"D:\CRB\FoodGut\positive_sample_db\data\processed\foodgut_modular_positive_candidate_pool_v2.csv"
MACRO2MOD = {"small_molecule_polyphenol": "A", "carbohydrate_glycan": "B",
             "protein_amino_acid": "C", "lipid_fat": "D"}


def ec_class(ec) -> int:
    s = "" if pd.isna(ec) else str(ec)
    return int(s[0]) if s[:1].isdigit() else 0


def load_module_rules(rng, n_rules):
    d = pd.read_csv(POOL, low_memory=False, usecols=["macro_module", "reaction_rule_smarts", "enzyme_ec"])
    d["module"] = d["macro_module"].map(MACRO2MOD)
    d = d.dropna(subset=["module", "reaction_rule_smarts"])
    d = d[d["reaction_rule_smarts"].astype(str).str.contains(">>")]
    d["ecc"] = d["enzyme_ec"].map(ec_class)
    out = {}
    for m, g in d.groupby("module"):
        gg = g.drop_duplicates("reaction_rule_smarts")[["reaction_rule_smarts", "ecc"]]
        if len(gg) > n_rules:
            gg = gg.sample(n=n_rules, random_state=int(rng.integers(0, 2**31 - 1)))
        out[m] = list(gg.itertuples(index=False, name=None))
    return out


def main() -> None:
    from rdkit import Chem
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-rules", type=int, default=600)
    ap.add_argument("--max-decoy", type=int, default=100)
    ap.add_argument("--weak-cap", type=int, default=1200)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    rng = np.random.default_rng(0)

    rxn = pd.read_parquet(LTR_OUT / "reactions.parquet")
    rxn["ecc"] = rxn["ev_ec_class"] if "ev_ec_class" in rxn.columns else 0
    true_prod = rxn.groupby(["module", "sb"])["pb"].agg(set).to_dict()
    true_ecc = rxn.groupby(["module", "sb"]).agg(ecc=("ecc", "max")).to_dict()["ecc"]
    rules = load_module_rules(rng, args.n_rules)
    kio.log.info("module rules: %s", {m: len(v) for m, v in rules.items()})

    rows = []
    for m in ["A", "B", "C", "D"]:
        cols = ["sb", "substrate_smiles", "substrate_scaffold", "is_clean", "year"]
        if "tier" in rxn.columns:
            cols.append("tier")
        sub = rxn[rxn["module"] == m].drop_duplicates("sb")[cols].copy()
        if "tier" not in sub.columns:
            sub["tier"] = np.where(sub["is_clean"], "gold", "weak")
        # gold+silver 全留（=测试+anchor）；weak 子采样（仅训练 breadth）
        prio = sub[sub["tier"].isin(["gold", "silver"])]
        weak = sub[sub["tier"] == "weak"]
        if len(weak) > args.weak_cap:
            weak = weak.sample(n=args.weak_cap, random_state=0)
        sub = pd.concat([prio, weak]).drop_duplicates("sb")
        kio.log.info("module %s: %d substrates × %d rules", m, len(sub), len(rules[m]))
        for si, (_, r) in enumerate(sub.iterrows()):
            s = r["substrate_smiles"]; sb = r["sb"]
            cmol = Chem.MolFromSmiles(s)
            if cmol is None:
                continue
            cmolh = Chem.AddHs(cmol)
            truth = true_prod.get((m, sb), set())
            cand = {}  # pb -> (smiles, ecc)
            for smarts, ecc in rules[m]:
                for psmi in run_reactants(smarts, cmol, cmolh):
                    pik = kio.smiles_to_inchikey(psmi)
                    pb = kio.inchikey_block1(pik)
                    if not pb or pb == sb:
                        continue
                    if pb not in cand:
                        cand[pb] = (psmi, ecc)
            # decoy = 非真产物，采样上限
            decoys = [(pb, v[0], v[1]) for pb, v in cand.items() if pb not in truth]
            if len(decoys) > args.max_decoy:
                idx = rng.choice(len(decoys), args.max_decoy, replace=False)
                decoys = [decoys[i] for i in idx]
            # 注入真产物（保证可被排序）
            tp_rows = rxn[(rxn.module == m) & (rxn.sb == sb)][["pb", "product_smiles", "ecc"]].values
            recs = [(pb, psm, ecc, 1) for pb, psm, ecc in tp_rows]
            recs += [(pb, psm, ecc, 0) for pb, psm, ecc in decoys]
            for pb, psm, ecc, y in recs:
                rows.append({"module": m, "sb": sb, "substrate_smiles": s,
                             "substrate_scaffold": r["substrate_scaffold"],
                             "pb": pb, "product_smiles": psm, "ec_class": int(ecc), "y": int(y),
                             "is_clean": bool(r["is_clean"]), "tier": r.get("tier", "weak"), "year": r["year"]})
            if (si + 1) % 300 == 0:
                kio.log.info("  %s %d/%d subs, rows=%d", m, si + 1, len(sub), len(rows))

    cand = pd.DataFrame(rows)
    # 每底物至少要有 1 正 1 负才有意义
    grp = cand.groupby("sb")["y"].agg(["sum", "count"])
    keep_subs = set(grp[(grp["sum"] >= 1) & (grp["count"] - grp["sum"] >= 1)].index)
    cand = cand[cand["sb"].isin(keep_subs)].copy()
    kio.write_table(cand, LTR_OUT / "candidates.parquet", force=args.force)
    kio.log.info("== candidates: %d rows, %d substrates (pos=%d decoy=%d) ==",
                 len(cand), cand["sb"].nunique(), int((cand.y == 1).sum()), int((cand.y == 0).sum()))
    kio.log.info("median candidates/substrate: %.0f | clean substrates with usable set: %d",
                 cand.groupby("sb").size().median(), cand[cand.is_clean]["sb"].nunique())


if __name__ == "__main__":
    main()
