"""
ltr_clean_candidates.py — 杜绝作弊的候选生成（部署忠实）

与之前的区别（堵死三类作弊）：
  - 候选**只用规则 RunReactants 生成**（部署时本就如此）；不注入"没生成出来的真产物"，不加 vocab-NN decoy。
  - 每个候选带**它生成规则的真实 EC**（真/假一视同仁）→ EC 不再泄漏标签。
  - 只保留"规则确实生成了真产物"的底物（generation hit）；其余=generation miss，单独计数。
真产物若没被任何规则生成 → 该底物丢弃（计入 miss），因为部署时也排不出来。

输出 clean_candidates.parquet + 每模块 generation recall。
用法：python scripts/ssrf/ml_ranking_kernel/modular/ltr_clean_candidates.py [--n-rules 800] [--weak-cap 600] [--force]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402
from ltr_build import LTR_OUT  # noqa: E402
from ltr_candidates import load_module_rules  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-rules", type=int, default=800)
    ap.add_argument("--weak-cap", type=int, default=600)
    ap.add_argument("--max-cand", type=int, default=120)
    ap.add_argument("--require-ec", action="store_true", help="只用有 EC(EC≠0)的规则 → 中和 ec_only confound")
    ap.add_argument("--out", default="clean_candidates.parquet")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    rng = np.random.default_rng(0)

    rxn = pd.read_parquet(LTR_OUT / "reactions.parquet")
    true_prod = rxn.groupby(["module", "sb"])["pb"].agg(set).to_dict()
    rules = load_module_rules(rng, args.n_rules)
    if args.require_ec:
        rules = {m: [(s, e) for s, e in v if e != 0] for m, v in rules.items()}
        kio.log.info("EC-neutralized: 只保留有 EC 的规则")
    kio.log.info("module rules: %s", {m: len(v) for m, v in rules.items()})

    rows = []; genrec = {}
    for m in ["A", "B", "C", "D"]:
        cols = ["sb", "substrate_smiles", "substrate_scaffold", "is_clean", "year"]
        if "tier" in rxn.columns:
            cols.append("tier")
        sub = rxn[rxn["module"] == m].drop_duplicates("sb")[cols].copy()
        if "tier" not in sub.columns:
            sub["tier"] = np.where(sub["is_clean"], "gold", "weak")
        prio = sub[sub["tier"].isin(["gold", "silver"])]
        weak = sub[sub["tier"] == "weak"]
        if len(weak) > args.weak_cap:
            weak = weak.sample(n=args.weak_cap, random_state=0)
        sub = pd.concat([prio, weak]).drop_duplicates("sb")
        hit = tot = 0
        for si, r in enumerate(sub.itertuples(index=False)):
            s = r.substrate_smiles; sb = r.sb
            cmol = Chem.MolFromSmiles(s)
            if cmol is None:
                continue
            cmolh = Chem.AddHs(cmol)
            truth = true_prod.get((m, sb), set())
            cand = {}   # pb -> (smiles, ec)  规则生成产物 + 其生成规则 EC
            for smarts, ecc in rules[m]:
                for psmi in run_reactants(smarts, cmol, cmolh):
                    pik = kio.smiles_to_inchikey(psmi); pb = kio.inchikey_block1(pik)
                    if not pb or pb == sb:
                        continue
                    if pb not in cand:
                        cand[pb] = (psmi, int(ecc))
            tot += 1
            gen_truth = truth & set(cand.keys())
            if not gen_truth:
                continue                       # generation miss → 部署时也排不出，丢弃
            hit += 1
            # 候选 = 全部规则生成产物（真/假都用规则 EC，无注入）
            items = list(cand.items())
            # 控量：保真产物 + 采样 decoy
            decoy_items = [(pb, v) for pb, v in items if pb not in truth]
            if len(decoy_items) > args.max_cand:
                idx = rng.choice(len(decoy_items), args.max_cand, replace=False)
                decoy_items = [decoy_items[i] for i in idx]
            keep = [(pb, cand[pb]) for pb in gen_truth] + decoy_items
            for pb, (psm, ecc) in keep:
                rows.append({"module": m, "sb": sb, "substrate_smiles": s,
                             "substrate_scaffold": r.substrate_scaffold,
                             "pb": pb, "product_smiles": psm, "ec_class": int(ecc),
                             "y": 1 if pb in truth else 0,
                             "is_clean": bool(r.is_clean), "tier": getattr(r, "tier", "weak"), "year": r.year})
            if (si + 1) % 300 == 0:
                kio.log.info("  %s %d/%d subs (hit=%d) rows=%d", m, si + 1, len(sub), hit, len(rows))
        genrec[m] = (hit, tot)
        kio.log.info("module %s generation recall: %d/%d = %.2f", m, hit, tot, hit / max(tot, 1))

    cand = pd.DataFrame(rows)
    # 每底物至少 1 正 1 负
    grp = cand.groupby("sb")["y"].agg(["sum", "count"])
    keep_subs = set(grp[(grp["sum"] >= 1) & (grp["count"] - grp["sum"] >= 1)].index)
    cand = cand[cand["sb"].isin(keep_subs)].copy()
    kio.write_table(cand, LTR_OUT / args.out, force=args.force)
    kio.log.info("== CLEAN candidates: %d rows, %d substrates (median/sub=%.0f) ==",
                 len(cand), cand["sb"].nunique(), cand.groupby("sb").size().median())
    kio.log.info("test-eligible (gold+silver) substrates per module: %s",
                 cand[cand.tier.isin(["gold", "silver"])].drop_duplicates("sb").groupby("module").size().to_dict())
    kio.log.info("generation recall per module: %s", {k: f"{v[0]}/{v[1]}" for k, v in genrec.items()})
    # EC 一致性自检：真/假 的 EC==0 比例应接近（不再泄漏）
    kio.log.info("[anti-cheat check] EC==0 rate  pos=%.2f  decoy=%.2f (应接近)",
                 (cand[cand.y == 1].ec_class == 0).mean(), (cand[cand.y == 0].ec_class == 0).mean())


if __name__ == "__main__":
    main()
