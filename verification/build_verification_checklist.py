"""
build_verification_checklist.py — 把"需要人工核验真伪"的样本汇成一张可填表

覆盖:
  - 29 条 curated 规则(curated_gut_smarts.csv):核 PMID↔化学 + SMARTS 有效性
  - gap7 的 18 条 NEEDS_REVIEW(curated_literature_gap7):核结构TBD + 单步 + PMID
  - gold_labels 多酚集(gpt_polyphenol_round1_raw.csv):SMILES空 + doi未核

输出 verification/sample_verification_checklist.csv，列含 issue_flag/what_to_verify/verify_result(留空给人填)。
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from rdkit.Chem import AllChem
from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import kio  # noqa: E402

OUT = kio.OUTPUTS / "verification"
COLS = ["source", "id", "module", "reaction_type", "substrate_name", "substrate_smiles",
        "product_name", "product_smiles", "claimed_transformation", "evidence_pmid", "doi",
        "confidence", "issue_flag", "what_to_verify", "verify_result", "reviewer_notes"]


def smarts_ok(s):
    try:
        r = AllChem.ReactionFromSmarts(str(s)); return r is not None and r.GetNumReactantTemplates() > 0
    except Exception:
        return False


def has_tbd(*vals):
    return any("tbd" in str(v).lower() or str(v).strip() in ("", "nan", "?") for v in vals)


def main():
    kio.setup_logging()
    rows = []

    # --- 29 curated 规则 ---
    p = kio.resolve_input("Dataset/curated_gut_smarts.csv")
    if p.exists():
        c = pd.read_csv(p)
        for r in c.itertuples(index=False):
            sm = getattr(r, "rule_smarts", "")
            bad = not smarts_ok(sm)
            rows.append({"source": "curated_rule_29", "id": getattr(r, "rule_id", ""),
                         "module": "", "reaction_type": getattr(r, "reaction_type", ""),
                         "substrate_name": "", "substrate_smiles": "",
                         "product_name": "", "product_smiles": "",
                         "claimed_transformation": getattr(r, "description", ""),
                         "evidence_pmid": getattr(r, "evidence_pmid", ""), "doi": "",
                         "confidence": getattr(r, "confidence", ""),
                         "issue_flag": "SMARTS_INVALID" if bad else "verify_pmid_chem",
                         "what_to_verify": "①PMID是否真报道此转化 ②SMARTS化学是否正确"
                         + ("（⚠️SMARTS当前无法解析，需修或弃）" if bad else ""),
                         "verify_result": "", "reviewer_notes": ""})

    # --- gap7 NEEDS_REVIEW ---
    p = kio.resolve_input("Dataset/curated_literature/2026_04_19_microbiome_reaction_gap_fill/curated_literature_gap7_import_ready_v1.csv")
    if p.exists():
        g = pd.read_csv(p)
        rev = g[g["curator_flag"].astype(str).str.upper().ne("OK")]
        for r in rev.itertuples(index=False):
            tbd = has_tbd(getattr(r, "substrate_smiles", ""), getattr(r, "product_smiles", ""))
            rows.append({"source": "gap7_needs_review", "id": getattr(r, "rule_id", ""),
                         "module": "", "reaction_type": getattr(r, "reaction_category", ""),
                         "substrate_name": getattr(r, "substrate_name", ""),
                         "substrate_smiles": getattr(r, "substrate_smiles", ""),
                         "product_name": getattr(r, "product_name", ""),
                         "product_smiles": getattr(r, "product_smiles", ""),
                         "claimed_transformation": getattr(r, "reaction_mechanism_note", ""),
                         "evidence_pmid": getattr(r, "primary_pmid", ""), "doi": "",
                         "confidence": getattr(r, "confidence_tier", ""),
                         "issue_flag": "SMILES_TBD" if tbd else "NEEDS_REVIEW",
                         "what_to_verify": "①补全/核验结构（TBD?）②确认单步 ③PMID是否匹配",
                         "verify_result": "", "reviewer_notes": ""})

    # --- gold_labels 多酚集（SMILES 多空、doi 未核）---
    p = kio.resolve_input("Dataset/gold_labels/gpt_polyphenol_round1_raw.csv")
    if p.exists():
        gl = pd.read_csv(p)
        sub_c = next((c for c in gl.columns if "substrate" in c.lower() and "name" in c.lower()), None)
        prod_c = next((c for c in gl.columns if "product" in c.lower() and "name" in c.lower()), None)
        pmid_c = next((c for c in gl.columns if "pmid" in c.lower()), None)
        doi_c = next((c for c in gl.columns if "doi" in c.lower()), None)
        for i, r in gl.head(60).iterrows():
            rows.append({"source": "gold_labels_polyphenol", "id": f"gpt_{i}",
                         "module": "", "reaction_type": "",
                         "substrate_name": r[sub_c] if sub_c else "", "substrate_smiles": "",
                         "product_name": r[prod_c] if prod_c else "", "product_smiles": "",
                         "claimed_transformation": "", "evidence_pmid": r[pmid_c] if pmid_c else "",
                         "doi": r[doi_c] if doi_c else "", "confidence": "",
                         "issue_flag": "STRUCTURE_MISSING+DOI_UNVERIFIED",
                         "what_to_verify": "①GPT抽取是否真实（核 PMID/DOI）②补结构",
                         "verify_result": "", "reviewer_notes": ""})

    df = pd.DataFrame(rows)[COLS]
    kio.write_table(df, OUT / "sample_verification_checklist.csv", force=True)
    kio.log.info("核验清单: %d 条", len(df))
    kio.log.info("按来源: %s", df.source.value_counts().to_dict())
    kio.log.info("按 issue: %s", df.issue_flag.value_counts().to_dict())


if __name__ == "__main__":
    main()
