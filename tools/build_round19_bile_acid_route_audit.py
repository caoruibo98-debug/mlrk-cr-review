from __future__ import annotations

import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "curation"
DOC = ROOT / "docs" / "reviews" / "bile_acid_route_round19.md"
POOL = Path(r"D:\CRB\FoodGut\positive_sample_db\data\processed\foodgut_modular_positive_candidate_pool_v2.csv")
LTR_OUT = ROOT / "outputs" / "modular" / "ltr"
RHEA_RAW = OUT / "rhea_round19_raw"
NCBI_RAW = OUT / "ncbi_round19_raw"

sys.path.insert(0, str(ROOT / "modular"))
sys.path.insert(0, str(ROOT / "src"))
import kio  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402
from ltr_candidates import load_module_rules  # noqa: E402


RDLogger.DisableLog("rdApp.*")

TARGETS = [
    {
        "round": 19,
        "source_id": "RHEA:16309",
        "priority_family": "bile_acid_deconjugation_and_lipid_context",
        "module": "D",
        "sb": "BHTRKEVKTKCXOH",
        "pb": "RUDATBOHQWOJDD",
        "substrate_name": "taurochenodeoxycholate",
        "product_name": "chenodeoxycholate",
        "co_product_name": "taurine",
        "expected_ec": "3.5.1.74",
        "expected_ec_class": 3,
        "expected_ec_source": "https://www.rhea-db.org/rhea/16309",
    },
    {
        "round": 19,
        "source_id": "RHEA:19353",
        "priority_family": "bile_acid_deconjugation_and_lipid_context",
        "module": "D",
        "sb": "RFDAIACWWDREDC",
        "pb": "BHQCQFFYRZLCQQ",
        "substrate_name": "glycocholate",
        "product_name": "cholate",
        "co_product_name": "glycine",
        "expected_ec": "3.5.1.24",
        "expected_ec_class": 3,
        "expected_ec_source": "https://enzyme.expasy.org/EC/3.5.1.24",
    },
]

MACRO2MOD = {
    "small_molecule_polyphenol": "A",
    "carbohydrate_glycan": "B",
    "protein_amino_acid": "C",
    "lipid_fat": "D",
}


def ec_class(ec: object) -> int:
    text = "" if pd.isna(ec) else str(ec)
    return int(text[0]) if text[:1].isdigit() else 0


def key_block1(smiles: str) -> str:
    inchikey = kio.smiles_to_inchikey(smiles)
    return kio.inchikey_block1(inchikey)


def read_rhea_summary() -> pd.DataFrame:
    rows = []
    for target in TARGETS:
        path = RHEA_RAW / f"{target['source_id'].replace(':', '_')}.json"
        if not path.exists():
            rows.append({**target, "rhea_status": "missing_raw", "rhea_equation": "", "balanced": None})
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        results = data.get("results", [])
        if not results:
            rows.append({**target, "rhea_status": "no_result", "rhea_equation": "", "balanced": None})
            continue
        rec = results[0]
        rows.append(
            {
                **target,
                "rhea_status": rec.get("status", ""),
                "rhea_equation": rec.get("equation", ""),
                "balanced": rec.get("balanced", None),
                "raw_file": str(path.relative_to(ROOT)),
            }
        )
    return pd.DataFrame(rows)


def read_pubmed_bsh_leads() -> pd.DataFrame:
    path = NCBI_RAW / "bile_salt_hydrolase_known_pmids_esummary.json"
    if not path.exists():
        return pd.DataFrame()
    data = json.loads(path.read_text(encoding="utf-8"))
    result = data.get("result", {})
    rows = []
    for uid in result.get("uids", []):
        rec = result.get(uid, {})
        rows.append(
            {
                "round": 19,
                "source_type": "PubMed",
                "pmid": uid,
                "title": rec.get("title", ""),
                "journal": rec.get("source", ""),
                "pubdate": rec.get("pubdate", ""),
                "training_use": "mechanism_or_context_only_until_exact_assay_extraction",
                "raw_file": str(path.relative_to(ROOT)),
            }
        )
    return pd.DataFrame(rows)


def load_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    return (
        pd.read_parquet(LTR_OUT / "reactions.parquet"),
        pd.read_parquet(LTR_OUT / "clean_candidates_full.parquet"),
        pd.read_parquet(LTR_OUT / "clean_candidates.parquet"),
    )


def joined(values: Iterable[object]) -> str:
    return ";".join(sorted({str(v) for v in values if pd.notna(v)}))


def route_rows(reactions: pd.DataFrame, clean_full: pd.DataFrame, clean: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target in TARGETS:
        module, sb, pb = target["module"], target["sb"], target["pb"]
        rxn_sub = reactions[(reactions["module"].eq(module)) & (reactions["sb"].eq(sb))]
        full_sub = clean_full[(clean_full["module"].eq(module)) & (clean_full["sb"].eq(sb))]
        clean_sub = clean[(clean["module"].eq(module)) & (clean["sb"].eq(sb))]
        rxn_pair = rxn_sub[rxn_sub["pb"].eq(pb)]
        full_pair = full_sub[full_sub["pb"].eq(pb)]
        clean_pair = clean_sub[clean_sub["pb"].eq(pb)]
        full_target_ec = joined(full_pair.get("ec_class", pd.Series(dtype=object)))
        clean_ec_classes = joined(clean_sub.get("ec_class", pd.Series(dtype=object)))
        target_loss_cause = "unknown"
        if len(full_pair) and not len(clean_pair) and full_target_ec in {"0", "0.0"} and "0" not in clean_ec_classes.split(";"):
            target_loss_cause = "target_reaches_clean_full_as_ec0_then_absent_from_ec_bearing_final_clean"
        rows.append(
            {
                **target,
                "pair_key": f"{sb}__{pb}",
                "reactions_pair_rows": len(rxn_pair),
                "reaction_tier": joined(rxn_pair.get("tier", pd.Series(dtype=object))),
                "reaction_is_clean": bool(rxn_pair["is_clean"].max()) if len(rxn_pair) else False,
                "reaction_ev_ec_class": joined(rxn_pair.get("ev_ec_class", pd.Series(dtype=object))),
                "reaction_ev_has_rule": joined(rxn_pair.get("ev_has_rule", pd.Series(dtype=object))),
                "clean_full_pair_rows": len(full_pair),
                "clean_full_target_ec_class": full_target_ec,
                "clean_full_positive_products_for_substrate": joined(full_sub.loc[full_sub["y"].eq(1), "pb"]),
                "clean_pair_rows": len(clean_pair),
                "clean_positive_products_for_substrate": joined(clean_sub.loc[clean_sub["y"].eq(1), "pb"]),
                "clean_ec_classes_for_substrate": clean_ec_classes,
                "target_loss_cause_round19": target_loss_cause,
                "route_repair_decision": "repair_rule_ec_provenance_before_retraining",
                "training_allowed_round19": False,
            }
        )
    return pd.DataFrame(rows)


def run_rule_set(rule_set: list[tuple[str, int]], substrate_smiles: str, target_pb: str) -> tuple[bool, int, str]:
    mol = Chem.MolFromSmiles(substrate_smiles)
    if mol is None:
        return False, 0, ""
    molh = Chem.AddHs(mol)
    generated: dict[str, set[int]] = {}
    for smarts, ecc in rule_set:
        for product_smiles in run_reactants(smarts, mol, molh):
            pb = key_block1(product_smiles)
            if not pb:
                continue
            generated.setdefault(pb, set()).add(int(ecc))
    return target_pb in generated, len(generated), joined(generated.get(target_pb, set()))


def sampled_rule_sets() -> dict[str, list[tuple[str, int]]]:
    rng = np.random.default_rng(0)
    sampled = load_module_rules(rng, 800)["D"]
    all_rules = load_module_rules(np.random.default_rng(0), 100000)["D"]
    return {
        "sampled_800_no_ec_filter": sampled,
        "sampled_800_require_ec": [(s, e) for s, e in sampled if e != 0],
        "all_pool_no_ec_filter": all_rules,
        "all_pool_require_ec": [(s, e) for s, e in all_rules if e != 0],
    }


def rule_variant_rows(reactions: pd.DataFrame) -> pd.DataFrame:
    rule_sets = sampled_rule_sets()
    rows = []
    for target in TARGETS:
        substrate_smiles = reactions.loc[
            (reactions["module"].eq(target["module"])) & (reactions["sb"].eq(target["sb"])),
            "substrate_smiles",
        ].iloc[0]
        for variant, rules in rule_sets.items():
            hit, generated_count, target_ec_classes = run_rule_set(rules, substrate_smiles, target["pb"])
            rows.append(
                {
                    **target,
                    "rule_variant": variant,
                    "rule_count": len(rules),
                    "target_generated": hit,
                    "generated_product_count": generated_count,
                    "target_generated_ec_classes": target_ec_classes,
                    "variant_gate": "passes_target_route" if hit else "drops_target",
                }
            )
    return pd.DataFrame(rows)


def pool_rule_frame() -> pd.DataFrame:
    wanted = [
        "macro_module",
        "reaction_rule_smarts",
        "enzyme_ec",
        "source_dataset",
        "source_origin_type",
        "training_use_recommendation",
        "reaction_category",
        "reaction_type",
        "reaction_rule_id",
    ]
    data = pd.read_csv(POOL, low_memory=False, usecols=lambda c: c in wanted)
    data["module"] = data["macro_module"].map(MACRO2MOD)
    data = data.dropna(subset=["module", "reaction_rule_smarts"])
    data = data[data["reaction_rule_smarts"].astype(str).str.contains(">>", regex=False)].copy()
    data["ecc"] = data["enzyme_ec"].map(ec_class)
    data["rule_hash"] = data["reaction_rule_smarts"].map(lambda s: hashlib.sha1(str(s).encode("utf-8")).hexdigest()[:12])
    return data


def target_generating_rule_rows(reactions: pd.DataFrame) -> pd.DataFrame:
    rules = pool_rule_frame()
    unique = rules[rules["module"].eq("D")].drop_duplicates("reaction_rule_smarts")
    out = []
    for target in TARGETS:
        substrate_smiles = reactions.loc[
            (reactions["module"].eq(target["module"])) & (reactions["sb"].eq(target["sb"])),
            "substrate_smiles",
        ].iloc[0]
        mol = Chem.MolFromSmiles(substrate_smiles)
        molh = Chem.AddHs(mol)
        for row in unique.itertuples(index=False):
            generated = set()
            for product_smiles in run_reactants(row.reaction_rule_smarts, mol, molh):
                generated.add(key_block1(product_smiles))
            if target["pb"] not in generated:
                continue
            same_rule = rules[rules["reaction_rule_smarts"].eq(row.reaction_rule_smarts)]
            out.append(
                {
                    **target,
                    "rule_hash": row.rule_hash,
                    "rule_ecc": int(row.ecc),
                    "source_datasets": joined(same_rule.get("source_dataset", [])),
                    "source_origin_types": joined(same_rule.get("source_origin_type", [])),
                    "training_recommendations": joined(same_rule.get("training_use_recommendation", [])),
                    "reaction_categories": joined(same_rule.get("reaction_category", [])),
                    "reaction_types": joined(same_rule.get("reaction_type", [])),
                    "reaction_rule_ids": joined(same_rule.get("reaction_rule_id", [])),
                    "pool_rows_for_rule": len(same_rule),
                    "smarts_preview": str(row.reaction_rule_smarts)[:180],
                    "ec_repair_needed": int(row.ecc) == 0,
                }
            )
    return pd.DataFrame(out)


def fix_recommendations(
    routes: pd.DataFrame, variants: pd.DataFrame, target_rules: pd.DataFrame, rhea: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    for target in TARGETS:
        mask = lambda df: df["source_id"].eq(target["source_id"])
        route = routes[mask(routes)].iloc[0]
        variant_sub = variants[mask(variants)]
        rule_sub = target_rules[mask(target_rules)] if not target_rules.empty else pd.DataFrame()
        rhea_sub = rhea[mask(rhea)]
        sampled_require_ec_hit = bool(
            variant_sub.loc[variant_sub["rule_variant"].eq("sampled_800_require_ec"), "target_generated"].any()
        )
        all_require_ec_hit = bool(
            variant_sub.loc[variant_sub["rule_variant"].eq("all_pool_require_ec"), "target_generated"].any()
        )
        all_no_ec_hit = bool(
            variant_sub.loc[variant_sub["rule_variant"].eq("all_pool_no_ec_filter"), "target_generated"].any()
        )
        only_ec0_generates = all_no_ec_hit and not all_require_ec_hit and (
            not rule_sub.empty and set(rule_sub["rule_ecc"].astype(int)) == {0}
        )
        if only_ec0_generates:
            production_blocker = "ec_provenance_gap_no_ec_bearing_rule_generates_target"
            fix_action = "add_source_backed_ec_to_existing_bsh_rule_or_overlay_then_rerun_require_ec_dryrun"
        elif all_require_ec_hit and not sampled_require_ec_hit:
            production_blocker = "sampled_require_ec_rule_selection_gap"
            fix_action = "deterministically_include_source_backed_bsh_rules_in_require_ec_generator"
        else:
            production_blocker = "target_route_gap_requires_manual_rule_trace"
            fix_action = "trace_rule_selection_and_generator_output_before_retraining"
        rows.append(
            {
                **target,
                "rhea_approved_exact": bool(len(rhea_sub) and rhea_sub.iloc[0]["rhea_status"] == "approved"),
                "present_in_clean_full": bool(route["clean_full_pair_rows"] > 0),
                "present_in_final_clean": bool(route["clean_pair_rows"] > 0),
                "sampled_require_ec_generates_target": sampled_require_ec_hit,
                "all_pool_require_ec_generates_target": all_require_ec_hit,
                "only_ec0_rules_generate_target": only_ec0_generates,
                "production_blocker_round19": production_blocker,
                "expected_ec": target["expected_ec"],
                "expected_ec_class": target["expected_ec_class"],
                "fix_action": fix_action,
                "do_not_do": "do_not_duplicate_import_positive_and_do_not_disable_ec_filter_globally",
                "success_criterion": "target_generated=True under sampled require-EC or targeted source-backed overlay, then pair appears in final clean dry run",
                "training_allowed_round19": False,
            }
        )
    return pd.DataFrame(rows)


def write_doc(
    routes: pd.DataFrame,
    variants: pd.DataFrame,
    target_rules: pd.DataFrame,
    rhea: pd.DataFrame,
    pubmed: pd.DataFrame,
    fixes: pd.DataFrame,
) -> None:
    DOC.parent.mkdir(parents=True, exist_ok=True)
    text = "# Round19 Bile Acid Route Audit\n\n"
    text += "## Working conclusion\n\n"
    text += (
        "The P0 bile-acid issue is not a request for more positive rows. "
        "Both target reactions are already source-backed and present before the final clean path. "
        "The final route drops them for two related reasons: one target has no EC-bearing generating rule, "
        "and the other has EC-bearing rules in the full pool but they are not deterministically included in the sampled require-EC generator.\n\n"
    )
    text += "## Route loss\n\n```text\n"
    text += routes[
        [
            "source_id",
            "pair_key",
            "clean_full_target_ec_class",
            "clean_pair_rows",
            "clean_ec_classes_for_substrate",
            "target_loss_cause_round19",
        ]
    ].to_string(index=False)
    text += "\n```\n\n"
    text += "## Rule variant dry run\n\n```text\n"
    text += variants[
        ["source_id", "rule_variant", "rule_count", "target_generated", "target_generated_ec_classes"]
    ].to_string(index=False)
    text += "\n```\n\n"
    text += "## Target-generating rule summary\n\n```text\n"
    if target_rules.empty:
        text += "no target-generating rules found\n"
    else:
        text += target_rules[
            ["source_id", "rule_hash", "rule_ecc", "source_datasets", "pool_rows_for_rule", "ec_repair_needed"]
        ].to_string(index=False)
        text += "\n"
    text += "```\n\n"
    text += "## External evidence\n\n"
    text += "- RHEA:16309: taurochenodeoxycholate + H2O = chenodeoxycholate + taurine; Rhea page reports EC 3.5.1.74.\n"
    text += "- RHEA:19353: glycocholate + H2O = cholate + glycine; Rhea/ENZYME reports EC 3.5.1.24.\n"
    text += "- PubMed BSH leads are context/mechanism evidence until exact assay tables are manually extracted.\n\n"
    text += "Rhea exact rows:\n\n```text\n"
    text += rhea[["source_id", "rhea_status", "rhea_equation", "balanced", "expected_ec"]].to_string(index=False)
    text += "\n```\n\n"
    text += "PubMed leads:\n\n```text\n"
    text += (
        pubmed[["pmid", "title", "training_use"]].to_string(index=False)
        if not pubmed.empty
        else "no PubMed lead rows"
    )
    text += "\n```\n\n"
    text += "## Round20 gate\n\n```text\n"
    text += fixes[
        [
            "source_id",
            "production_blocker_round19",
            "fix_action",
            "success_criterion",
            "training_allowed_round19",
        ]
    ].to_string(index=False)
    text += "\n```\n\n"
    text += "## Written artifacts\n\n"
    text += "- `data/curation/round19_bile_acid_route_loss.csv`\n"
    text += "- `data/curation/round19_bile_acid_rule_variant_dryrun.csv`\n"
    text += "- `data/curation/round19_bile_acid_target_generating_rules.csv`\n"
    text += "- `data/curation/round19_bile_acid_external_evidence.csv`\n"
    text += "- `data/curation/round19_bile_acid_pubmed_leads.csv`\n"
    text += "- `data/curation/round19_bile_acid_fix_recommendations.csv`\n"
    DOC.write_text(text, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    reactions, clean_full, clean = load_tables()
    rhea = read_rhea_summary()
    pubmed = read_pubmed_bsh_leads()
    routes = route_rows(reactions, clean_full, clean)
    variants = rule_variant_rows(reactions)
    target_rules = target_generating_rule_rows(reactions)
    fixes = fix_recommendations(routes, variants, target_rules, rhea)

    outputs = {
        "round19_bile_acid_route_loss.csv": routes,
        "round19_bile_acid_rule_variant_dryrun.csv": variants,
        "round19_bile_acid_target_generating_rules.csv": target_rules,
        "round19_bile_acid_external_evidence.csv": rhea,
        "round19_bile_acid_pubmed_leads.csv": pubmed,
        "round19_bile_acid_fix_recommendations.csv": fixes,
    }
    for name, df in outputs.items():
        path = OUT / name
        df.to_csv(path, index=False)
        print(f"wrote {path.relative_to(ROOT)} rows={len(df)}")
    write_doc(routes, variants, target_rules, rhea, pubmed, fixes)
    print(f"wrote {DOC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
