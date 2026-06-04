from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
CUR = ROOT / "data" / "curation"
LTR = ROOT / "outputs" / "modular" / "ltr"
RHEA_RAW = CUR / "rhea_round25_raw"
NCBI_RAW = CUR / "ncbi_round25_raw"
ROUND = "round25"

RULE92_ID = "EnzymeMap:rule_id_92"
RULE92_SMARTS = "[#6:1]-[#6:2]-[#8:3]>>[#6:1].[#6:2]=[#8:3]"
RULE92_EC = "4.1.1.59/4.1.1.61/4.1.1.63"

sys.path.insert(0, str(ROOT / "src"))
from generate_candidates import run_reactants  # noqa: E402
import kio  # noqa: E402


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, low_memory=False)


def write_csv(df: pd.DataFrame, rel: str) -> Path:
    path = ROOT / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def load_json(path: Path) -> Any:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def first_pubchem_property(name: str) -> dict[str, Any]:
    path = NCBI_RAW / f"pubchem_{name}.json"
    data = load_json(path) or {}
    props = data.get("PropertyTable", {}).get("Properties", [])
    return props[0] if props else {}


def aromatic_carboxylate_scope(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return "invalid_smiles"
    patt = Chem.MolFromSmarts("[c][CX3](=O)[O;H1,-1]")
    if patt is not None and mol.HasSubstructMatch(patt):
        return "aromatic_carboxylate_like"
    patt2 = Chem.MolFromSmarts("[CX3](=O)[O;H1,-1]")
    if patt2 is not None and mol.HasSubstructMatch(patt2):
        return "carboxylate_like"
    return "other_positive_substrate"


def block1(smiles: Any) -> str:
    ik = kio.smiles_to_inchikey(smiles)
    return kio.inchikey_block1(ik)


def build_rule92_overgeneration_audit() -> tuple[pd.DataFrame, pd.DataFrame]:
    reactions = read_table(LTR / "reactions.parquet")
    clean_full = read_table(LTR / "clean_candidates_full.parquet")
    clean = read_table(LTR / "clean_candidates.parquet")

    pos_pairs = {
        (str(r.sb), str(r.pb))
        for r in reactions[["sb", "pb"]].dropna().itertuples(index=False)
    }
    clean_full_pairs = {
        (str(r.sb), str(r.pb), int(r.y))
        for r in clean_full[["sb", "pb", "y"]].dropna().itertuples(index=False)
    }
    clean_pairs = {
        (str(r.sb), str(r.pb), int(r.y))
        for r in clean[["sb", "pb", "y"]].dropna().itertuples(index=False)
    }
    clean_full_pair_y = {
        (str(r.sb), str(r.pb)): int(r.y)
        for r in clean_full[["sb", "pb", "y"]].dropna().itertuples(index=False)
    }
    clean_pair_y = {
        (str(r.sb), str(r.pb)): int(r.y)
        for r in clean[["sb", "pb", "y"]].dropna().itertuples(index=False)
    }

    subs = (
        reactions[
            ["module", "sb", "substrate_smiles", "substrate_scaffold", "tier"]
        ]
        .drop_duplicates(["module", "sb"])
        .dropna(subset=["substrate_smiles"])
        .reset_index(drop=True)
    )

    rows: list[dict[str, Any]] = []
    substrate_rows: list[dict[str, Any]] = []
    for r in subs.itertuples(index=False):
        smi = str(r.substrate_smiles)
        c_mol = Chem.MolFromSmiles(smi)
        c_mol_h = Chem.AddHs(c_mol) if c_mol is not None else None
        products = sorted(set(run_reactants(RULE92_SMARTS, c_mol, c_mol_h)))
        scope = aromatic_carboxylate_scope(smi)
        target_hits = 0
        clean_full_hits = 0
        clean_hits = 0
        unknown_hits = 0

        for psmi in products:
            pb = block1(psmi)
            if not pb or pb == str(r.sb):
                continue
            pair = (str(r.sb), pb)
            known_positive = pair in pos_pairs
            cf_y = clean_full_pair_y.get(pair)
            c_y = clean_pair_y.get(pair)
            in_clean_full = any((str(r.sb), pb, y) in clean_full_pairs for y in (0, 1))
            in_clean = any((str(r.sb), pb, y) in clean_pairs for y in (0, 1))
            if known_positive:
                target_hits += 1
                status = "known_positive_recovered"
                fn_status = "not_negative"
            elif in_clean_full and cf_y == 0:
                clean_full_hits += 1
                status = "unlabeled_current_decoy_overlap"
                fn_status = "possible_false_negative_do_not_harden"
            elif in_clean and c_y == 0:
                clean_hits += 1
                status = "unlabeled_final_decoy_overlap"
                fn_status = "possible_false_negative_do_not_harden"
            else:
                unknown_hits += 1
                status = "new_unlabeled_generation"
                fn_status = "unknown_not_negative"

            rows.append(
                {
                    "round": ROUND,
                    "rule_source": "EnzymeMap",
                    "rule_id": RULE92_ID,
                    "rule_smarts": RULE92_SMARTS,
                    "rule_ec_context": RULE92_EC,
                    "module": r.module,
                    "substrate_block1": r.sb,
                    "substrate_smiles": smi,
                    "substrate_scope": scope,
                    "substrate_tier": r.tier,
                    "generated_product_smiles": psmi,
                    "generated_product_block1": pb,
                    "known_positive_same_pair": known_positive,
                    "current_clean_candidates_full_same_pair": in_clean_full,
                    "current_clean_candidates_full_y": cf_y if cf_y is not None else "",
                    "current_clean_candidates_same_pair": in_clean,
                    "current_clean_candidates_y": c_y if c_y is not None else "",
                    "overgeneration_status": status,
                    "false_negative_screen_status": fn_status,
                    "training_allowed_round25": False,
                    "rule_promotion_allowed_round25": False,
                }
            )

        substrate_rows.append(
            {
                "round": ROUND,
                "rule_id": RULE92_ID,
                "module": r.module,
                "substrate_block1": r.sb,
                "substrate_scope": scope,
                "substrate_tier": r.tier,
                "generated_product_count": len(products),
                "known_positive_generated_count": target_hits,
                "current_clean_full_decoy_overlap_count": clean_full_hits,
                "current_clean_decoy_overlap_count": clean_hits,
                "unknown_generated_count": unknown_hits,
            }
        )

    audit = pd.DataFrame(rows)
    sub_summary = pd.DataFrame(substrate_rows)
    total_pairs = len(audit)
    known = int(audit["known_positive_same_pair"].sum()) if total_pairs else 0
    current_decoy_overlap = int(
        audit["false_negative_screen_status"]
        .eq("possible_false_negative_do_not_harden")
        .sum()
    ) if total_pairs else 0
    scope_counts = (
        sub_summary["substrate_scope"].value_counts(dropna=False).to_dict()
        if not sub_summary.empty
        else {}
    )
    product_scope = (
        audit.groupby("substrate_scope").size().to_dict() if not audit.empty else {}
    )
    decision = pd.DataFrame(
        [
            {
                "round": ROUND,
                "rule_id": RULE92_ID,
                "rule_source": "EnzymeMap",
                "rule_smarts": RULE92_SMARTS,
                "rule_ec_context": RULE92_EC,
                "screening_unit": "module_substrate_row",
                "module_substrate_rows_screened": len(sub_summary),
                "unique_substrate_blocks_screened": int(sub_summary["substrate_block1"].nunique()),
                "module_substrate_rows_generating_any_product": int(
                    (sub_summary["generated_product_count"] > 0).sum()
                ),
                "generated_pair_rows": total_pairs,
                "generated_pairs_known_positive": known,
                "generated_pairs_current_decoy_overlap": current_decoy_overlap,
                "known_positive_proxy_rate": known / total_pairs if total_pairs else 0,
                "substrate_scope_counts": json.dumps(scope_counts, sort_keys=True),
                "generated_product_scope_counts": json.dumps(
                    product_scope, sort_keys=True
                ),
                "rule92_decision": (
                    "do_not_promote_broad_template_without_family_specific_scope"
                ),
                "why": (
                    "The rule can recover the protocatechuate->catechol type, but it is a broad "
                    "decarboxylation-like template that emits many unlabeled products; unlabeled "
                    "products are not biological negatives."
                ),
                "minimum_promotion_gate": (
                    "Restrict to source-backed aromatic carboxylate subfamily, require exact "
                    "Rhea/EnzymeMap/PMID evidence, atom mapping, and a false-negative audit before "
                    "candidate generation overlay."
                ),
                "training_allowed_round25": False,
                "rule_promotion_allowed_round25": False,
            }
        ]
    )
    return audit, decision


def build_urolithin_exact_extraction_queue() -> pd.DataFrame:
    pubchem_a = first_pubchem_property("urolithin_A")
    pubchem_c = first_pubchem_property("urolithin_C")
    pubchem_m6 = first_pubchem_property("urolithin_M6")
    rhea_urolithin = load_json(RHEA_RAW / "rhea_urolithin_dehydroxylase.json") or {}
    rhea_hits = int(rhea_urolithin.get("count", 0) or 0)

    rows = [
        {
            "round": ROUND,
            "pmid": "39856097",
            "year": 2025,
            "title": "Diet-derived urolithin A is produced by a dehydroxylase encoded by human gut Enterocloster species",
            "journal": "Nature Communications",
            "organism_or_system": "Enterocloster species; E. bolteae and E. asparagiformis ucd operon",
            "suspected_exact_pair": "urolithin C -> urolithin A",
            "substrate_pubchem_cid": pubchem_c.get("CID", ""),
            "substrate_smiles_pubchem": pubchem_c.get("CanonicalSMILES")
            or pubchem_c.get("SMILES", ""),
            "substrate_inchikey_pubchem": pubchem_c.get("InChIKey", ""),
            "product_pubchem_cid": pubchem_a.get("CID", ""),
            "product_smiles_pubchem": pubchem_a.get("CanonicalSMILES")
            or pubchem_a.get("SMILES", ""),
            "product_inchikey_pubchem": pubchem_a.get("InChIKey", ""),
            "rhea_exact_hit_count_for_urolithin_dehydroxylase": rhea_hits,
            "evidence_status": "strong_paper_specific_operon_but_exact_table_extraction_pending",
            "exact_pair_table_status": "figure_and_source_data_need_extraction",
            "structure_status": "pubchem_structures_for_urolithin_A_and_C_present",
            "atom_mapping_required": True,
            "training_allowed_round25": False,
            "next_action": "Extract source data/figure table, verify exact names and positions, map atoms, then dry-run a restricted 9-dehydroxylation rule.",
        },
        {
            "round": ROUND,
            "pmid": "41298472",
            "year": 2025,
            "title": "The presence and induction of regioselective dehydroxylases dictate urolithin metabolism by Enterocloster species",
            "journal": "NPJ Biofilms Microbiomes",
            "organism_or_system": "Enterocloster species; regioselective urolithin dehydroxylases",
            "suspected_exact_pair": "urolithin M6/intermediates -> lower-hydroxylated urolithins",
            "substrate_pubchem_cid": pubchem_m6.get("CID", ""),
            "substrate_smiles_pubchem": pubchem_m6.get("CanonicalSMILES")
            or pubchem_m6.get("SMILES", ""),
            "substrate_inchikey_pubchem": pubchem_m6.get("InChIKey", ""),
            "product_pubchem_cid": "",
            "product_smiles_pubchem": "",
            "product_inchikey_pubchem": "",
            "rhea_exact_hit_count_for_urolithin_dehydroxylase": rhea_hits,
            "evidence_status": "paper_relevant_but_pair_direction_and_products_need_source_data_extraction",
            "exact_pair_table_status": "supplementary_table_and_source_data_required",
            "structure_status": "partial_pubchem_mapping_urolithin_M6_present",
            "atom_mapping_required": True,
            "training_allowed_round25": False,
            "next_action": "Extract all listed urolithins and LC-MS products; do not infer products from names alone.",
        },
        {
            "round": ROUND,
            "pmid": "41797252",
            "year": 2026,
            "title": "Urolithin 9-dehydroxylase from Enterocloster bolteae JCM 12243T catalyzing regiospecific dehydroxylation of urolithins",
            "journal": "Enzyme and Microbial Technology",
            "organism_or_system": "Enterocloster bolteae JCM 12243T recombinant/purified enzyme",
            "suspected_exact_pair": "five urolithins with 9-position dehydroxylation",
            "substrate_pubchem_cid": "",
            "substrate_smiles_pubchem": "",
            "substrate_inchikey_pubchem": "",
            "product_pubchem_cid": "",
            "product_smiles_pubchem": "",
            "product_inchikey_pubchem": "",
            "rhea_exact_hit_count_for_urolithin_dehydroxylase": rhea_hits,
            "evidence_status": "high_value_recent_enzyme_paper_but_exact_five_pairs_need_table_extraction",
            "exact_pair_table_status": "article_or_supplement_required",
            "structure_status": "not_mapped_in_round25",
            "atom_mapping_required": True,
            "training_allowed_round25": False,
            "next_action": "Retrieve open article/supplement, enumerate the five substrate-product pairs, then PubChem/ChEBI-map each molecule.",
        },
        {
            "round": ROUND,
            "pmid": "37494568",
            "year": 2023,
            "title": "NMR Spectroscopic Identification of Urolithin G, a Novel Trihydroxy Urolithin Produced by Human Intestinal Enterocloster Species",
            "journal": "Journal of Agricultural and Food Chemistry",
            "organism_or_system": "Human intestinal Enterocloster species",
            "suspected_exact_pair": "urolithin D -> urolithin G and related trihydroxy urolithins",
            "substrate_pubchem_cid": "",
            "substrate_smiles_pubchem": "",
            "substrate_inchikey_pubchem": "",
            "product_pubchem_cid": "",
            "product_smiles_pubchem": "",
            "product_inchikey_pubchem": "",
            "rhea_exact_hit_count_for_urolithin_dehydroxylase": rhea_hits,
            "evidence_status": "metabolite_identity_paper_but_not_yet_rule_level_training_data",
            "exact_pair_table_status": "table_1_and_structure_mapping_required",
            "structure_status": "urolithin_G_pubchem_name_lookup_failed_round25",
            "atom_mapping_required": True,
            "training_allowed_round25": False,
            "next_action": "Map Uro-D/Uro-G structures from paper or ChEBI/PubChem synonym search before sample import.",
        },
        {
            "round": ROUND,
            "pmid": "36840624",
            "year": 2023,
            "title": "Gut Bacteria Involved in Ellagic Acid Metabolism To Yield Human Urolithin Metabotypes Revealed",
            "journal": "Journal of Agricultural and Food Chemistry",
            "organism_or_system": "Gut bacteria involved in ellagic acid to urolithin metabotypes",
            "suspected_exact_pair": "ellagic acid/urolithin intermediates -> urolithin A/B/intermediates",
            "substrate_pubchem_cid": "",
            "substrate_smiles_pubchem": "",
            "substrate_inchikey_pubchem": "",
            "product_pubchem_cid": "",
            "product_smiles_pubchem": "",
            "product_inchikey_pubchem": "",
            "rhea_exact_hit_count_for_urolithin_dehydroxylase": rhea_hits,
            "evidence_status": "context_source_for_metabotypes_not_direct_rule_source_yet",
            "exact_pair_table_status": "requires_manual_table_extraction",
            "structure_status": "not_mapped_in_round25",
            "atom_mapping_required": True,
            "training_allowed_round25": False,
            "next_action": "Use as source-discovery and organism context; import exact pairs only after table-backed structure mapping.",
        },
    ]
    return pd.DataFrame(rows)


def build_negative_decoy_gate(rule_audit: pd.DataFrame) -> pd.DataFrame:
    if rule_audit.empty:
        possible_false = 0
        unknown = 0
    else:
        possible_false = int(
            rule_audit["false_negative_screen_status"]
            .eq("possible_false_negative_do_not_harden")
            .sum()
        )
        unknown = int(
            rule_audit["false_negative_screen_status"].eq("unknown_not_negative").sum()
        )
    return pd.DataFrame(
        [
            {
                "round": ROUND,
                "negative_strategy": "same_substrate_rule_generated_unlabeled_products",
                "can_be_y0_hard_negative": False,
                "can_be_ranking_decoy": True,
                "requires_source_absence_check": True,
                "requires_family_false_negative_screen": True,
                "observed_possible_false_negative_overlap_round25": possible_false,
                "observed_unknown_generation_round25": unknown,
                "reason": (
                    "Absence from current positive databases is not evidence that a gut microbe cannot "
                    "produce the metabolite. Treat these as decoys/unknowns unless an explicit negative "
                    "assay or incompatible biology source exists."
                ),
                "minimum_fields_for_decoy_csv": (
                    "substrate/product SMILES and InChIKey; generation rule; source family; "
                    "known-positive exclusion; source query manifest; split group; label_scope=decoy_not_biological_negative"
                ),
            },
            {
                "round": ROUND,
                "negative_strategy": "explicit_literature_or_database_no_activity_assay",
                "can_be_y0_hard_negative": True,
                "can_be_ranking_decoy": True,
                "requires_source_absence_check": False,
                "requires_family_false_negative_screen": True,
                "observed_possible_false_negative_overlap_round25": "",
                "observed_unknown_generation_round25": "",
                "reason": (
                    "A hard negative is defensible only when a source reports no conversion under "
                    "specified organism/enzyme/assay conditions."
                ),
                "minimum_fields_for_decoy_csv": (
                    "organism/enzyme; assay condition; substrate; absent product; detection method; "
                    "PMID/DOI/database id; condition scope; not-generalized-to-all-microbes flag"
                ),
            },
            {
                "round": ROUND,
                "negative_strategy": "cofactor_or_mass_balance_invalid_products",
                "can_be_y0_hard_negative": False,
                "can_be_ranking_decoy": False,
                "requires_source_absence_check": False,
                "requires_family_false_negative_screen": False,
                "observed_possible_false_negative_overlap_round25": "",
                "observed_unknown_generation_round25": "",
                "reason": (
                    "Invalid generated chemistry should be filtered, not used as a negative label; otherwise "
                    "the ranker learns RDKit/template artifacts."
                ),
                "minimum_fields_for_decoy_csv": "filter_reason; invalidation rule; no training label",
            },
        ]
    )


def build_external_query_manifest() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "round": ROUND,
                "source": "Rhea",
                "query": "RHEA:22416",
                "raw_path": str(RHEA_RAW / "rhea_22416.json"),
                "result_summary": "approved exact protocatechuate decarboxylase reaction",
                "url": "https://www.rhea-db.org/rhea/22416",
            },
            {
                "round": ROUND,
                "source": "Rhea",
                "query": "urolithin dehydroxylase",
                "raw_path": str(RHEA_RAW / "rhea_urolithin_dehydroxylase.json"),
                "result_summary": "zero exact Rhea hits in round25 query",
                "url": "https://www.rhea-db.org/",
            },
            {
                "round": ROUND,
                "source": "Rhea",
                "query": "dopamine 3-tyramine",
                "raw_path": str(RHEA_RAW / "rhea_dopamine_3_tyramine.json"),
                "result_summary": "approved exact dopamine + AH2 = 3-tyramine + A + H2O",
                "url": "https://www.rhea-db.org/",
            },
            {
                "round": ROUND,
                "source": "PubMed",
                "query": "urolithin dehydroxylase Enterocloster",
                "raw_path": str(
                    NCBI_RAW
                    / "pubmed_urolithin_dehydroxylase_enterocloster_esearch.json"
                ),
                "result_summary": "PMIDs 41797252, 41298472, 39856097, 37494568",
                "url": "https://pubmed.ncbi.nlm.nih.gov/",
            },
            {
                "round": ROUND,
                "source": "PubChem",
                "query": "urolithin A/C/M6 property lookups",
                "raw_path": str(NCBI_RAW),
                "result_summary": "PubChem structures found for urolithin A, C, and M6; urolithin G name lookup failed",
                "url": "https://pubchem.ncbi.nlm.nih.gov/",
            },
            {
                "round": ROUND,
                "source": "GitHub",
                "query": "rxn4chemistry/biocatalysis-model",
                "raw_path": "",
                "result_summary": "ECReact is drawn from Rhea, BRENDA, PathBank, MetaNetX and covers all 7 EC classes.",
                "url": "https://github.com/rxn4chemistry/biocatalysis-model",
            },
            {
                "round": ROUND,
                "source": "GitHub",
                "query": "hesther/enzymemap",
                "raw_path": "",
                "result_summary": "EnzymeMap processed_reactions.csv.gz is the current database file used for rule-level audit.",
                "url": "https://github.com/hesther/enzymemap",
            },
            {
                "round": ROUND,
                "source": "RetroRules",
                "query": "reaction rule generation docs",
                "raw_path": "",
                "result_summary": "RetroRules uses AAM reaction-center templates and mono-substrate decomposition; not all plausible food-gut reactions are precomputed.",
                "url": "https://retrorules.org/docs",
            },
        ]
    )


def write_review(
    rule_audit: pd.DataFrame,
    decision: pd.DataFrame,
    queue: pd.DataFrame,
    neg_gate: pd.DataFrame,
    manifest: pd.DataFrame,
) -> Path:
    dec = decision.iloc[0].to_dict()
    known_rate = float(dec["known_positive_proxy_rate"])
    lines = [
        "# Round25 rule92 overgeneration and urolithin extraction queue",
        "",
        "## Working conclusion",
        "",
        "Ray's current judgment is correct in the important sense: the production blocker is not the existence of only 22 samples. The 22-row file is a challenge panel. The blocking issue is earlier in the pipeline: source-backed reaction-family coverage and the label contract for generated negatives.",
        "",
        "## Rule 92 audit",
        "",
        f"- Rule: `{RULE92_SMARTS}` from EnzymeMap rule id 92 context ({RULE92_EC}).",
        f"- Module-substrate rows screened: {dec['module_substrate_rows_screened']}.",
        f"- Unique substrate blocks screened: {dec['unique_substrate_blocks_screened']}.",
        f"- Module-substrate rows with at least one generated product: {dec['module_substrate_rows_generating_any_product']}.",
        f"- Generated substrate-product rows: {dec['generated_pair_rows']}.",
        f"- Generated rows already known as positive: {dec['generated_pairs_known_positive']} ({known_rate:.4f} proxy rate; this is not precision because most rows are unlabeled).",
        f"- Generated rows overlapping current decoy space: {dec['generated_pairs_current_decoy_overlap']}.",
        "",
        "Interpretation: rule 92 can generate the protocatechuate to catechol product, but as a broad template it creates many unlabeled products. These unlabeled products cannot be hardened as biological negatives.",
        "",
        "## Urolithin evidence queue",
        "",
        f"- Urolithin extraction queue contains {len(queue)} high-value source rows; the targeted PubMed Enterocloster query directly returned four dehydroxylase-focused PMIDs.",
        "- Rhea query for `urolithin dehydroxylase` returned zero hits in this round.",
        "- PubChem mappings are present for urolithin A, C, and M6; urolithin G needs a synonym/structure-specific lookup.",
        "- No urolithin row is training-allowed in round25; every row needs exact table extraction and atom mapping first.",
        "",
        "## Negative sample policy",
        "",
        "- Missing database evidence is not a negative label.",
        "- Rule-generated unknowns can be ranking decoys only if the CSV says `label_scope=decoy_not_biological_negative`.",
        "- Hard negatives require explicit no-conversion assay evidence under a named organism/enzyme/condition.",
        "",
        "## Files",
        "",
        "- `data/curation/round25_rule92_overgeneration_audit.csv`",
        "- `data/curation/round25_rule92_candidate_decision.csv`",
        "- `data/curation/round25_urolithin_exact_extraction_queue.csv`",
        "- `data/curation/round25_negative_decoy_gate.csv`",
        "- `data/curation/round25_external_query_manifest.csv`",
        "",
        "## Sources recorded",
        "",
    ]
    for r in manifest.to_dict("records"):
        lines.append(f"- {r['source']}: {r['url']} ({r['result_summary']})")
    path = ROOT / "docs" / "reviews" / "round25_rule92_overgeneration_and_urolithin_queue.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> None:
    rule_audit, decision = build_rule92_overgeneration_audit()
    queue = build_urolithin_exact_extraction_queue()
    neg_gate = build_negative_decoy_gate(rule_audit)
    manifest = build_external_query_manifest()

    write_csv(rule_audit, "data/curation/round25_rule92_overgeneration_audit.csv")
    write_csv(decision, "data/curation/round25_rule92_candidate_decision.csv")
    write_csv(queue, "data/curation/round25_urolithin_exact_extraction_queue.csv")
    write_csv(neg_gate, "data/curation/round25_negative_decoy_gate.csv")
    write_csv(manifest, "data/curation/round25_external_query_manifest.csv")
    review = write_review(rule_audit, decision, queue, neg_gate, manifest)
    print(f"rule_audit_rows={len(rule_audit)}")
    print(f"decision_rows={len(decision)}")
    print(f"urolithin_queue_rows={len(queue)}")
    print(f"negative_gate_rows={len(neg_gate)}")
    print(f"review={review}")


if __name__ == "__main__":
    main()
