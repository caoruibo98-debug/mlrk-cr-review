#!/usr/bin/env python
"""Round17 reaction-family coverage audit.

Round16 showed that some source-backed positives are already present locally
but are lost before final clean candidates. Round17 moves up one level:
it quantifies reaction-family coverage across the raw pool, reactions,
clean_candidates_full, and clean_candidates, then attaches external Rhea/PubMed
evidence leads and a conservative negative-sampling policy.

This script is audit-only. It does not modify labels, rules, models, or
candidate-generation code.
"""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
CURATION = REPO / "data" / "curation"
DOCS = REPO / "docs" / "reviews"

RAW_POOL = Path(r"D:\CRB\FoodGut\positive_sample_db\data\processed\foodgut_modular_positive_candidate_pool_v2.csv")
REACTIONS = REPO / "outputs" / "modular" / "ltr" / "reactions.parquet"
CLEAN_FULL = REPO / "outputs" / "modular" / "ltr" / "clean_candidates_full.parquet"
CLEAN = REPO / "outputs" / "modular" / "ltr" / "clean_candidates.parquet"
METRICS = REPO / "outputs" / "modular" / "ltr" / "clean2_metrics.csv"
ROUND16_LOCAL = CURATION / "round16_clean_path_trace.csv"
ROUND16_DECISIONS = CURATION / "round16_training_decisions.csv"
RHEA17_RAW = CURATION / "rhea_round17_raw"
NCBI17_RAW = CURATION / "ncbi_round17_raw"

OUT_COVERAGE = CURATION / "round17_reaction_family_coverage_matrix.csv"
OUT_RHEA = CURATION / "round17_rhea_query_triage.csv"
OUT_PUBMED = CURATION / "round17_pubmed_query_triage.csv"
OUT_POSITIVE = CURATION / "round17_positive_expansion_queue.csv"
OUT_NEGATIVE = CURATION / "round17_negative_sampling_logic.csv"
OUT_GITHUB = CURATION / "round17_github_external_coverage_contract.csv"
OUT_DECISIONS = CURATION / "round17_training_decisions.csv"
OUT_DOC = DOCS / "reaction_family_coverage_round17.md"

MACRO2MOD = {
    "small_molecule_polyphenol": "A",
    "carbohydrate_glycan": "B",
    "protein_amino_acid": "C",
    "lipid_fat": "D",
}

LIT_SOURCES = {
    "benchmark_ssrf_rclss_2026-05-22",
    "curated_literature_gap7_import_ready_v1",
    "curated_literature_non7_extension_candidates_v1",
}

FAMILY_ORDER = [
    "bile_acid_deconjugation_and_lipid_context",
    "bile_acid_secondary_transformations",
    "glucuronide_sulfate_deconjugation",
    "glycoside_and_hmo_hydrolysis",
    "hydroxycinnamate_reduction_and_hydrolysis",
    "polyphenol_ring_fission",
    "urolithin_dehydroxylation",
    "isoflavone_reductive_and_glycoside_metabolism",
    "lignan_redox_and_deglycosylation",
    "prenylflavonoid_o_demethylation",
    "tryptophan_indole_metabolism",
    "tma_choline_carnitine_metabolism",
    "azoreductase_nitroreductase_xenobiotic",
    "protein_peptide_amino_acid_generic",
    "fatty_lipid_generic",
    "carbohydrate_generic",
    "polyphenol_generic",
    "other_or_unclassified",
]

FAMILY_MODULE_HINT = {
    "bile_acid_deconjugation_and_lipid_context": "D",
    "bile_acid_secondary_transformations": "D",
    "glucuronide_sulfate_deconjugation": "A",
    "glycoside_and_hmo_hydrolysis": "B",
    "hydroxycinnamate_reduction_and_hydrolysis": "A",
    "polyphenol_ring_fission": "A",
    "urolithin_dehydroxylation": "A",
    "isoflavone_reductive_and_glycoside_metabolism": "B",
    "lignan_redox_and_deglycosylation": "A",
    "prenylflavonoid_o_demethylation": "A",
    "tryptophan_indole_metabolism": "C",
    "tma_choline_carnitine_metabolism": "D",
    "azoreductase_nitroreductase_xenobiotic": "A",
    "protein_peptide_amino_acid_generic": "C",
    "fatty_lipid_generic": "D",
    "carbohydrate_generic": "B",
    "polyphenol_generic": "A",
    "other_or_unclassified": "",
}

GENERIC_FAMILIES = {
    "protein_peptide_amino_acid_generic",
    "fatty_lipid_generic",
    "carbohydrate_generic",
    "polyphenol_generic",
    "other_or_unclassified",
}

QUERY_FAMILIES = {
    "beta_glucuronidase": "glucuronide_sulfate_deconjugation",
    "glucuronide_hydrolase": "glucuronide_sulfate_deconjugation",
    "arylsulfatase": "glucuronide_sulfate_deconjugation",
    "cholate_7_alpha_dehydroxylase": "bile_acid_secondary_transformations",
    "bile_acid_7_alpha_dehydroxylation": "bile_acid_secondary_transformations",
    "tryptophan_indole_lyase": "tryptophan_indole_metabolism",
    "indolelactate_dehydrogenase": "tryptophan_indole_metabolism",
    "choline_trimethylamine_lyase": "tma_choline_carnitine_metabolism",
    "carnitine_trimethylamine": "tma_choline_carnitine_metabolism",
    "ferulate_reductase": "hydroxycinnamate_reduction_and_hydrolysis",
    "ferulic_acid_demethylation": "hydroxycinnamate_reduction_and_hydrolysis",
    "nitroreductase": "azoreductase_nitroreductase_xenobiotic",
    "azoreductase": "azoreductase_nitroreductase_xenobiotic",
}

PUBMED_QUERY_FAMILIES = {
    "gut_microbiota_beta_glucuronidase_polyphenol_glucuronide_deconjugation": "glucuronide_sulfate_deconjugation",
    "gut_microbiota_sulfatase_polyphenol_sulfate_deconjugation": "glucuronide_sulfate_deconjugation",
    "gut_microbiota_bile_acid_7_alpha_dehydroxylation_deoxycholic_acid": "bile_acid_secondary_transformations",
    "gut_microbiota_tryptophan_indole_indolepropionic_acid": "tryptophan_indole_metabolism",
    "gut_microbiota_choline_trimethylamine_lyase": "tma_choline_carnitine_metabolism",
    "gut_microbiota_carnitine_trimethylamine": "tma_choline_carnitine_metabolism",
    "gut_microbiota_ferulic_acid_dihydroferulic_acid_reduction": "hydroxycinnamate_reduction_and_hydrolysis",
    "gut_microbiota_azoreductase_nitroreductase_xenobiotic_metabolism": "azoreductase_nitroreductase_xenobiotic",
}


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def norm(value: object) -> str:
    return re.sub(r"\s+", " ", "" if pd.isna(value) else str(value)).strip()


def block1(value: object) -> str:
    text = norm(value)
    return text.split("-", 1)[0][:14] if text else ""


def pair_key(sb: object, pb: object) -> str:
    return f"{sb}__{pb}"


def text_blob(row: pd.Series) -> str:
    fields = [
        "substrate_name",
        "product_name",
        "reaction_type",
        "reaction_category",
        "reaction_step_class",
        "macro_module",
        "submodule",
        "module_rule",
        "enzyme_name",
        "enzyme_ec",
        "notes",
        "supporting_sentence",
        "database_name",
        "source_dataset",
    ]
    return " | ".join(norm(row.get(f, "")) for f in fields).lower()


def classify_family_from_text(blob: str, module: str = "") -> str:
    if any(k in blob for k in ["taurocholate", "taurochenodeoxycholate", "glycocholate", "glycochenodeoxycholate", "bile salt hydrolase", "bsh", "deconjugat"]):
        return "bile_acid_deconjugation_and_lipid_context"
    if any(k in blob for k in ["deoxycholic", "deoxycholate", "lithocholic", "lithocholate", "7 alpha", "7-alpha", "secondary bile", "hydroxysteroid dehydrogenase"]):
        return "bile_acid_secondary_transformations"
    if any(k in blob for k in ["glucuronide", "glucuronoside", "glucuronidase", "sulfate", "sulfatase", "aryl sulfate"]):
        return "glucuronide_sulfate_deconjugation"
    if any(k in blob for k in ["fucosyllactose", "fucosidase", "hmo", "human milk oligosaccharide", "glycoside", "glucoside", "rhamnoside", "disaccharide", "raffinose", "starch and sucrose"]):
        return "glycoside_and_hmo_hydrolysis"
    if any(k in blob for k in ["chlorogen", "caffeic", "caffeate", "ferulic", "ferulate", "hydroxycinnamate", "dihydrocaffeic", "dihydroferulic"]):
        return "hydroxycinnamate_reduction_and_hydrolysis"
    if any(k in blob for k in ["quercetin", "catechin", "epicatechin", "valerolactone", "ring fission", "flavonol", "flavan", "dihydroxyphenylacetic"]):
        return "polyphenol_ring_fission"
    if any(k in blob for k in ["urolithin", "ellagic", "ellagitannin"]):
        return "urolithin_dehydroxylation"
    if any(k in blob for k in ["daidzein", "genistein", "equol", "isoflavone", "dihydrodaidzein", "o-desmethylangolensin"]):
        return "isoflavone_reductive_and_glycoside_metabolism"
    if any(k in blob for k in ["pinoresinol", "lariciresinol", "enterodiol", "enterolactone", "lignan"]):
        return "lignan_redox_and_deglycosylation"
    if any(k in blob for k in ["isoxanthohumol", "prenylnaringenin", "prenylflavonoid", "o-demethylation", "demethylation"]):
        return "prenylflavonoid_o_demethylation"
    if any(k in blob for k in ["tryptophan", "tryptamine", "indole", "indolelactate", "indole-3", "ipa", "indolepropionic"]):
        return "tryptophan_indole_metabolism"
    if any(k in blob for k in ["trimethylamine", "tma", "choline", "carnitine", "cutc", "cnta"]):
        return "tma_choline_carnitine_metabolism"
    if any(k in blob for k in ["azoreductase", "nitroreductase", "azo dye", "nitrobenzene"]):
        return "azoreductase_nitroreductase_xenobiotic"
    if module == "C":
        return "protein_peptide_amino_acid_generic"
    if module == "D":
        return "fatty_lipid_generic"
    if module == "B":
        return "carbohydrate_generic"
    if module == "A":
        return "polyphenol_generic"
    return "other_or_unclassified"


def assign_raw_family(row: pd.Series) -> str:
    module = MACRO2MOD.get(norm(row.get("macro_module", "")), "")
    return classify_family_from_text(text_blob(row), module)


def tier_for_raw(row: pd.Series) -> str:
    is_gold = norm(row.get("source_origin_type", "")).eq("manual_literature_curated") if hasattr(norm(row.get("source_origin_type", "")), "eq") else False
    is_gold = norm(row.get("source_origin_type", "")) == "manual_literature_curated" or norm(row.get("source_dataset", "")) in LIT_SOURCES
    if is_gold:
        return "gold"
    rec = norm(row.get("training_use_recommendation", ""))
    if rec in {"database_positive_training_with_leakage_guard", "review_step_scope_before_strict_training"}:
        return "silver"
    return "weak"


def family_priority_from_values(values: list[str]) -> str:
    counts = Counter(v for v in values if v)
    if not counts:
        return "other_or_unclassified"
    for fam in FAMILY_ORDER:
        if fam in counts:
            return fam
    return counts.most_common(1)[0][0]


def build_raw_pair_map() -> pd.DataFrame:
    raw = pd.read_csv(RAW_POOL, low_memory=False)
    raw["sb"] = raw["substrate_inchikey"].map(block1)
    raw["pb"] = raw["product_inchikey"].map(block1)
    raw = raw[raw["sb"].ne("") & raw["pb"].ne("") & raw["sb"].ne(raw["pb"])].copy()
    raw["module"] = raw["macro_module"].map(MACRO2MOD)
    raw["family"] = raw.apply(assign_raw_family, axis=1)
    raw["tier"] = raw.apply(tier_for_raw, axis=1)
    raw["pair_key"] = raw.apply(lambda r: pair_key(r["sb"], r["pb"]), axis=1)

    rows = []
    for pk, g in raw.groupby("pair_key"):
        families = g["family"].tolist()
        tiers = set(g["tier"].tolist())
        rows.append(
            {
                "pair_key": pk,
                "sb": g["sb"].iloc[0],
                "pb": g["pb"].iloc[0],
                "family": family_priority_from_values(families),
                "all_family_votes": ";".join(f"{k}:{v}" for k, v in Counter(families).most_common()),
                "module_votes": ";".join(f"{k}:{v}" for k, v in Counter(g["module"].dropna().tolist()).most_common()),
                "raw_rows": len(g),
                "raw_tier": "gold" if "gold" in tiers else ("silver" if "silver" in tiers else "weak"),
                "source_datasets": ";".join(sorted(set(g["source_dataset"].dropna().astype(str)))),
                "example_substrate_name": norm(g["substrate_name"].dropna().iloc[0]) if g["substrate_name"].notna().any() else "",
                "example_product_name": norm(g["product_name"].dropna().iloc[0]) if g["product_name"].notna().any() else "",
            }
        )
    return pd.DataFrame(rows)


def add_pair_key(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["pair_key"] = out["sb"].astype(str) + "__" + out["pb"].astype(str)
    return out


def merge_family(df: pd.DataFrame, raw_map: pd.DataFrame, positive_only: bool = False) -> pd.DataFrame:
    if positive_only and "y" in df.columns:
        df = df[df["y"].eq(1)].copy()
    d = add_pair_key(df)
    d = d.merge(raw_map[["pair_key", "family", "raw_tier"]], on="pair_key", how="left")
    d["family"] = d.apply(
        lambda r: r["family"]
        if isinstance(r.get("family"), str) and r["family"]
        else classify_family_from_text("", str(r.get("module", ""))),
        axis=1,
    )
    return d


def summarize_by_family(df: pd.DataFrame, key: str) -> pd.DataFrame:
    rows = []
    for fam, g in df.groupby("family"):
        rows.append(
            {
                "family": fam,
                f"{key}_pairs": g["pair_key"].nunique(),
                f"{key}_substrates": g["sb"].nunique() if "sb" in g.columns else "",
                f"{key}_gold_pairs": g[g.get("tier", g.get("raw_tier", "")).eq("gold")]["pair_key"].nunique()
                if ("tier" in g.columns or "raw_tier" in g.columns)
                else "",
                f"{key}_silver_pairs": g[g.get("tier", g.get("raw_tier", "")).eq("silver")]["pair_key"].nunique()
                if ("tier" in g.columns or "raw_tier" in g.columns)
                else "",
                f"{key}_weak_pairs": g[g.get("tier", g.get("raw_tier", "")).eq("weak")]["pair_key"].nunique()
                if ("tier" in g.columns or "raw_tier" in g.columns)
                else "",
            }
        )
    return pd.DataFrame(rows)


def parse_rhea17() -> pd.DataFrame:
    rows = []
    for path in sorted(RHEA17_RAW.glob("*.json")):
        q = path.stem
        raw = read_json(path)
        records = raw.get("results", [])
        family = QUERY_FAMILIES.get(q, "other_or_unclassified")
        if not records:
            rows.append(
                {
                    "round": 17,
                    "query_slug": q,
                    "priority_family": family,
                    "rhea_id": "",
                    "equation": "",
                    "status": "no_rhea_result",
                    "balanced": "",
                    "candidate_class": "no_direct_rhea_hit",
                    "training_use": "not_a_sample",
                    "raw_file": str(path.relative_to(REPO)),
                }
            )
            continue
        for rec in records:
            equation = norm(rec.get("equation", ""))
            generic = equation.lower().startswith("an ") or " a " in equation.lower() or " an " in equation.lower()
            if rec.get("status") == "approved" and rec.get("balanced") is True and generic:
                candidate_class = "generic_rule_template_context"
                training_use = "rule_family_context_only_until exact compounds are mapped"
            elif rec.get("status") == "approved" and rec.get("balanced") is True:
                candidate_class = "approved_exact_or_specific_reaction_lead"
                training_use = "candidate_after ChEBI/PubChem mapping, local de-dup, and generator route check"
            else:
                candidate_class = "weak_context"
                training_use = "not_a_sample"
            rows.append(
                {
                    "round": 17,
                    "query_slug": q,
                    "priority_family": family,
                    "rhea_id": f"RHEA:{rec.get('id')}",
                    "equation": equation,
                    "status": rec.get("status", ""),
                    "balanced": rec.get("balanced", ""),
                    "candidate_class": candidate_class,
                    "training_use": training_use,
                    "raw_file": str(path.relative_to(REPO)),
                }
            )
    return pd.DataFrame(rows)


def parse_pubmed17() -> pd.DataFrame:
    summary = read_json(NCBI17_RAW / "round17_pubmed_esummary.json")
    result = summary.get("result", {})
    title_by_pmid = {pmid: result.get(pmid, {}).get("title", "") for pmid in result.get("uids", [])}
    journal_by_pmid = {
        pmid: result.get(pmid, {}).get("fulljournalname", result.get(pmid, {}).get("source", ""))
        for pmid in result.get("uids", [])
    }
    rows = []
    for path in sorted(NCBI17_RAW.glob("*_esearch.json")):
        q = re.sub(r"_esearch$", "", path.stem)
        raw = read_json(path)
        ids = raw.get("esearchresult", {}).get("idlist", [])
        family = PUBMED_QUERY_FAMILIES.get(q, "other_or_unclassified")
        if not ids:
            rows.append(
                {
                    "round": 17,
                    "query_slug": q,
                    "priority_family": family,
                    "pmid": "",
                    "title": "",
                    "journal": "",
                    "candidate_class": "search_miss",
                    "training_use": "not_a_sample",
                    "raw_file": str(path.relative_to(REPO)),
                }
            )
            continue
        for pmid in ids:
            title = title_by_pmid.get(pmid, "")
            low = title.lower()
            if any(k in low for k in ["lyase", "dehydrogenase", "azoreductase activity", "hydrolase", "characterization", "metabolism of oxo-bile"]):
                candidate_class = "enzyme_or_mechanism_lead"
                training_use = "candidate_after exact substrate/product extraction and assay/source verification"
            elif any(k in low for k in ["review", "axis", "biomarker", "levels", "supplementation", "prospective", "multi-omics", "dysbiosis"]):
                candidate_class = "context_only"
                training_use = "not_a_sample"
            else:
                candidate_class = "literature_lead_for_manual_exact_pair_review"
                training_use = "candidate_only_after exact substrate/product/evidence extraction"
            rows.append(
                {
                    "round": 17,
                    "query_slug": q,
                    "priority_family": family,
                    "pmid": pmid,
                    "title": title,
                    "journal": journal_by_pmid.get(pmid, ""),
                    "candidate_class": candidate_class,
                    "training_use": training_use,
                    "raw_file": str(path.relative_to(REPO)),
                }
            )
    return pd.DataFrame(rows)


def build_coverage_matrix(raw_map: pd.DataFrame, rhea: pd.DataFrame, pubmed: pd.DataFrame) -> pd.DataFrame:
    rxn = merge_family(pd.read_parquet(REACTIONS), raw_map)
    clean_full = merge_family(pd.read_parquet(CLEAN_FULL), raw_map, positive_only=True)
    clean = merge_family(pd.read_parquet(CLEAN), raw_map, positive_only=True)

    raw_sum = summarize_by_family(raw_map.rename(columns={"raw_tier": "tier"}), "raw")
    rxn_sum = summarize_by_family(rxn, "reactions")
    cf_sum = summarize_by_family(clean_full, "clean_full")
    clean_sum = summarize_by_family(clean, "clean")

    frames = [raw_sum, rxn_sum, cf_sum, clean_sum]
    out = pd.DataFrame({"family": FAMILY_ORDER})
    for frame in frames:
        out = out.merge(frame, on="family", how="left")
    for col in out.columns:
        if col != "family":
            out[col] = out[col].fillna(0).astype(int)

    rhea_counts = rhea.groupby(["priority_family", "candidate_class"]).size().unstack(fill_value=0)
    pubmed_counts = pubmed.groupby(["priority_family", "candidate_class"]).size().unstack(fill_value=0)
    route = pd.read_csv(ROUND16_LOCAL) if ROUND16_LOCAL.exists() else pd.DataFrame()
    route_counts = route.groupby("priority_family")["pair_key"].count().to_dict() if not route.empty else {}
    route_blockers = route[~route["loss_stage_round16"].eq("present_in_final_clean")].groupby("priority_family")["pair_key"].count().to_dict() if not route.empty else {}

    rows = []
    for r in out.itertuples(index=False):
        data = r._asdict()
        fam = data["family"]
        data["module_hint"] = FAMILY_MODULE_HINT.get(fam, "")
        data["lost_reactions_to_clean_full_pairs"] = max(data["reactions_pairs"] - data["clean_full_pairs"], 0)
        data["lost_clean_full_to_clean_pairs"] = max(data["clean_full_pairs"] - data["clean_pairs"], 0)
        data["round16_route_review_pairs"] = int(route_counts.get(fam, 0))
        data["round16_route_blocker_pairs"] = int(route_blockers.get(fam, 0))
        data["rhea17_specific_or_exact_leads"] = int(rhea_counts.get("approved_exact_or_specific_reaction_lead", pd.Series()).get(fam, 0)) if not rhea_counts.empty else 0
        data["rhea17_generic_template_context"] = int(rhea_counts.get("generic_rule_template_context", pd.Series()).get(fam, 0)) if not rhea_counts.empty else 0
        data["rhea17_no_direct_hit"] = int(rhea_counts.get("no_direct_rhea_hit", pd.Series()).get(fam, 0)) if not rhea_counts.empty else 0
        data["pubmed17_enzyme_or_mechanism_leads"] = int(pubmed_counts.get("enzyme_or_mechanism_lead", pd.Series()).get(fam, 0)) if not pubmed_counts.empty else 0
        data["pubmed17_manual_exact_pair_leads"] = int(pubmed_counts.get("literature_lead_for_manual_exact_pair_review", pd.Series()).get(fam, 0)) if not pubmed_counts.empty else 0
        data["pubmed17_context_only"] = int(pubmed_counts.get("context_only", pd.Series()).get(fam, 0)) if not pubmed_counts.empty else 0
        data["priority_score"] = priority_score(data)
        data["priority_level"] = "P0" if data["priority_score"] >= 7 else ("P1" if data["priority_score"] >= 4 else "P2")
        data["dominant_gap_type"] = dominant_gap(data)
        data["next_action"] = next_action(data)
        rows.append(data)
    return pd.DataFrame(rows)


def priority_score(data: dict[str, object]) -> int:
    if data["family"] in GENERIC_FAMILIES:
        return 0
    score = 0
    if int(data["round16_route_blocker_pairs"]) > 0:
        score += 4
    if int(data["clean_pairs"]) < 10 and (int(data["rhea17_specific_or_exact_leads"]) + int(data["pubmed17_enzyme_or_mechanism_leads"]) > 0):
        score += 3
    if int(data["clean_gold_pairs"]) < 5 and int(data["raw_gold_pairs"]) > 0:
        score += 2
    if int(data["lost_reactions_to_clean_full_pairs"]) > 20:
        score += 2
    if int(data["pubmed17_context_only"]) > 10 and int(data["clean_pairs"]) < 20:
        score += 1
    return score


def dominant_gap(data: dict[str, object]) -> str:
    if data["family"] in GENERIC_FAMILIES:
        return "generic_bucket_needs_subfamily_split"
    if int(data["round16_route_blocker_pairs"]) > 0:
        return "route_repair_before_sample_import"
    if int(data["clean_pairs"]) < 10 and (int(data["rhea17_specific_or_exact_leads"]) + int(data["pubmed17_enzyme_or_mechanism_leads"]) > 0):
        return "source_backed_family_undercoverage"
    if int(data["clean_gold_pairs"]) < 5 and int(data["raw_gold_pairs"]) > 0:
        return "gold_test_visibility_gap"
    if int(data["lost_reactions_to_clean_full_pairs"]) > 20:
        return "generator_recall_or_rule_pool_gap"
    if int(data["pubmed17_context_only"]) > 10:
        return "literature_context_needs_exact_pair_extraction"
    return "lower_priority_or_generic_family"


def next_action(data: dict[str, object]) -> str:
    gap = dominant_gap(data)
    if gap == "generic_bucket_needs_subfamily_split":
        return "split this broad bucket into named biochemical subfamilies before using it for production scoring"
    if gap == "route_repair_before_sample_import":
        return "repair clean route for existing positives; do not import duplicate labels"
    if gap == "source_backed_family_undercoverage":
        return "extract exact Rhea/PubMed substrate-product pairs, map structures, de-dup, then run generator gate"
    if gap == "gold_test_visibility_gap":
        return "promote already curated gold/silver positives into family-balanced evaluation where leakage-safe"
    if gap == "generator_recall_or_rule_pool_gap":
        return "trace why reactions are present but not generated; add source-traceable rule overlays only"
    if gap == "literature_context_needs_exact_pair_extraction":
        return "screen papers for explicit reactions; treat associations/pathway context as non-training evidence"
    return "monitor; keep as generic background unless production use case requires it"


def build_positive_queue(coverage: pd.DataFrame, rhea: pd.DataFrame, pubmed: pd.DataFrame) -> pd.DataFrame:
    rows = []
    high = set(coverage[coverage["priority_level"].isin(["P0", "P1"])]["family"])
    for rec in rhea.itertuples(index=False):
        if rec.priority_family not in high:
            continue
        if rec.candidate_class not in {"approved_exact_or_specific_reaction_lead", "generic_rule_template_context"}:
            continue
        rows.append(
            {
                "round": 17,
                "source_type": "Rhea",
                "priority_family": rec.priority_family,
                "source_id": rec.rhea_id,
                "candidate_text": rec.equation,
                "candidate_status": rec.candidate_class,
                "training_allowed_round17": False,
                "required_next_verification": "Map ChEBI/PubChem structures; de-dup local pairs; check clean generator route; only then decide import or rule repair.",
                "raw_file": rec.raw_file,
            }
        )
    for rec in pubmed.itertuples(index=False):
        if rec.priority_family not in high:
            continue
        if rec.candidate_class not in {"enzyme_or_mechanism_lead", "literature_lead_for_manual_exact_pair_review"}:
            continue
        rows.append(
            {
                "round": 17,
                "source_type": "PubMed",
                "priority_family": rec.priority_family,
                "source_id": f"PMID:{rec.pmid}",
                "candidate_text": rec.title,
                "candidate_status": rec.candidate_class,
                "training_allowed_round17": False,
                "required_next_verification": "Read full paper/abstract for explicit single-step substrate-product conversion, microbe/strain, assay condition, and structure mapping.",
                "raw_file": rec.raw_file,
            }
        )
    return pd.DataFrame(rows)


def build_negative_logic() -> pd.DataFrame:
    rows = []
    for fam in FAMILY_ORDER:
        if fam in {"other_or_unclassified", "polyphenol_generic", "carbohydrate_generic", "fatty_lipid_generic", "protein_peptide_amino_acid_generic"}:
            continue
        rows.extend(
            [
                {
                    "round": 17,
                    "priority_family": fam,
                    "negative_type": "hard_decoy_rule_generated_nontruth",
                    "safe_for_training_scope": "ranking_negative_only",
                    "source_requirements": "Generated by the same source-backed module/family rule set and absent from local positives after leakage check.",
                    "forbidden_shortcut": "Do not call a generator miss a biological negative.",
                    "example_use": "Within-substrate candidate ranking where all candidates are chemically plausible outputs.",
                },
                {
                    "round": 17,
                    "priority_family": fam,
                    "negative_type": "explicit_assay_no_conversion",
                    "safe_for_training_scope": "biological_negative_after_manual_review",
                    "source_requirements": "Paper/database must explicitly report substrate, organism/strain/community, condition/time, analytical method, and no detected product.",
                    "forbidden_shortcut": "Do not infer no-conversion from absence in Rhea/PubMed or absence from a metabolic model.",
                    "example_use": "Separate calibration/evaluation set, not mixed silently with hard decoys.",
                },
            ]
        )
    return pd.DataFrame(rows)


def build_github_contract() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "round": 17,
                "source": "github_web_verified",
                "repo_or_dataset": "rxn4chemistry/biocatalysis-model",
                "url": "https://github.com/rxn4chemistry/biocatalysis-model",
                "coverage_contract": "ECREACT uses reaction SMILES with EC/source fields, aggregated from Rhea, BRENDA, PathBank, and MetaNetX.",
                "foodgut_gap": "FoodGut needs explicit EC/source provenance per candidate and product-disjoint/family-balanced evaluation.",
            },
            {
                "round": 17,
                "source": "github_web_verified",
                "repo_or_dataset": "hesther/enzymemap",
                "url": "https://github.com/hesther/enzymemap",
                "coverage_contract": "EnzymeMap emphasizes atom mapping, correction, validation, and processed reaction tables.",
                "foodgut_gap": "FoodGut should not import pathway text as labels; it needs structure validation and reaction correction before training.",
            },
            {
                "round": 17,
                "source": "github_web_verified",
                "repo_or_dataset": "jotech/gapseq",
                "url": "https://github.com/jotech/gapseq",
                "coverage_contract": "gapseq uses curated bacterial metabolic reaction databases, sequence evidence, gap-filling trace flags, and phenotype benchmarks.",
                "foodgut_gap": "FoodGut needs trace flags explaining every reaction/rule admission and phenotype/literature validation beyond ranker metrics.",
            },
        ]
    )


def build_decisions(coverage: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in coverage.itertuples(index=False):
        if r.priority_level == "P0":
            decision = "block_production_until_family_gate_resolved"
            allowed = False
        elif r.priority_level == "P1":
            decision = "curate_next_before_retraining"
            allowed = False
        else:
            decision = "monitor_or_generic_background"
            allowed = False
        rows.append(
            {
                "round": 17,
                "priority_family": r.family,
                "priority_level": r.priority_level,
                "dominant_gap_type": r.dominant_gap_type,
                "training_allowed_round17": allowed,
                "decision": decision,
                "recommended_next_action": r.next_action,
            }
        )
    return pd.DataFrame(rows)


def write_doc(coverage: pd.DataFrame, rhea: pd.DataFrame, pubmed: pd.DataFrame, positive: pd.DataFrame) -> None:
    metrics = pd.read_csv(METRICS)
    npts = sorted(metrics["n_pts"].unique().tolist()) if "n_pts" in metrics.columns else []
    top = coverage.sort_values(["priority_level", "priority_score"], ascending=[True, False]).head(8)
    rhea_counts = rhea.groupby("candidate_class").size().to_dict()
    pubmed_counts = pubmed.groupby("candidate_class").size().to_dict()

    doc = f"""# Round17 Reaction Family Coverage Audit

## Working conclusion

The production blocker is not just total sample count. It is a combined family-coverage and route-recall problem:

1. Some source-backed positives already exist locally but are not visible to final clean candidates.
2. Several gut-microbiome reaction families have external evidence but weak final clean/gold representation.
3. Current metrics are still too small for production interpretation (`clean2_metrics.csv` n_pts={npts}).

## Highest-priority families

```text
{top[['family', 'priority_level', 'priority_score', 'dominant_gap_type', 'clean_pairs', 'clean_gold_pairs', 'next_action']].to_string(index=False)}
```

## External retrieval summary

Rhea Round17:

```text
{json.dumps(rhea_counts, ensure_ascii=False, indent=2)}
```

PubMed Round17:

```text
{json.dumps(pubmed_counts, ensure_ascii=False, indent=2)}
```

Rhea gives exact or generic reaction support for glucuronide/sulfate deconjugation, tryptophan/indole transformations, and choline/carnitine to TMA. PubMed gives many leads for tryptophan/indole and TMA biology, but most are association/pathway contexts and must not become labels without exact reaction extraction.

## Positive-label gate

Round17 keeps all new external leads as `training_allowed_round17=False`. A lead can move forward only after:

1. exact substrate and product are extracted,
2. ChEBI/PubChem/InChIKey mapping passes,
3. the pair is de-duplicated against the raw pool and `reactions.parquet`,
4. the clean generator can route the candidate or a source-traceable rule overlay is justified,
5. the label is assigned to a family-balanced split without product leakage.

## Negative-label gate

No-hit is not a biological negative. Use rule-generated non-truth products as ranking decoys only. Use assay negatives only when a paper or database explicitly reports no conversion under defined conditions.

## Written artifacts

- `data/curation/round17_reaction_family_coverage_matrix.csv`
- `data/curation/round17_rhea_query_triage.csv`
- `data/curation/round17_pubmed_query_triage.csv`
- `data/curation/round17_positive_expansion_queue.csv`
- `data/curation/round17_negative_sampling_logic.csv`
- `data/curation/round17_github_external_coverage_contract.csv`
- `data/curation/round17_training_decisions.csv`
"""
    DOCS.mkdir(parents=True, exist_ok=True)
    OUT_DOC.write_text(doc, encoding="utf-8")


def main() -> None:
    CURATION.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    raw_map = build_raw_pair_map()
    rhea = parse_rhea17()
    pubmed = parse_pubmed17()
    coverage = build_coverage_matrix(raw_map, rhea, pubmed)
    positive = build_positive_queue(coverage, rhea, pubmed)
    negative = build_negative_logic()
    github = build_github_contract()
    decisions = build_decisions(coverage)

    coverage.to_csv(OUT_COVERAGE, index=False)
    rhea.to_csv(OUT_RHEA, index=False)
    pubmed.to_csv(OUT_PUBMED, index=False)
    positive.to_csv(OUT_POSITIVE, index=False)
    negative.to_csv(OUT_NEGATIVE, index=False)
    github.to_csv(OUT_GITHUB, index=False)
    decisions.to_csv(OUT_DECISIONS, index=False)
    write_doc(coverage, rhea, pubmed, positive)

    print(f"wrote {OUT_COVERAGE.relative_to(REPO)} rows={len(coverage)}")
    print(f"wrote {OUT_RHEA.relative_to(REPO)} rows={len(rhea)}")
    print(f"wrote {OUT_PUBMED.relative_to(REPO)} rows={len(pubmed)}")
    print(f"wrote {OUT_POSITIVE.relative_to(REPO)} rows={len(positive)}")
    print(f"wrote {OUT_NEGATIVE.relative_to(REPO)} rows={len(negative)}")
    print(f"wrote {OUT_DOC.relative_to(REPO)}")


if __name__ == "__main__":
    main()
