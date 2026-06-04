from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CUR = ROOT / "data" / "curation"
RUNTIME = ROOT / "runtime" / "external_round24"
ROUND = "round24"


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def write_csv(df: pd.DataFrame, rel: str) -> Path:
    path = ROOT / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def build_external_rule_level_coverage() -> pd.DataFrame:
    p0 = read_csv(CUR / "round23_p0_source_verification.csv")
    ecreact = read_csv(CUR / "round23_ecreact_p0_exact_pair_screen.csv")
    enzymemap = read_csv(RUNTIME / "enzymemap_p0_exact_pair_screen_round24.csv")
    dryrun = read_csv(RUNTIME / "enzymemap_rule_dryrun_round24.csv")

    cols = [
        "candidate_id",
        "ecreact_forward_hits",
        "ecreact_reverse_hits",
    ]
    p0 = p0.merge(ecreact[cols], on="candidate_id", how="left")
    if not enzymemap.empty:
        keep = [
            "candidate_id",
            "enzymemap_forward_hits",
            "enzymemap_reverse_hits",
            "example_forward_rule_id",
            "example_reverse_rule_id",
            "example_forward_ec",
            "example_reverse_ec",
            "example_forward_source",
            "example_reverse_source",
            "example_forward_quality",
            "example_reverse_quality",
        ]
        p0 = p0.merge(enzymemap[keep], on="candidate_id", how="left")
    else:
        p0["enzymemap_forward_hits"] = "not_screened_runtime_missing"
        p0["enzymemap_reverse_hits"] = "not_screened_runtime_missing"

    dry_summary = []
    if not dryrun.empty:
        for cid, g in dryrun.groupby("candidate_id"):
            dry_summary.append(
                {
                    "candidate_id": cid,
                    "enzymemap_rule_dryrun_checks": len(g),
                    "enzymemap_dryrun_target_generated_count": int(g["target_generated"].sum()),
                    "enzymemap_dryrun_rule_ids_that_generate_target": ";".join(
                        sorted(
                            set(
                                g.loc[g["target_generated"], "enzymemap_rule_id"]
                                .dropna()
                                .astype(str)
                                .tolist()
                            )
                        )
                    ),
                    "enzymemap_dryrun_ecs_that_generate_target": ";".join(
                        sorted(
                            set(
                                g.loc[g["target_generated"], "enzymemap_ec"]
                                .dropna()
                                .astype(str)
                                .tolist()
                            )
                        )
                    ),
                }
            )
    p0 = p0.merge(pd.DataFrame(dry_summary), on="candidate_id", how="left")
    for col in [
        "enzymemap_rule_dryrun_checks",
        "enzymemap_dryrun_target_generated_count",
    ]:
        p0[col] = p0[col].fillna(0).astype(int)

    rows = []
    for r in p0.to_dict("records"):
        cid = r["candidate_id"]
        has_rule_generating_target = r.get("enzymemap_dryrun_target_generated_count", 0) > 0
        if cid == "R21-P0-009" and has_rule_generating_target:
            rule_level_status = "external_rule_generates_target_but_blocked_by_biology_and_generality_gate"
            next_gate = (
                "recurate gut-specific PMID/organism evidence, quantify rule overgeneration on aromatic acids, "
                "then run isolated overlay dry-run with false-negative screen"
            )
        elif cid == "R21-P0-012":
            rule_level_status = "enzymemap_reverse_reaction_found_but_wrong_biological_context_and_no_target_generation"
            next_gate = "recurate gut bacterial 21-dehydroxylase evidence; do not use Bos taurus hydroxylase rule"
        elif "urolithin" in str(r.get("reaction_family", "")):
            rule_level_status = "literature_supported_but_no_rhea_ecreact_enzymemap_rule_level_coverage"
            next_gate = "extract exact urolithin substrate/product table from Enterocloster papers and atom-map before rule derivation"
        elif cid == "R21-P0-002":
            rule_level_status = "exact_rhea_reaction_but_no_ecreact_enzymemap_rule_level_coverage"
            next_gate = "query/download RetroRules or derive mapped template from RHEA:61520 only after atom mapping"
        elif cid == "R21-P0-013":
            rule_level_status = "pubmed_supported_but_no_exact_enzymemap_ecreact_pair"
            next_gate = "verify exact hydrocaffeic acid dehydroxylation route and atom-mapped rule source"
        else:
            rule_level_status = "not_rule_level_supported_in_round24"
            next_gate = "targeted source query required before training"

        rows.append(
            {
                "round": ROUND,
                "candidate_id": cid,
                "module": r.get("module", ""),
                "reaction_family": r.get("reaction_family", ""),
                "substrate_name": r.get("substrate_name", ""),
                "product_name": r.get("product_name", ""),
                "round22_target_generated": r.get("round22_target_generated", ""),
                "rhea_ids_round23": r.get("rhea_ids_round23", ""),
                "ecreact_forward_hits_round23": r.get("ecreact_forward_hits", ""),
                "ecreact_reverse_hits_round23": r.get("ecreact_reverse_hits", ""),
                "enzymemap_forward_hits_round24": r.get("enzymemap_forward_hits", ""),
                "enzymemap_reverse_hits_round24": r.get("enzymemap_reverse_hits", ""),
                "enzymemap_example_forward_rule_id": r.get("example_forward_rule_id", ""),
                "enzymemap_example_reverse_rule_id": r.get("example_reverse_rule_id", ""),
                "enzymemap_example_forward_ec": r.get("example_forward_ec", ""),
                "enzymemap_example_reverse_ec": r.get("example_reverse_ec", ""),
                "enzymemap_rule_dryrun_checks": r.get("enzymemap_rule_dryrun_checks", 0),
                "enzymemap_dryrun_target_generated_count": r.get(
                    "enzymemap_dryrun_target_generated_count", 0
                ),
                "enzymemap_dryrun_rule_ids_that_generate_target": r.get(
                    "enzymemap_dryrun_rule_ids_that_generate_target", ""
                ),
                "enzymemap_dryrun_ecs_that_generate_target": r.get(
                    "enzymemap_dryrun_ecs_that_generate_target", ""
                ),
                "rule_level_status_round24": rule_level_status,
                "training_allowed_round24": False,
                "rule_promotion_allowed_round24": False,
                "next_gate_round24": next_gate,
            }
        )
    return pd.DataFrame(rows)


def build_rule_promotion_candidates(coverage: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in coverage.to_dict("records"):
        if r["candidate_id"] == "R21-P0-009":
            rows.append(
                {
                    "round": ROUND,
                    "candidate_id": r["candidate_id"],
                    "candidate_rule_source": "EnzymeMap",
                    "candidate_rule_ids": r["enzymemap_dryrun_rule_ids_that_generate_target"],
                    "candidate_ecs": r["enzymemap_dryrun_ecs_that_generate_target"],
                    "target_generated_in_current_run_reactants": bool(
                        r["enzymemap_dryrun_target_generated_count"]
                    ),
                    "why_not_promoted": (
                        "Rule is chemically broad decarboxylation; prior gut-specific PMID mismatch remains; "
                        "needs overgeneration and false-negative screen before any overlay."
                    ),
                    "minimum_promotion_gate": (
                        "Exact Rhea/EnzymeMap pair plus gut organism/assay/source PMID re-curation, "
                        "then dry-run against aromatic acid challenge set with decoy audit."
                    ),
                    "training_allowed_round24": False,
                    "rule_promotion_allowed_round24": False,
                }
            )
        elif r["candidate_id"] == "R21-P0-012":
            rows.append(
                {
                    "round": ROUND,
                    "candidate_id": r["candidate_id"],
                    "candidate_rule_source": "EnzymeMap",
                    "candidate_rule_ids": r["enzymemap_example_reverse_rule_id"],
                    "candidate_ecs": r["enzymemap_example_reverse_ec"],
                    "target_generated_in_current_run_reactants": False,
                    "why_not_promoted": (
                        "The EnzymeMap hit is a Bos taurus hydroxylase-like reverse context and did not generate "
                        "the target from corticosterone in current run_reactants."
                    ),
                    "minimum_promotion_gate": (
                        "Find gut bacterial 21-dehydroxylase source with exact substrate/product and mapped rule; "
                        "do not reuse mammalian hydroxylase as a gut dehydroxylation rule."
                    ),
                    "training_allowed_round24": False,
                    "rule_promotion_allowed_round24": False,
                }
            )
    return pd.DataFrame(rows)


def build_keyword_screen_summary() -> pd.DataFrame:
    kw = read_csv(RUNTIME / "enzymemap_keyword_screen_round24.csv")
    if kw.empty:
        return pd.DataFrame(
            [
                {
                    "round": ROUND,
                    "source": "EnzymeMap",
                    "term": "all",
                    "hit_count": "not_screened_runtime_missing",
                    "interpretation": "Runtime keyword screen was missing.",
                }
            ]
        )
    kw.insert(0, "round", ROUND)
    interpretations = []
    for row in kw.to_dict("records"):
        term = row["term"]
        count = int(row["hit_count"])
        if term == "urolithin" and count == 0:
            interp = "No EnzymeMap keyword coverage for urolithin; literature extraction is still required."
        elif term in {"dopamine", "tyramine"} and count > 0:
            interp = "Related dopamine/tyramine hydroxylation chemistry exists, but exact dopamine -> m-tyramine was not covered."
        elif term == "protocatechu" and count > 0:
            interp = "Protocatechuate chemistry is well represented; exact decarboxylation requires gut-context re-curation."
        elif term == "dehydroxyl" and count == 0:
            interp = "The word dehydroxylation is not represented in EnzymeMap text fields; exact structure/rule search is needed."
        else:
            interp = "Keyword evidence only; not a training label."
        interpretations.append(interp)
    kw["interpretation"] = interpretations
    return kw


def build_family_decision_matrix(coverage: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for family, g in coverage.groupby("reaction_family"):
        exact_rule_count = int((g["enzymemap_dryrun_target_generated_count"] > 0).sum())
        rhea_exact_count = int(g["rhea_ids_round23"].fillna("").astype(str).str.contains("RHEA:").sum())
        ecreact_hits = int(
            pd.to_numeric(g["ecreact_forward_hits_round23"], errors="coerce").fillna(0).sum()
            + pd.to_numeric(g["ecreact_reverse_hits_round23"], errors="coerce").fillna(0).sum()
        )
        enzymemap_hits = int(
            pd.to_numeric(g["enzymemap_forward_hits_round24"], errors="coerce").fillna(0).sum()
            + pd.to_numeric(g["enzymemap_reverse_hits_round24"], errors="coerce").fillna(0).sum()
        )
        if family == "decarboxylation":
            decision = "source_backed_rule_candidate_exists_but_not_promoted"
        elif family == "polyphenol_urolithin":
            decision = "literature_supported_rule_source_missing"
        elif family == "catechol_dehydroxylation":
            decision = "exact_rhea_reaction_rule_source_missing"
        elif family == "functional_group_removal":
            decision = "mixed_evidence_mostly_rule_source_missing"
        else:
            decision = "requires_targeted_rescreen"
        rows.append(
            {
                "round": ROUND,
                "reaction_family": family,
                "p0_candidate_count": len(g),
                "rhea_exact_candidate_count": rhea_exact_count,
                "ecreact_pair_hit_total": ecreact_hits,
                "enzymemap_pair_hit_total": enzymemap_hits,
                "enzymemap_rule_dryrun_target_generating_candidate_count": exact_rule_count,
                "round24_decision": decision,
                "training_allowed_round24": False,
                "next_best_action": (
                    "Promote nothing automatically; use the most specific source-backed rule candidate "
                    "only after family-level overgeneration, false-negative, and biological-context review."
                ),
            }
        )
    return pd.DataFrame(rows)


def build_negative_mining_update() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "round": ROUND,
                "negative_source": "EnzymeMap/ECReact exact-pair absence",
                "allowed_use": "unknown_only",
                "blocked_use": "hard_negative",
                "reason": "Absence from a public reaction dataset can reflect curation scope, not biological impossibility.",
                "next_step": "Use absence only as a false-negative screen flag for generated decoys.",
            },
            {
                "round": ROUND,
                "negative_source": "EnzymeMap reverse or non-gut organism hit",
                "allowed_use": "direction_or_context_warning",
                "blocked_use": "gut_microbe_negative_or_positive",
                "reason": "A mammalian or reverse-direction enzyme hit does not establish gut microbial production or non-production.",
                "next_step": "Record organism, EC, direction, and source before label assignment.",
            },
            {
                "round": ROUND,
                "negative_source": "Rule-generated aromatic acid decoys from EnzymeMap rule 92",
                "allowed_use": "ranking_decoy_after_false_negative_screen",
                "blocked_use": "biological_false_reaction",
                "reason": "Rule 92 is broad decarboxylation and may generate plausible but uncurated metabolites.",
                "next_step": "Screen each generated decoy against Rhea, EnzymeMap, ECReact, current positives, and PubMed before y=0 hardening.",
            },
            {
                "round": ROUND,
                "negative_source": "Assay non-detection in literature",
                "allowed_use": "condition_specific_negative",
                "blocked_use": "universal_negative",
                "reason": "Non-detection depends on strain, medium, assay, time, and detection limit.",
                "next_step": "Mine papers for tested substrate-product pairs with explicit not-detected result fields.",
            },
        ]
    )


def build_query_manifest() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "round": ROUND,
                "source": "EnzymeMap GitHub raw",
                "query_or_file": "data/processed_reactions.csv.gz",
                "url": "https://raw.githubusercontent.com/hesther/enzymemap/main/data/processed_reactions.csv.gz",
                "local_runtime_path": "runtime/external_round24/enzymemap_processed_reactions.csv.gz",
                "tracked_summary_path": "data/curation/round24_enzymemap_p0_exact_pair_screen.csv",
                "status": "downloaded_and_screened",
                "result_summary": "349,458 rows; exact pair coverage found for R21-P0-009 and reverse/context hit for R21-P0-012.",
            },
            {
                "round": ROUND,
                "source": "RetroRules API",
                "query_or_file": "/api/templates?ec=1.17.99&radius=4 and /api/templates?ec=4.1.1&radius=4",
                "url": "https://retrorules.org/docs",
                "local_runtime_path": "",
                "tracked_summary_path": "",
                "status": "timed_out_round24",
                "result_summary": "Docs support SMARTS/template/radius/AAM workflow, but live API requests timed out in this run.",
            },
            {
                "round": ROUND,
                "source": "ECReact Zenodo cache from Round23",
                "query_or_file": "ecreact-nofilter-1.0.csv.gz exact P0 pair screen",
                "url": "https://zenodo.org/records/8318231",
                "local_runtime_path": "runtime/external_round23/ecreact-nofilter-1.0.csv.gz",
                "tracked_summary_path": "data/curation/round23_ecreact_p0_exact_pair_screen.csv",
                "status": "screened_round23_reused",
                "result_summary": "Exact ECReact pair coverage only for R21-P0-009.",
            },
        ]
    )


def build_review_md() -> str:
    return "\n".join(
        [
            "# Round24 Rule-level Coverage Review",
            "",
            "## Working conclusion",
            "",
            "Round24 moved from reaction-level evidence to rule-level evidence. The main production blocker still sits before ranking: most P0 families do not have a source-backed rule/template that current `run_reactants` can use.",
            "",
            "## Rule-level findings",
            "",
            "- EnzymeMap `processed_reactions.csv.gz` was downloaded into runtime and screened at exact substrate/product level.",
            "- `protocatechuic acid -> catechol` has EnzymeMap exact pair coverage and rule 92 generates the target in current `run_reactants`.",
            "- Rule 92 is broad aromatic acid decarboxylation and is not promoted because the gut-specific evidence and overgeneration risk are unresolved.",
            "- Urolithin P0 reactions have no Rhea/ECReact/EnzymeMap exact rule-level coverage in this run despite strong recent PubMed support.",
            "- `dopamine -> m-tyramine` has exact Rhea support but no EnzymeMap/ECReact exact pair or rule-level support in this run.",
            "- `corticosterone -> 11beta-hydroxyprogesterone` has an EnzymeMap reverse/context hit from Bos taurus hydroxylase chemistry, but it did not generate the target and is not microbiome-appropriate.",
            "",
            "## Training gate",
            "",
            "No Round24 sample or rule is admitted to training. The one target-generating rule candidate must first pass biological re-curation, aromatic acid overgeneration review, and false-negative screening.",
            "",
            "## Generated artifacts",
            "",
            "* `data/curation/round24_external_rule_level_coverage.csv`",
            "* `data/curation/round24_rule_promotion_candidates.csv`",
            "* `data/curation/round24_enzymemap_p0_exact_pair_screen.csv`",
            "* `data/curation/round24_enzymemap_rule_dryrun.csv`",
            "* `data/curation/round24_enzymemap_keyword_screen.csv`",
            "* `data/curation/round24_family_decision_matrix.csv`",
            "* `data/curation/round24_negative_mining_update.csv`",
            "* `data/curation/round24_external_query_manifest.csv`",
            "",
        ]
    )


def main() -> None:
    paths = {}
    coverage = build_external_rule_level_coverage()
    paths["coverage"] = write_csv(coverage, "data/curation/round24_external_rule_level_coverage.csv")
    paths["promote"] = write_csv(
        build_rule_promotion_candidates(coverage),
        "data/curation/round24_rule_promotion_candidates.csv",
    )
    paths["keyword"] = write_csv(
        build_keyword_screen_summary(), "data/curation/round24_enzymemap_keyword_screen.csv"
    )
    exact = read_csv(RUNTIME / "enzymemap_p0_exact_pair_screen_round24.csv")
    paths["exact"] = write_csv(exact, "data/curation/round24_enzymemap_p0_exact_pair_screen.csv")
    dry = read_csv(RUNTIME / "enzymemap_rule_dryrun_round24.csv")
    paths["dry"] = write_csv(dry, "data/curation/round24_enzymemap_rule_dryrun.csv")
    paths["family"] = write_csv(
        build_family_decision_matrix(coverage), "data/curation/round24_family_decision_matrix.csv"
    )
    paths["negatives"] = write_csv(
        build_negative_mining_update(), "data/curation/round24_negative_mining_update.csv"
    )
    paths["queries"] = write_csv(
        build_query_manifest(), "data/curation/round24_external_query_manifest.csv"
    )
    review = ROOT / "docs/reviews/round24_rule_level_coverage.md"
    review.write_text(build_review_md(), encoding="utf-8")
    paths["review"] = review

    for name, path in paths.items():
        print(f"{name}: {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
