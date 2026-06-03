# Generation 15 Review

## Production question

Can repository layout be standardized without breaking legacy scientific scripts?

## Changes

1. Added `docs/REPOSITORY_GUIDE.md`.
2. Added `tools/repo_doctor.py` and `scripts/repo_doctor.py`.
3. Added `outputs/appraisal/repo_doctor.json`.
4. Added layout contract tests for:
   - repository doctor readiness,
   - canonical reviewer commands,
   - script/tool entrypoint pairing,
   - preserved legacy `modular/` and `src/` boundaries,
   - repository guide coverage.
5. Added a `repository_contract` section to `production_artifact_manifest.json`.
6. Integrated the repository doctor report into the production scorecard.

## Results

- Repository doctor: `ready`
- Repository doctor checks: `33`
- Repository doctor failures: `0`
- Contract tests: `PASSED 43`
- Production scorecard status: `internal_mvp_only`
- Core strict top-5 hit rate: `6 / 6`
- Challenge strict top-5 hit rate: `8 / 22`
- Readiness: `internal_mvp_only`

## Audit interpretation

This generation makes the repository easier to review and maintain. It does not move legacy scientific code and does not claim improved biological prediction.

The key production improvement is that reviewer-facing commands now have an explicit policy:

> `scripts/` contains the commands a reviewer runs; `tools/` contains implementation; legacy `modular/` and `src/` code stays in place unless migrated in a dedicated compatibility round.

This makes the codebase look less like a one-off generated bundle and more like a maintained project with stable entrypoints.

## Remaining gap

The layout is documented and checked, but the legacy modules are still large and domain-specific. Future migration should be gradual and compatibility-tested.

## Next round

Generation 16 should make the README reviewer-first and GitHub-ready: installation, quickstart, expected outputs, limitations, and reviewer commands should be visible without reading internal history.
