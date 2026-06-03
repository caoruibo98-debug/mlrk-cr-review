# Scientific Positioning

## Narrow Safe Position

The current model should be positioned as a food-polyphenol metabolite candidate generation and ranking system for internal research prioritization. Its strongest demonstrated use is glycoside-to-aglycone recovery and ranking for common dietary flavonoids and isoflavones.

It should not yet be positioned as a broad gut microbiome metabolism simulator, strain-aware predictor, or clinical/consumer health product.

## Comparison Lanes

| Lane | Representative systems | What they test | How this model compares today |
| --- | --- | --- | --- |
| Small-molecule biotransformation prediction | BioTransformer 3/4 | Product generation, metabolism modules, enzyme annotations | Comparable in spirit only; needs same-panel execution and InChIKey recall comparison. |
| Enzyme/reaction-based gut metabolite prediction | MicrobeRX, GutBug | Enzyme or reaction plausibility in gut microbes | Our evidence overlay is lighter and not genome-scale; compare product recall and EC/PubMed traceability separately. |
| Community metabolic potential | MIMOSA2, AGORA, AGREDA | Microbiome sample context, taxa/reactions, metabolite potential | Out of current scope; useful as a coverage and pathway-evidence reference, not as a rank@k competitor. |
| Food phenolic pathway reconstructions | AGREDA and enzyme-promiscuity extensions | Dietary phenolic pathway breadth | Strong target for gap analysis because our current panel is mostly glycoside hydrolysis. |

## Evidence Anchors

- BioTransformer 3.0 is described as a web server for small-molecule metabolism prediction with gut microbial and other modules; BioTransformer 4.0 is reported as a successor with seven modules.
- MicrobeRX is described as an enzymatic-reaction-based gut microbiome metabolite prediction tool using genome-scale metabolic model derived reaction rules and drug metabolic reactions.
- MIMOSA2 builds community metabolic potential from microbiome features and reaction databases, but it addresses paired microbiome-metabolome interpretation rather than single-substrate product ranking.
- AGREDA extends human gut microbiota diet metabolism reconstruction, including phenolic compounds, making it a useful pathway-coverage reference for food-compound metabolism.

## Current Biological Logic Check

The existing model logic is biologically plausible when:

1. the expected metabolite is reachable by encoded reaction rules or the conservative glycoside rescue,
2. ranking is interpreted only within the generated candidate set,
3. evidence fields are treated as explanatory support rather than training features,
4. output quality filters reject chemically implausible or unsupported products,
5. the claim remains internal research prioritization.

The logic becomes weak when:

1. the true pathway requires multi-step ring cleavage, dehydroxylation, reduction, demethylation, or cross-feeding,
2. the reaction depends on strain-level genes or community context,
3. a product is expected because of literature evidence but absent from the candidate set,
4. a user interprets normalized ranking scores as wet-lab probabilities,
5. the evaluation panel stays too small or chemically homogeneous.

## Production Direction

The next production-safe direction is not to broaden claims first. It is to broaden the evaluation panel and failure taxonomy first, then decide whether the model can support:

1. food-polyphenol glycoside triage,
2. broader phenolic ester/lignan/isoflavone triage,
3. multi-step microbiome transformation exploration,
4. strain- or gene-aware prediction.

Only the first category is currently supported with strong internal evidence.
