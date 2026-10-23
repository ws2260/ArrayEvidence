# ArrayEvidence Development Roadmap

Status: Planning document
Updated: 2026-10-21

## Milestone 1: Technical reconciliation foundation

- Manifest and VCF parsing
- Reference and strand validation
- SNV and multiallelic reconciliation
- Indel event reconstruction and normalization
- Normalized allele candidate discovery
- Synthetic tests and AIDA development audits

Historical development milestone; Module 1 v0.1.0 has subsequently been released.

## Milestone 2: Allele-aware reconciliation integration

- Explicit VCF record and ALT allele identity
- Allele-level candidate indexing
- Genotype projection with provenance
- Ambiguous and missing-evidence states
- Integration with reconcile()
- End-to-end synthetic regression tests

Development approach: Red-green-refactor TDD.

## Milestone 3: Standalone Module 1 outputs

- Versioned evidence schema
- Technical QC
- Immutable evidence bundles
- Bundle integrity verification
- Technical QC and provenance reporting

## Milestone 4: Downstream utility (under consideration)

- General-purpose VCF/BCF export under explicit policies
- Demultiplexing-ready export
- Imputation-ready export
- sc-eQTL/association-ready export
- Application-specific benchmarking

## Milestone 5: Independent validation

- Freeze technical methods and policies
- Evaluate held-out authorized GSA-WGS data
- Quantify concordance, discordance, coverage,
  missingness, and limitations

## Milestone 6: Optional biomedical extensions (under consideration)

- Purpose-scoped variant annotation
- Phenotype-specific prioritization
- Evidence adjudication
- Research-oriented reporting

These modules are not prerequisites for Module 1.

## Development rules

1. Define invariants before implementation.
2. Write synthetic tests before production code.
3. Use AIDA to evaluate assumptions, not optimize blindly.
4. Record architectural decisions.
5. Preserve research-data governance boundaries.
6. Keep modules independently usable.
7. Freeze methods before held-out validation.
