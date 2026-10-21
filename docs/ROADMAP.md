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

Current status: Indel candidate discovery checkpoint reached.

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
- Export policies
- Reconciled VCF/BCF
- QC and provenance reports

## Milestone 4: Downstream utility

- Demultiplexing-ready export
- Imputation-ready export
- sc-eQTL/association-ready export
- Application-specific benchmarking

## Milestone 5: Independent validation

- Freeze technical methods and policies
- Evaluate held-out authorized GSA-WGS data
- Quantify concordance, discordance, coverage,
  missingness, and limitations

## Milestone 6: Optional biomedical extensions

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
