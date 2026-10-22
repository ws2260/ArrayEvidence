# ArrayEvidence

**ArrayEvidence** is a deterministic, privacy-first, assay-aware
pipeline for evaluating genetic evidence from genotyping arrays.

The initial focus is Illumina Global Screening Array (GSA)
data and DRAGEN Array VCF outputs.

## Goals

ArrayEvidence aims to:

1. Reconcile array assay definitions with genotype calls.
2. Validate allele identity against a reference genome.
3. Preserve genotype uncertainty, assayability, and provenance.
4. Produce technically interpretable genotype datasets.
5. Support downstream genetic analyses through independent modules.
6. Enable optional, purpose-scoped biomedical interpretation.

ArrayEvidence starts downstream of genotype calling. It does not
replace DRAGEN Array or other primary genotype callers.

## Modular architecture

### Module 1 — Technical genotype evidence

The standalone core provides:

- Manifest-to-VCF reconciliation
- Reference and strand validation
- SNV and indel allele matching
- Multiallelic genotype projection
- Technical QC and provenance
- Policy-filtered genotype exports

Planned outputs include evidence tables, VCF/BCF, and QC reports.

Module 1 does not require clinical or disease annotation.

### Module 2A — Genotype utility

Optional application-specific exports for:

- Single-cell donor demultiplexing
- Genotype imputation
- sc-eQTL and genetic association analyses

### Module 2B — Variant annotation

Optional, authorization-dependent biomedical annotations.

### Module 3 — Biomedical interpretation

Optional phenotype-specific prioritization, evidence adjudication,
and research-oriented reporting.

## Purpose-limited analysis

Not every study permits broad genetic interpretation.

ArrayEvidence is designed to support technical-only workflows
and explicitly scoped biomedical analyses.

Selecting a module does not itself establish ethics, consent,
or data-use authorization.

## Module 1 v0.1.0 — Implemented functionality

The standalone Module 1 release candidate provides:

- Assay manifest and single-sample VCF parsing
- Reference-aware SNV and indel allele reconciliation
- Supported multiallelic SNV and validated indel genotype projection
- Deterministic, immutable evidence-bundle export
- Bundle integrity verification
- Read-only aggregate technical QC

The evidence bundle contains:

- `assay_evidence.tsv`
- `allele_candidates.tsv`
- `bundle_manifest.json`

Module 1 starts downstream of genotype calling. It does not call
genotypes from IDAT files.

General-purpose genotype VCF/BCF export, imputation, sc-eQTL analysis,
biomedical annotation, and clinical interpretation are not implemented
in this release.

## Installation

Requires Python 3.10 or newer.

From a clean clone, create a virtual environment and install:

    python -m venv .venv
    source .venv/bin/activate
    python -m pip install -e .

For synthetic regression tests:

    python -m pip install pytest
    python -m pytest -q

The package requires `pysam>=0.23,<0.24`.

## Command-line interface

Build a new immutable evidence bundle:

    arrayevidence build --manifest MANIFEST.csv --vcf SAMPLE.vcf --reference GRCh38.fa --sample-id PSEUDONYM --output NEW_BUNDLE

The reference FASTA must have an existing compatible index.
The output directory must not already exist.

Verify an existing bundle:

    arrayevidence verify NEW_BUNDLE

Produce deterministic aggregate QC JSON:

    arrayevidence qc NEW_BUNDLE

The QC command reports aggregate evidence categories, not
individual-level genotypes.

## Interpretation limitations

- `VCF_ABSENT` does not establish homozygous reference.
- `NOT_EVALUATED` does not imply a negative genotype result.
- Allele `MATCH` is not genotype concordance.
- Ordinary biallelic SNV genotype projection is not evaluated.
- Supported multiallelic SNV and validated indel scenarios have
  more specific genotype-projection behavior.
- Technical evidence does not establish clinical validity.

## Technical validation

Module 1 currently passes 373 synthetic regression tests.

Technical integration validation covered 48 independently generated
and verified AIDA evidence bundles:

- 31,393,296 assay evidence rows
- 31,013,136 allele candidate rows
- 98.428% allele `MATCH` evidence
- 97.844% `VCF_CALLED` coordinate evidence

These are technical evidence distributions, not estimates of
genotype accuracy or concordance against an independent truth set.

See [AIDA technical validation](docs/aida_technical_validation.md).

Qdenga vaccination-cohort data were not used to develop Module 1
and are reserved for a subsequent independent application.

## Development and validation

- Synthetic fixtures: software correctness and regression tests
- AIDA: authorized development and technical characterization
- Held-out paired GSA-WGS data: planned independent validation
  after method freeze

Participant-derived genetic data must not be committed to
the publication source repository.

## Documentation

- [Architecture](docs/ARCHITECTURE.md)
- [Technical evidence contract](docs/EVIDENCE_CONTRACT.md)
- [Development roadmap](docs/ROADMAP.md)
- [Data policy](docs/DATA_POLICY.md)
- [ADR-0001: Indel candidate discovery](docs/decisions/0001-indel-candidate-discovery.md)
- [ADR-0002: Modular purpose-scoped pipeline](docs/decisions/0002-modular-purpose-scoped-pipeline.md)

## Project status

Research software under active development.

Technical reconciliation is the first independently releasable
component. Biomedical annotation and interpretation are
optional future extensions.
