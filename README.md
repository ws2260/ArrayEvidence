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
- Structured technical evidence tables
- Bundle integrity verification

Module 1 produces technical evidence, not downstream genotype-format exports.

Module 1 does not require clinical or disease annotation.

### Module 2A — Genotype utility (under consideration)

Potential genotype-format export capabilities, including VCF/BCF, for:

- Single-cell donor demultiplexing
- Preparation of inputs for external genotype imputation
- sc-eQTL and genetic association analyses

These capabilities are not implemented or committed to a release.

### Module 2B — Variant annotation (under consideration)

Optional, authorization-dependent biomedical annotations.

### Module 3 — Biomedical interpretation (under consideration)

Optional phenotype-specific prioritization, evidence adjudication,
and research-oriented reporting.

## Purpose-limited analysis

Not every study permits broad genetic interpretation.

ArrayEvidence is designed to support technical-only workflows
and explicitly scoped biomedical analyses.

Selecting a module does not itself establish ethics, consent,
or data-use authorization.

## Module 1 v0.1.0 — Implemented functionality

The released Module 1 v0.1.0 provides:

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

## Recommended upstream genotype calling

For original Illumina IDAT data, see the [Local DRAGEN Array IDAT-to-VCF vignette](docs/vignettes/dragen_array_from_idat.md).

DRAGEN Array is separate software; ArrayEvidence starts downstream of genotype calling.

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

## Optional Nextflow workflow

`workflows/from_vcf.nf` runs BUILD, VERIFY, QC, and atomic staging from an existing assay manifest, single-sample VCF, and indexed FASTA.

Required parameters: `--manifest`, `--vcf`, `--reference`, `--arrayevidence_bin`, and `--staging_dir`. Optional: `--sample_id` (default `PSEUDONYM_001`).

Nextflow and Java are separate dependencies. Python 3 and the ArrayEvidence CLI must also be available to the workflow. Use `--python_bin /absolute/path/to/python3` to select the interpreter for staging reservation and finalization; the default is `python3`.

The staging parent must exist; a fresh staging destination must not. Successful runs publish `validated/bundle` and `validated/qc.json` without overwriting an existing destination.

Resume (`-resume SESSION_ID`) is permitted only for the same session while `validated` is absent. Completed sessions require a new staging destination.

Atomic publication requires Linux `renameat2` and a compatible filesystem. It does not guarantee power-loss durability or protection against malicious concurrent filesystem modification.

Run synthetic Nextflow tests with `ARRAYEVIDENCE_RUN_NEXTFLOW=1 python -m pytest -q`. The test runner currently expects a conda environment named `arrayevidence-nextflow` containing Nextflow.

## Interpretation limitations

- `VCF_ABSENT` does not establish homozygous reference.
- `NOT_EVALUATED` does not imply a negative genotype result.
- Allele `MATCH` is not genotype concordance.
- Ordinary biallelic SNV genotype projection is not evaluated.
- Supported multiallelic SNV and validated indel scenarios have
  more specific genotype-projection behavior.
- Technical evidence does not establish clinical validity.

## Technical validation

Module 1 currently passes 402 automated tests (including opt-in Nextflow tests).

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
- [Validation status](docs/VALIDATION_STATUS.md)
- [Qdenga independent application](docs/qdenga_independent_application.md)
- [Data policy](docs/DATA_POLICY.md)
- [ADR-0001: Indel candidate discovery](docs/decisions/0001-indel-candidate-discovery.md)
- [ADR-0002: Modular purpose-scoped pipeline](docs/decisions/0002-modular-purpose-scoped-pipeline.md)

## Project status

Research software under active development.

Technical reconciliation is the first independently releasable
component. Biomedical annotation and interpretation are
optional future extensions.
