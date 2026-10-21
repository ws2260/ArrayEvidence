# ArrayEvidence Data Policy

Status: Active
Updated: 2026-10-21

## Core rule

The publication source repository must not contain
participant-derived genotype, sequencing, or raw array data.

This prohibition applies even when records are de-identified
or sample identifiers have been removed.

## Data that must remain outside the repository

Examples include:

- IDAT and GTC files
- Participant-derived VCF and BCF files
- PLINK genotype datasets
- Sequencing reads and participant-level sequencing outputs
- Participant-derived genotype or phenotype tables
- Real participant-derived test fixtures
- Local configuration containing sensitive dataset paths

Public reference resources and proprietary assay resources
must be handled according to their respective licenses
and distribution permissions.

## Synthetic testing

Automated tests must use deliberately synthetic data
generated independently of participant records.

Synthetic fixtures must not be constructed by copying,
subsetting, anonymizing, or modifying real participant data.

## Development datasets

AIDA may be used for authorized local development and
technical characterization.

Its participant-derived data must remain outside the
publication repository.

Development results are not independent validation.

## Held-out validation

Authorized paired GSA-WGS data are reserved for independent
evaluation after method and analysis-policy freeze.

Held-out data must not be used to tune matching rules
or optimize development thresholds.

## Purpose limitation

Possession of genetic data does not imply permission for
unrestricted genetic analysis.

The pipeline must support technical-only workflows and
appropriately restricted analysis scopes.

Scope configuration alone is not a substitute for
institutional authorization or access controls.

## Repository safeguards

Before committing:

1. Inspect all staged files.
2. Verify that no participant-derived data are included.
3. Verify that no credentials or sensitive local paths appear.
4. Review newly added fixtures for synthetic provenance.
5. Avoid blanket staging commands when file safety is uncertain.

Git ignore rules are supplementary safeguards, not a
replacement for reviewing repository contents.

## Publication separation

The publication-track implementation must remain independent
of legacy datasets or code whose reuse or publication
authorization has not been established.
