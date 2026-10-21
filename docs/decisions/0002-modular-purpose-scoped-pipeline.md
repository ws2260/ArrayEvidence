# ADR-0002: Modular, purpose-scoped ArrayEvidence pipeline

Date: 2026-10-21
Status: Accepted

## Context

Genotyping-array data have multiple legitimate research uses,
including single-cell demultiplexing, imputation, sc-eQTL
analysis, and phenotype-specific genetic research.

These applications do not necessarily require biomedical
variant annotation or clinical interpretation.

Furthermore, research consent, ethics approval, and data-use
agreements may restrict permitted genetic analyses.

A mandatory end-to-end clinical interpretation pipeline would
therefore be scientifically unnecessary for many applications
and potentially incompatible with approved research purposes.

## Decision

ArrayEvidence will use independently executable modules.

The technical evidence and genotype reconciliation module
will be a standalone product.

Genotype utility/export, biomedical annotation, and
biomedical interpretation will be separate optional modules.

Analysis scope and output permissions will be explicit.

Module interfaces will preserve provenance, uncertainty,
and evidence separation.

The pipeline will not assume that possession of genetic
data implies authorization for unrestricted analysis.

## Consequences

Benefits:

- Independent utility of technical genotype reconciliation
- Support for demultiplexing, imputation, and sc-eQTL workflows
- Reduced unnecessary biomedical analysis
- Easier validation and testing
- Clearer publication and release boundaries
- Better support for purpose-limited research

Costs:

- Versioned interfaces and schemas are required
- Export policies must be maintained
- Scope enforcement needs dedicated design and testing
- Different modules may require independent validation

## Implementation implications

The existing reconciliation components remain the foundation
of Module 1.

Before integrating new functionality, define interfaces
that do not depend on biomedical annotation.

Clinical annotation must not become a dependency of
genotype reconciliation.

Purpose-scoped processing and restricted output generation
must be treated as architectural requirements rather than
optional reporting filters.

See docs/ARCHITECTURE.md.
