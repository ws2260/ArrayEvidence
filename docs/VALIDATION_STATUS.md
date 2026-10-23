# ArrayEvidence Validation Status

This document distinguishes software correctness, technical
integration, independent application, and genotype accuracy.

## Synthetic regression testing

Automated synthetic tests assess defined software behavior and
regression safety. Passing tests does not establish genotype
accuracy in real samples.

## AIDA technical validation

The AIDA study evaluated Module 1 technical integration using
48 independently generated and verified evidence bundles.

See [AIDA technical validation](aida_technical_validation.md).

## Qdenga independent cohort application

UAT-001 evaluated the frozen v0.1.0 workflow on 21 held-out
Qdenga samples.

All 21 samples completed build, verify, and QC successfully.
A separate clean installation of the published wheel reproduced
one sample's evidence tables byte-identically and generated
identical aggregate QC.

See [Qdenga independent application](qdenga_independent_application.md).

## Outstanding validation

Independent genotype concordance and accuracy assessment
against an orthogonal truth set, such as matched WGS, remains
outstanding.

Technical allele MATCH, VCF_CALLED, and successful bundle
verification must not be interpreted as genotype accuracy.
