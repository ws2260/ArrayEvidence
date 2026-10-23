# Qdenga Independent Application — ArrayEvidence v0.1.0

## Scope

Independent functional user acceptance testing (UAT-001) of the
frozen ArrayEvidence v0.1.0 release on 21 Qdenga cohort samples.

This is a technical application and reproducibility assessment,
not independent genotype accuracy or clinical validation.

## Workflow

Original Illumina GSA IDAT data were processed using DRAGEN Array
1.4.0 to generate single-sample GRCh38 VCFs. ArrayEvidence Module 1
was then run using its documented `build`, `verify`, and `qc` commands.

The assay resources were GSA-24v3-0_A2 BPM/CSV and
GSA-24v3-0_A1 EGT.

## Functional acceptance

| Criterion | Result |
|---|---|
| Samples | 21 |
| Successful build/verify/QC | 21/21 |
| Complete evidence bundles | 21/21 |
| Valid QC JSON outputs | 21/21 |
| Assay evidence rows | 13,734,567 |
| Allele candidate rows | 13,568,247 |
| Assay rows per sample | 654,027 |
| Candidate rows per sample | 646,107 |
| Median VCF_CALLED assay rows | 641,576 |
| Median processing time, 20 cohort runs | 323.63 seconds |

All evidence tables had consistent schemas across samples.

## Clean-install reproducibility

The published v0.1.0 wheel was downloaded and its SHA256 digest
verified against the GitHub release metadata.

A fresh virtual environment successfully installed the wheel.
Reprocessing one sample using the clean installation produced:

- Successful build, verify, and QC commands
- Byte-identical `assay_evidence.tsv`
- Byte-identical `allele_candidates.tsv`
- Identical aggregate QC JSON

Published wheel SHA256:
`67e700af922c950111d28ffb1f663c7d2854106e674a6ae404ce3e8bdbcc74da`

## Technical evidence distributions

The following allele evidence counts were identical across all
21 samples (654,027 assay rows per sample).

| Allele evidence | Rows per sample |
|---|---:|
| MATCH | 643,743 |
| PARTIAL_MATCH | 2,160 |
| MISMATCH | 10 |
| UNRESOLVED | 3,357 |
| NOT_EVALUATED | 4,757 |

The median number of VCF_CALLED assay rows was 641,576
(range 637,469–641,905). VCF_NO_CALL ranged from 762 to
5,198 assay rows.

All 10 MISMATCH cases per sample were classified as
ALLELE_SET_DIFFERENT. The 2,160 PARTIAL_MATCH cases were
associated with multiallelic VCF representation.

These classifications describe technical evidence, not
genotype accuracy or clinical validity.

## Release provenance

- Published release: v0.1.0
- Release tag commit: `dda3e5e`
- Documentation/UAT checkout: `dda3e5e`
- Published wheel SHA256 verified before clean installation
- Clean-install evidence tables reproduced byte-identically

The Qdenga cohort was held out from Module 1 development.
The assay resources and upstream DRAGEN workflow had previously
been established using AIDA. Accordingly, this is an independent
cohort application, not independent validation of every upstream
processing component.

## Interpretation and limitations

- Allele MATCH does not establish genotype concordance.
- VCF_ABSENT does not establish homozygous reference.
- Technical QC does not establish clinical validity.
- Clean-install reproducibility was tested on one sample.
- Independent genotype accuracy validation requires an orthogonal
  truth set, such as matched whole-genome sequencing.
- Results apply to the tested software, assay, and environment.

## Conclusion

UAT-001 met its prespecified functional acceptance criteria
without modifying the frozen ArrayEvidence v0.1.0 implementation.
