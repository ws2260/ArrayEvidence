# ArrayEvidence Technical Evidence Contract

Status: Implemented Module 1 evidence schema and TSV export
Updated: 2026-10-21
Schema version: `1.0.0`

## Scope

This document defines the current Module 1 technical evidence
contract, including internal reconciliation results and the
versioned TSV export.

ArrayEvidence preserves coordinate evidence, allele identity,
source-record provenance, and genotype projection as distinct
dimensions.

The export is technical evidence, not a replacement for primary
genotype calling, a clinical interpretation, or an anonymization
mechanism.

## Internal reconciliation evidence

The internal `ReconciliationResult` retains the following
coordinate-based fields:

- `coordinate`: manifest-derived genomic coordinate, or `None`.
- `status`: legacy coordinate-based reconciliation status.
- `candidate_record_count`: number of VCF records at that coordinate.
- `coordinate_evidence`: coordinate lookup outcome.
- `genotype_evidence`: call status of coordinate-matched records,
  not necessarily the record used for normalized allele matching.

These fields do not establish allele equivalence.

In particular, `VCF_ABSENT` means no record was found at the
manifest coordinate. It does not establish homozygous-reference
genotype, lack of assay coverage, or biological absence.

Allele-based fields include:

- `allele_evidence`: outcome of allele identity assessment.
- `allele_reason`: reason associated with that assessment.
- `indel_evidence`: normalized indel candidate-discovery outcome.
- `matched_alt_indices`: 1-based VCF ALT indices associated with
  matched allele candidates.
- `allele_sources`: original VCF record and 1-based ALT index for
  each retained allele candidate.

Indel candidate discovery operates on reference-normalized events,
not solely on manifest-to-VCF coordinate equality.

A shifted representation can therefore produce both
`coordinate_evidence=ABSENT` and `indel_evidence=UNIQUE`.

For indels, `UNIQUE` refers to one matching VCF ALT allele candidate,
not merely one VCF record. Multiple equivalent ALT entries within
a single VCF record remain ambiguous.

`indel_evidence=NOT_EVALUATED` does not mean that an indel was
searched and found absent. It can also occur when reference-based
indel processing was not performed.

## Genotype projection

The internal and exported projection fields are:

- `projection_status`: result of projecting the genotype onto the
  established assay-matched ALT allele.
- `projection_reason`: reason for that projection result.
- `assay_alt_dosage`: number of matched ALT copies when evaluable.
- `other_alt_present`: whether another nonreference ALT is present.
- `ploidy`: number of genotype alleles when evaluable.

A unique allele match does not imply a nonreference genotype.

`ASSAY_REF_ONLY` means the matched assay ALT was not present in
an evaluable genotype. It must not be inferred from a missing
VCF record or an unresolved allele match.

An ambiguous or unresolved indel candidate must not produce
assay ALT dosage.

`OTHER_ALT_PRESENT` takes precedence over `ASSAY_ALT_PRESENT`
when both the matched ALT and another ALT occur in the genotype.
The matched ALT dosage is still retained when evaluable.

Projection statuses supported by the current schema are:

- `ASSAY_ALT_PRESENT`
- `ASSAY_REF_ONLY`
- `OTHER_ALT_PRESENT`
- `NO_CALL`
- `NOT_EVALUATED`
- `UNRESOLVED`

For evaluable projections, dosage, ploidy, and the other-ALT flag
must be present and internally consistent. For nonevaluable
projections, those three fields must be null.

## Versioned external schema

Schema version `1.0.0` exports two tab-separated files:

- `assay_evidence.tsv`: one technical evidence row per sample-assay key.
- `allele_candidates.tsv`: zero or more source VCF allele candidates
  associated with each sample-assay key.

The version is included as `schema_version` in both tables.

### assay_evidence.tsv

The fixed column order is:

| Column | Meaning |
| --- | --- |
| `schema_version` | Evidence schema version |
| `sample_id` | Caller-supplied sample identifier |
| `assay_id` | Manifest assay identifier |
| `manifest_chrom` | Manifest chromosome, if mapped |
| `manifest_pos` | Manifest 1-based position, if mapped |
| `coordinate_status` | Legacy coordinate reconciliation status |
| `coordinate_evidence` | Coordinate lookup evidence |
| `coordinate_candidate_count` | VCF record count at manifest coordinate |
| `coordinate_genotype_evidence` | Coordinate-matched genotype evidence |
| `allele_evidence` | Reference-aware allele matching result |
| `allele_reason` | Reason for allele matching result |
| `indel_evidence` | Reference-normalized indel discovery result |
| `projection_status` | Assay-specific genotype projection status |
| `projection_reason` | Reason for projection result |
| `assay_alt_dosage` | Matched ALT copy count, if evaluable |
| `other_alt_present` | Presence of another nonreference ALT, if evaluable |
| `ploidy` | Number of genotype alleles, if evaluable |

The primary key is (`sample_id`, `assay_id`).

The `coordinate_candidate_count` field counts records at the
manifest coordinate, not reference-normalized allele candidates.

For example, a shifted indel can have
`coordinate_candidate_count=0`, `coordinate_status=VCF_ABSENT`,
and `indel_evidence=UNIQUE`, with a corresponding candidate row
at another VCF coordinate.

### allele_candidates.tsv

The fixed column order is:

| Column | Meaning |
| --- | --- |
| `schema_version` | Evidence schema version |
| `sample_id` | Parent sample identifier |
| `assay_id` | Parent assay identifier |
| `candidate_index` | Consecutive 1-based candidate index |
| `vcf_chrom` | Original VCF chromosome |
| `vcf_pos` | Original VCF 1-based position |
| `vcf_id` | Original VCF variant identifier |
| `vcf_ref` | Original VCF REF allele |
| `vcf_alt` | Matched original VCF ALT allele |
| `vcf_alt_index` | 1-based index in original VCF ALT list |
| `vcf_gt` | Original VCF genotype string |
| `candidate_status` | Candidate evidence status |

The candidate key is
(`sample_id`, `assay_id`, `candidate_index`).

Every candidate must reference an existing assay row. Candidate
indices must be unique and consecutive from 1 within each
sample-assay pair.

Candidate rows preserve source VCF representations. They are
not normalized VCF records and are not independent genotype calls.

### Parent-candidate consistency

The current exporter enforces these relationships:

| Parent evidence | Candidate requirement |
| --- | --- |
| Indel `UNIQUE` | Exactly one candidate, status `UNIQUE` |
| Indel `MULTIPLE` | At least two candidates, all status `MULTIPLE` |
| Indel `ABSENT` or `UNRESOLVED` | No candidates |
| Indel `NOT_EVALUATED`, allele `MATCH` | Exactly one candidate, status `MATCH` |
| Indel `NOT_EVALUATED`, allele `PARTIAL_MATCH` | Exactly one candidate, status `PARTIAL_MATCH` |
| Indel `NOT_EVALUATED`, other allele evidence | No candidates |

An indel with `MULTIPLE` evidence must retain unresolved allele
evidence and must not have an evaluable genotype projection.

These rules reflect the current reconciliation implementation.
They do not independently validate the underlying allele identity.

## TSV serialization contract

The exporter uses:

- UTF-8 encoding.
- Tab-delimited fields.
- Unix newline characters.
- Fixed column order from the versioned dataclasses.
- Deterministic assay sorting by sample ID and assay ID.
- Deterministic candidate sorting by sample ID, assay ID,
  and numeric candidate index.
- Complete serialized rows as deterministic tie-breakers.

`NULL` is the reserved token for missing values.

The following values remain distinct:

| Value | Serialized representation |
| --- | --- |
| Missing value (`None`) | `NULL` |
| Integer zero | `0` |
| Boolean false | `False` |
| Boolean true | `True` |
| Empty string | Empty TSV field |

Literal string `NULL` is rejected to avoid ambiguity.

Tabs, carriage returns, and newlines within field values are
rejected rather than silently escaped.

No genotype or allele inference is performed during serialization.

## Immutable evidence-bundle publication

The exporter interface is:

`export_evidence_tsv(output_dir, assay_rows, candidate_rows)`

It returns the paths to the two published TSV files.

The destination must not already exist, and its parent directory
must exist.

Publication proceeds as follows:

1. Materialize and validate both evidence tables.
2. Check schema values, relational integrity, and evidence consistency.
3. Serialize both TSVs in memory.
4. Create a temporary staging directory beside the destination.
5. Write both TSV files into staging.
6. Publish the completed directory using Linux `renameat2`
   with `RENAME_NOREPLACE`.
7. Clean up staging after a failed write or publication attempt.

An existing destination is never intentionally overwritten.
Concurrent destination creation is handled by the no-replace
publication operation.

On a supported local filesystem, publication provides atomic
visibility of the completed directory. It does not guarantee
power-loss durability.

The current implementation requires Linux/WSL support for
`renameat2` and the filesystem's no-replace rename operation.

The exporter does not implement an overwrite or update-in-place
operation for published bundles.

## Privacy and authorization

`sample_id` is supplied by the caller and preserved unchanged.

The exporter does not anonymize, pseudonymize, encrypt, redact,
authorize access, or enforce consent and permitted-use restrictions.

Even an apparently de-identified genotype evidence table may
contain sensitive participant-derived genetic information.

Participant-derived evidence bundles must remain in authorized
environments and must not be committed to the publication source
repository.

Synthetic tests must be constructed independently of real
participant records.

Public reference resources and proprietary assay resources must
be handled under their applicable permissions.

## Interpretation boundaries

Reconciliation does not replace primary genotype calling.

Technical allele equivalence does not establish genotype accuracy,
clinical significance, assay sensitivity, or suitability for a
particular downstream analysis.

Reference assembly compatibility, source permissions, and
authorized analysis scope remain caller and workflow requirements.

The exported schema does not yet include a full assay QC report,
reconciled VCF/BCF output, clinical annotation, or biomedical
interpretation.

## Versioning and future changes

The current external schema is `1.0.0`.

Changes to field definitions, column order, missing-value
conventions, or evidence semantics require explicit schema
version review.

Potential future additions include:

- Machine-readable bundle metadata and checksums.
- Reference assembly and resource provenance.
- Additional assayability and technical QC evidence.
- Reconciled VCF/BCF export under explicit policies.
- Optional purpose-scoped downstream module interfaces.

These are not implemented by the current TSV export milestone.
