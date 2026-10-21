# ArrayEvidence Technical Evidence Contract

Status: Initial implementation contract
Updated: 2026-10-21

## Scope

This document describes the current Module 1 reconciliation result.

The implementation preserves coordinate evidence, allele identity,
and genotype projection as distinct dimensions.

This is not yet a versioned external export schema.

## Coordinate-based fields

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

## Allele-based fields

- `allele_evidence`: outcome of allele identity assessment.
- `allele_reason`: reason associated with that assessment.
- `indel_evidence`: normalized indel candidate-discovery outcome.
- `matched_alt_indices`: 1-based VCF ALT indices associated with
  matched allele candidates.

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

## Genotype projection fields

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

## Interpretation boundaries

Reconciliation does not replace primary genotype calling.

Technical allele equivalence does not establish genotype accuracy,
clinical significance, assay sensitivity, or suitability for a
particular downstream analysis.

Reference assembly compatibility, source permissions, and
authorized analysis scope remain caller and workflow requirements.

## Future schema work

A future versioned export schema should make coordinate-based
and allele-based evidence explicit, preserve source-record
provenance, and distinguish not evaluated, unresolved, absent,
and no-call states without collapsing them.
