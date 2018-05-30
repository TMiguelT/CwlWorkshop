class: Workflow

cwlVersion: v1.0

inputs:
  tumour_reads:
    type: 'File[]'
  normal_reads:
    type: 'File[]'
  reference:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa

outputs: []

steps:
  - id: normal_alignment
    in:
      reference: reference
      reads: normal_reads
    out:
      - alignment_with_index
    run: alignment_workflow.cwl

  - id: tumour_alignment
    in:
      reference: reference
      reads: tumour_reads
    out:
      - alignment_with_index
    run: alignment_workflow.cwl

  - id: somatic_sniper
    in:
      normal: normal_alignment/alignment_with_index
      reference: reference
      tumour: tumour_alignment/alignment_with_index
    out:
      - vcf
    run: somatic-sniper.cwl

requirements:
  - class: SubworkflowFeatureRequirement
