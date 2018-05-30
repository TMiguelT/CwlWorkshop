class: Workflow
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
inputs:
  - id: normal_reads
    type: 'File[]'
    'sbg:x': 0
    'sbg:y': 214
  - id: reference
    type: File
    'sbg:x': -194
    'sbg:y': 95
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
  - id: tumour_reads
    type: 'File[]'
    'sbg:x': 0
    'sbg:y': 0
outputs:
  - id: vcf
    outputSource:
      - somatic_sniper/vcf
    type: stdout
    'sbg:x': 540.7967529296875
    'sbg:y': 121
steps:
  - id: normal_alignment
    in:
      - id: reads
        source:
          - normal_reads
      - id: reference
        source: reference
    out:
      - id: alignment_with_index
    run: alignment_workflow.cwl
    'sbg:x': 177
    'sbg:y': 193
  - id: tumour_alignment
    in:
      - id: reads
        source:
          - tumour_reads
      - id: reference
        source: reference
    out:
      - id: alignment_with_index
    run: alignment_workflow.cwl
    'sbg:x': 173.359375
    'sbg:y': 19.5
  - id: somatic_sniper
    in:
      - id: normal
        source: normal_alignment/alignment_with_index
      - id: reference
        source: reference
      - id: tumour
        source: tumour_alignment/alignment_with_index
    out:
      - id: vcf
    run: ./somatic-sniper.cwl
    'sbg:x': 397.796875
    'sbg:y': 114
requirements:
  - class: SubworkflowFeatureRequirement
