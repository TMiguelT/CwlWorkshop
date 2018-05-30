class: Workflow
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
inputs:
  - id: reference
    type: File
    'sbg:x': -454
    'sbg:y': -152
  - id: reads
    type: 'File[]'
    'sbg:x': -524
    'sbg:y': 24
outputs:
  - id: output
    outputSource:
      - freebayes/output
    type: File
    'sbg:x': 400.2967529296875
    'sbg:y': 33
steps:
  - id: bwa_mem
    in:
      - id: reads
        source:
          - reads
      - id: reference
        source: reference
    out:
      - id: alignment
    run: ./bwa-mem.cwl
    'sbg:x': -325.54144287109375
    'sbg:y': 6.75
  - id: samtools_index
    in:
      - id: alignment
        source: samtools_sort/sorted_alignment
    out:
      - id: alignment_with_index
    run: ./samtools-index.cwl
    'sbg:x': 98
    'sbg:y': 44
  - id: samtools_sort
    in:
      - id: alignment
        source: bwa_mem/alignment
    out:
      - id: sorted_alignment
    run: ./samtools-sort.cwl
    'sbg:x': -109.703125
    'sbg:y': 28
  - id: freebayes
    in:
      - id: reference
        source: reference
      - id: bam
        source: samtools_index/alignment_with_index
    out:
      - id: output
    run: ./freebayes.cwl
    'sbg:x': 270.296875
    'sbg:y': 41
requirements: []
