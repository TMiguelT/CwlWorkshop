class: Workflow
cwlVersion: v1.0
id: alignment
label: alignment
$namespaces:
  sbg: 'https://www.sevenbridges.com'
inputs:
  - id: reference
    type: File
    'sbg:x': -710.003662109375
    'sbg:y': -179.5
  - id: reads
    type: 'File[]'
    'sbg:x': -713.3142700195312
    'sbg:y': -26.5
outputs:
  - id: sorted_alignment
    outputSource:
      - samtools_index/sorted_alignment
    type: File
    'sbg:x': -34
    'sbg:y': -108
steps:
  - id: bwa_mem
    in:
      - id: output_filename
        default: alignment.bam
      - id: reads
        source:
          - reads
      - id: reference
        source: reference
    out:
      - id: output
    run: workflows/tools/bwa-mem.cwl
    'sbg:x': -570
    'sbg:y': -105
  - id: samtools_sort
    in:
      - id: alignment
        source: bwa_mem/output
    out:
      - id: sorted_bam
    run: ./samtools-sort.cwl
    label: samtools-sort
    'sbg:x': -383
    'sbg:y': -106
  - id: samtools_index
    in:
      - id: alignment
        source: samtools_sort/sorted_bam
    out:
      - id: sorted_alignment
    run: ./samtools-index.cwl
    label: samtools-index
    'sbg:x': -207
    'sbg:y': -108
requirements: []
