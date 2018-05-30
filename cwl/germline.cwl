class: Workflow
cwlVersion: v1.0
id: germline
label: germline
$namespaces:
  sbg: 'https://www.sevenbridges.com'
inputs:
  - id: reference
    type: File
    'sbg:x': -735
    'sbg:y': -238
  - id: reads
    type: 'File[]'
    'sbg:x': -817
    'sbg:y': -89
outputs:
  - id: variants
    outputSource:
      - freebayes/variants
    type: File
    'sbg:x': 234.93377685546875
    'sbg:y': -111.5
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
    'sbg:x': -549.203125
    'sbg:y': -113.5
  - id: samtools_sort
    in:
      - id: alignment
        source: bwa_mem/output
    out:
      - id: sorted_bam
    run: ./samtools-sort.cwl
    label: samtools-sort
    'sbg:x': -386
    'sbg:y': -85
  - id: freebayes
    in:
      - id: reference
        source: reference
      - id: bam
        source: samtools_index/sorted_alignment
    out:
      - id: variants
    run: ./freebayes.cwl
    label: freebayes
    'sbg:x': -27.066675186157227
    'sbg:y': -95.5
  - id: samtools_index
    in:
      - id: alignment
        source: samtools_sort/sorted_bam
    out:
      - id: sorted_alignment
    run: ./samtools-index.cwl
    label: samtools-index
    'sbg:x': -216.703125
    'sbg:y': -59.5
requirements: []
