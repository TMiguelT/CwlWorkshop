class: Workflow
cwlVersion: v1.0

inputs:
  - id: reference
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
  - id: reads
    type: 'File[]'

outputs:
  variants:
    outputSource: freebayes/output
    type: File

steps:
  bwa:
    in:
      output_filename:
        valueFrom: 'alignment.bam'
      reads: reads
      reference: reference
    out:
      - alignment
    run: bwa-mem.cwl

  sort:
    in:
      alignment: bwa/alignment
    out:
      - sorted_alignment
    run: samtools-sort.cwl

  index:
    in:
      alignments: sort/sorted_alignment
    out:
      - alignments_with_index
    run: samtools-index.cwl

  freebayes:
    in:
      reference: reference
      bam: index/alignments_with_index
    out:
      - output
    run: freebayes.cwl
