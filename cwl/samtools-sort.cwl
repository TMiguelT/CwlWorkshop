class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: samtools_sort
baseCommand:
  - samtools
  - sort
inputs:
  - id: alignment
    type: File
    inputBinding:
      position: 0
outputs:
  - id: sorted_alignment
    type: File
    outputBinding:
      glob: sorted_alignment.bam
label: samtools-sort
arguments:
  - position: 0
stdout: sorted_output.bam
