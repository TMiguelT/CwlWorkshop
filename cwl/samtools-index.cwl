class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
baseCommand:
  - samtools
  - index
inputs:
  - id: alignment
    type: File
    inputBinding:
      position: 0
    label: Input bam file
outputs:
  - id: alignment_with_index
    doc: The index file
    type: File
    outputBinding:
      glob: $(inputs.alignment.basename)
    secondaryFiles:
      - .bai
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.alignment)
  - class: InlineJavascriptRequirement
