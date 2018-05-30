class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: samtools_index
baseCommand:
  - samtools
  - index
inputs:
  - id: alignment
    type: File
    inputBinding:
      position: 0
outputs:
  - id: sorted_alignment
    type: File
    outputBinding:
      glob: $(inputs.alignment.basename)
    secondaryFiles:
      - .bai
label: samtools-index
requirements:
  - class: DockerRequirement
    dockerPull: biocontainers/samtools
  - class: InitialWorkDirRequirement
    listing:
      - entryname: $(inputs.alignment.basename)
  - class: InlineJavascriptRequirement
