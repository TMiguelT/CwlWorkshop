class: CommandLineTool

cwlVersion: v1.0

baseCommand:
  - freebayes

inputs:
  - id: reference
    type: File
    inputBinding:
      position: 0
      prefix: '--fasta-reference'

  - id: bam
    type: File
    inputBinding:
      position: 0

outputs:
  - id: output
    type: File
    outputBinding:
      glob: output.vcf

stdout: output.vcf
