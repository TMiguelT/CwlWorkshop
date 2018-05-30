class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: freebayes
baseCommand:
  - freebayes
inputs:
  - id: reference
    type: File
    inputBinding:
      position: 0
      prefix: '--fasta-reference'
    secondaryFiles:
      - .fai
  - id: bam
    type: File
    inputBinding:
      position: 0
    secondaryFiles:
      - .bai
outputs:
  - id: variants
    type: File
    outputBinding:
      glob: variants.vcf
label: freebayes
requirements:
  - class: DockerRequirement
    dockerPull: maxulysse/freebayes
stdout: variants.vcf
