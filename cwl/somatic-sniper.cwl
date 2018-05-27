cwlVersion: v1.0
class: CommandLineTool

baseCommand: bam-somaticsniper

arguments:
    - position: 2
      valueFrom: output.vcf

inputs:
    reference:
        inputBinding:
            prefix: -f
        type: File
        label: Reference genome

    tumour:
        inputBinding:
            position: 0
        type: File
        label: Tumour BAM

    normal:
        inputBinding:
            position: 1
        type: File
        label: Normal BAM

outputs:
    vcf:
        type: stdout

stdout: variants.vcf
