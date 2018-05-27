#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

inputs:

  reference:
    type: File
    secondaryFiles:
        - .amb
        - .ann
        - .bwt
        - .pac
        - .sa
    inputBinding:
      position: 0

  reads:
    type:
      type: array
      items: File
    inputBinding:
      position: 1

outputs:
  alignment:
    type: stdout

stdout: output.bam

baseCommand:
- bwa
- mem
