#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

inputs:
  alignment:
    type: File
    inputBinding:
      position: 2
      valueFrom: $(self.basename)
    label: Input bam file

baseCommand: [samtools, index, -b]

outputs:
  alignment_with_index:
    type: File
    secondaryFiles: .bai
    outputBinding:
      glob: $(inputs.alignment.basename)
    doc: The index file


