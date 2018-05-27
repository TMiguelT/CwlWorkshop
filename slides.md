layout: true
class: content
---
# Part 1: Introduction
.center[
![](images/cwl.png)
]
---
## Motivation

* Most bioinformatics involves running many command-line tools; aligners like `bwa` variant callers like `gatk`, and RNA
Seq tools like `cuffdiff`
* However, once this list of tools reaches a certain quantity and complexity, it becomes hard to  reproduce exactly what you ran and with what parameters
* A bash script may help with this, but proper workflows...

    * Run fully parallel to speed up execution
    * Work automatically with batch systems like SLURM
    * Are written declaratively, allowing the system to work out the optimal order of execution for you
    * Save you having to hard-code input parameters and temporary files
    * Are much more readable than bash
---
## CWL Structure
* Tools
    * Are wrappers that describe to the CWL engine how a tool works
    * Include its inputs and outputs, and their format
    * Any given tool can have a "correct" tool definition, unlike a workflow
    * Some already exist for commonly used tools
* Workflows
    * Explain how tools are connected to each other and in what order
    * Are generally project-specific
    * Can be nested inside each other
---
## Tooling for CWL
* Rabix
    * <http://rabix.io/>
    * Is an open source toolkit for creating and running CWL
    * Its most useful tool is the Rabix Composer - an application for writing CWL tools and workflows graphically
* Toil
    * <http://toil.ucsc-cgl.org/>
    * CWL has a number of "executors" - applications that can actually run CWL
    * Toil is the best supported executor - it can run workflows on your laptop, on a cluster, on NeCTAR, on AWS, or an a number of other platforms
    * To install Toil for CWL, run `pip install toil[cwl]`
---
## Workshop Goal
* By the end of the workshop we intend to have a fully-functioning somatic variant calling pipeline
* This means, it should take sequencing reads, as if from a cancer patient, and determine the DNA mutations that have
    occurred in their tumour
* Most of the exercises will contribute to this goal, so don't delete what you've written!

---
# Part 2: Tools
.center[
![](images/tool.png)
]
---
## Obtaining Tool Definitions

There are a few useful sources of CWL tool definitions:
* Dockstore
    * <https://dockstore.org>
    * Dockstore - a database of CWL and WDL workflows and tools
    * Once you find a tool definition you like, click "Files" → "Descriptor Files" → Download
.center[
![](images/dockstore_circled.png)
]
---
## Obtaining Tool Definitions

There are a few useful sources of CWL tool definitions:
* Official CWL Workflows Repository
    * <https://github.com/common-workflow-language/workflows>
    * Once you find a tool definition you like, right-click on the "Raw" button and click "Save link as" section
.center[
![](images/workflow_repo.png)
]
---
## Obtaining Tool Definitions
.alert.alert-primary[
.alert-heading[
### Exercise
]
* Try to find a simple wrapper for the tool `bwa mem` from the
* Download that tool definition, and run it on the provided data with `toil-cwl-runner bwa.cwl` (it should prompt you
    on how to specify the input files)
* Save this tool definition - we'll use it in our pipeline later
]
* This didn't do much beyond just running the `bwa` tool
* However, it did ensure the tool ran in a Docker container with BWA installed
---
## More on CWL Tools

* At minimum, a CWL tool definition must have three things
    * A command (to to run
    * A list of inputs (command line arguments and stdin)
    * A list of outputs (files and stdout)
* We will investigate how to make these tool definitions first using Rabix, and then from scratch
---
# Wrapping Samtools
.alert.alert-primary[
.alert-heading[
### Exercise
]
Follow along with the instructions to make a tool wrapper for `samtools sort`
]
---
## Exercise - Wrapping BWA
1\. Start by making a new tool definition in Rabix

.center[
![](images/rabix_new_tool.png)
]
---
## Exercise - Wrapping BWA
2\. Name it after the tool you're wrapping

.center[
![](images/rabix_samtools_name.png)
]
---
## Exercise - Wrapping BWA
3\. Add the "base command" - the fixed part of the command that will never change

.center[
![](images/rabix_samtools_base.png)
]
---
## Exercise - Wrapping BWA

4\. Define the inputs(s)

.center[
![](images/rabix_samtools_inputs.png)
]
---
## Exercise - Wrapping BWA

5\. Define the output(s)

.center[
![](images/rabix_samtools_outputs.png)
]
---
## Exercise - Wrapping BWA

6\. If the command produces output from stdout, you must specify that in the "Other" section

.center[
![](images/rabix_samtools_other.png)
]
---
## Wrapping Freebayes
.alert.alert-primary[
.alert-heading[
### Exercise
]
* [Freebayes](https://github.com/ekg/freebayes) is a variant caller, which takes a BAM alignment and determines how this alignment differs from the "normal"
reference genome
* Freebayes is called on the command-line as follows:
    ```bash
    freebayes --fasta-reference h.sapiens.fasta NA20504.bam
    ```
* Write a new CWL tool wrapper for Freebayes that supports this command
]

---
## Docker
* We have given CWL instructions on how to *run* these tools, but not how to *get* these tools
* For this we can use Docker
* Docker images are tiny virtual machines that have applications pre-installed inside of them
* You can find docker images of many common bioinformatics tools in [Biocontainers](https://biocontainers.pro/registry/)
* Once you've found a Docker image, you can plug it into the "Docker Image" section in Rabix:

    ![](images/docker_container_section.png)

---
## Updating our tool to use Docker
.alert.alert-primary[
.alert-heading[
### Exercise
]
* Find an appropriate Docker image for `bwa`, `samtools sort` and `freebayes`, using Biocontainers
* Once you have found the right images, plug them into the "Docker Image" section
]

---
## Secondary Files
* Some files, like indexes, are never considered a main input file, but are instead designed to accompany another file,
for example `.bai` files which accompany `bam` alignments, and `.tbi` indices which accompany `vcf` variant calls.
* These are called secondary files:

.center[
![](images/secondary_file.png)
]
---
## JavaScript Expressions
* Sometimes, some of the values in our CWL need to be calculated dynamically
* For example, a command might create an output file whose name is based on the input file
* For example, `gzip file.txt` produces `file.txt.gz`
* In order to do this, we can embed a JavaScript expression in some CWL fields

---
## JavaScript Expressions
### Syntax
* `$(...)` denotes a single JavaScript expression that will be evaluated, and the result used for this field
    * e.g.
        ```javascript
        $(1 + 2)
        ```
* `${...}` denotes a scope in which functions, variables, loops etc may be used, and a `return` statement is needed
    to return the final value
    * e.g.
        ```javascript
        ${
            var x = 1;
            var y = 2;
            return x + y;
        }
        ```
* Refer to the [expressions](https://www.commonwl.org/v1.0/CommandLineTool.html#Expressions) section of the CWL spec
---
## JavaScript Expressions
### Variables
* Within the scope of a JavaScript expression, the following variables are available:
    * `inputs`: a dictionary of all the inputs provided to this tool. If these inputs have the type
        [`File`](https://www.commonwl.org/v1.0/CommandLineTool.html#File) or
        [`Dictionary`](https://www.commonwl.org/v1.0/CommandLineTool.html#Directory), they have special properties
    * `self`: value depends on the specific context of this expression. e.g. when used in `secondaryFiles`, `self`
        is set to the main input or output, so that the function can calculate the secondary files for this file
    * `runtime`:
        * `runtime.outdir`: an absolute path to the designated output directory
        * `runtime.tmpdir`: an absolute path to the designated temporary directory
        * `runtime.cores`: number of CPU cores reserved for the tool process
        * `runtime.ram`: amount of RAM in mebibytes (2**20) reserved for the tool process
        * `runtime.outdirSize`: reserved storage space available in the designated output directory
        * `runtime.tmpdirSize`: reserved storage space available in the designated temporary directory


---
## Wrapping Samtools Index

* Next, we need to wrap `samtools` - a utility for working with alignment files
.alert.alert-primary[
.alert-heading[
### Exercise
]
* Use what you have learned from wrapping `bwa` to make a wrapper for the `samtools index` subcommand
* You can find the samtools manual, including all command-line flags for `samtools index` here: <http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS>
* The output from `samtools index` will be a BAM file, with its `.bai` index as a secondary file
]
---
## Tool YAML

* CWL tools are written in a structure called `YAML`
* To view the raw YAML in the Rabix viewer, click the "Code" tab at the top
* In YAML, key-value pairs in a dictionary are indicated by a colon `:` and list elements are indicated by a dash `-`:
* The sections we previously worked on correspond to sections in the yaml file:
    * "Base Command" → `baseCommand`
    * "Output Ports" → `outputs`
    * "Input Ports" → `inputs`
---
## Tool YAML
.row[
.col-sm[
```yaml
class: CommandLineTool
cwlVersion: v1.0
id: bwa
baseCommand:
  - bwa
  - mem
inputs:
  - id: reference
    type: File
    inputBinding:
      position: 0
  - id: reads
    type: 'File[]'
    inputBinding:
      position: 1
  - id: read_group
    type: string?
    inputBinding:
      position: 0
      prefix: '-R'
```
]

.col-sm[
```yaml
outputs:
  - id: alignment
    type: File
    outputBinding:
      glob: alignment.bam
label: bwa
stdout: alignment.bam
```
]
]
---
## Tool YAML

* This is a tool (as opposed to a workflow)
    ```yaml
    class: CommandLineTool
    ```
* This follows version 1.0 of the CWL standard:
    ```yaml
    cwlVersion: v1.0
    ```
* This tool is called `bwa`:
    ```yaml
    id: bwa
    label: bwa
    ```
* The base command is `bwa mem`:
    ```yaml
    baseCommand:
      - bwa
      - mem
      ```
---

## Tool YAML
.row[
.col-sm[
* The inputs:
    ```yaml
    inputs:
      - id: reference
        type: File
        inputBinding:
          position: 0
      - id: reads
        type: 'File[]'
        inputBinding:
          position: 1
      - id: read_group
        type: string?
        inputBinding:
          position: 0
          prefix: '-R'
    ```
]
.col-sm[
* The outputs (including stdout):
    ```yaml
    outputs:
      - id: alignment
        type: File
        outputBinding:
          glob: alignment.bam
    stdout: alignment.bam
```
]
]
---
## Wrapping Somatic Sniper
.alert.alert-primary[
.alert-heading[
### Exercise
]
* Somatic Sniper is a (bad) somatic variant caller
* Its command line usage can be found [here](https://github.com/genome/somatic-sniper/blob/master/gmt/documentation.md#usage)
* Use what you've learned about YAML tool definitions to write a tool definition for `bam-somaticsniper`
* Note that we want VCF output!
]
---
# Part 3: Writing Workflows
![](images/workflow.svg)
---
## Workflows in the Rabix Composer

* Refresher: workflows define how your tools connect to each other to form a data flow
* In Rabix, you add tools to your workflow by dragging and dropping from the sidebar
.center[
![](images/drag_tool.png)
]
* To connect the tools, you then drag a line between input ports and output ports
.center[
<video autoplay loop>
    <source src="images/connect.mp4" type="video/mp4">
</video>
]
---
* To specify an input that the user must provide, drag it onto empty space
.center[
<video autoplay loop>
    <source src="images/rabix_workflow_inputs2.mp4" width="518" height="282" type="video/mp4">
</video>
]
* To specify an output that is not used by another tool, drag it onto empty space
.center[
<video autoplay loop>
    <source src="images/rabix_workflow_output2.mp4" width="518" height="282" type="video/mp4">
</video>
]
---
## Making a Variant Calling Pipeline in Rabix
.alert.alert-primary[
.alert-heading[
### Exercise
]
* Make a basic workflow that connects:
    * `bwa` → `samtools sort` → `samtools index` → `freebayes`
]
---
## Subworkflows
* Often when making workflows you will want to encapsulate a whole pipeline of tools and re-use that together
as a single unit, for example "alignment"
* To do this, simply make a workflow that encapsulates this functionality, and use as a step in another workflow
* In the Rabix composer, you can drag workflows onto a workflow in the same way you add tools, however they have this
logo:
.fas.fa-share-alt.fa-rotate-180[]
.center[
![](images/subworkflow.png)
]

---
# A Tumour-Normal Variant Caller
.alert.alert-primary[
.alert-heading[
### Exercise
]
* First, make an alignment workflow that connects:
    * `bwa` → `samtools sort` → `samtools index`
* Next, make a tumour-normal workflow that uses this workflow to align both the tumour and the normal reads, and then
feeds the alignments into Somatic Sniper
]

---
## Workflows in CWL Files
* Like CWL tools, workflows are represented in a YAML format
* TODO
---
.alert.alert-primary[
.alert-heading[
## Exercise - Manual Workflow
]
* Using YAML, re-implement the basic germline variant calling workflow
* As before, this should connect:
    * `bwa` → `samtools sort` → `samtools index` → `freebayes`
]
