# SpliceVariantRNA
---

**SpliceVariantRNA** is a Snakemake pipeline designed to study alternative splicing in RNA sequencing data.
It automates and enhances the SpliceLauncher tool by applying additional filtering steps.
For each detected and sorted splicing event, a Sashimi plot is generated to facilitate interpretation.

## Table of Contents

- [SpliceVariantRNA](#splicevariantrna)
- [Installation & Setup](#installation--setup)
    - [Managing the SpliceLauncher Submodule](#managing-the-splicelauncher-submodule)
    - [Cloning the repository with its submodule](#cloning-the-repository-with-its-submodule)
    - [Initializing the submodule (if not cloned with `--recursive`)](#initializing-the-submodule-if-not-cloned-with---recursive)
    - [Updating the submodule](#updating-the-submodule)
    - [Installing Dependencies](#installing-dependencies)
        - [Python dependencies](#python-dependencies)
        - [R dependencies](#r-dependencies)
- [Requirements](#requirements)
    - [Checking versions](#checking-versions)
    - [Configuration File (config_example.yaml)](#configuration-file-config_exampleyaml)
- [Input Files & Data Structure](#input-files--data-structure)
    - [Required Input Files](#required-input-files)
    - [Input Directory Structure](#input-directory-structure)
    - [Sample Metadata (samples_metadata_example.tsv)](#sample-metadata-samples_metadata_exampletsv)
    - [Genes of Interest (genes_of_interest_example.txt)](#genes-of-interest-genes_of_interest_exampletxt)
    - [Color Palette (palette.txt)](#color-palette-palettetxt)
- [Running the SpliceVariantRNA](#running-the-splicevariantrna)
    - [Dry Run (Simulation Mode)](#dry-run-simulation-mode)
    - [Running the pipeline](#running-the-pipeline)
    - [Generating Dependency Graphs](#generating-dependency-graphs)
- [Processed Data Structure](#processed-data-structure)
    - [Reference data for mapping (`references/`)](#reference-data-for-mapping-references)
    - [Filtered FASTQ files **if trimming is enabled** (`fastq_trimmed/`)](#filtered-fastq-files-if-trimming-is-enabled-fastq_trimmed)
    - [Read alignments (`BAM/`)](#read-alignments-bam)
    - [Quality Control Structure](#quality-control-structure)
- [Result Data Structure](#result-data-structure)
    - [Directory Structure](#directory-structure)
    - [SpliceLauncher Output](#splicelauncher-output)
    - [Statistical Junctions](#statistical-junctions)
        - [Sashimi Plots for Statistical Junctions](#sashimi-plots-for-statistical-junctions)
    - [Non Statistical Junctions](#non-statistical-junctions)
        - [Filtering Criteria](#filtering-criteria)
        - [Files Generated](#files-generated)
        - [Sashimi Plots for Non Statistical Junctions](#sashimi-plots-for-non-statistical-junctions)
    - [Recap file](#recap-file)
        - [HGVS Annotation](#hgvs-annotation)
- [Authors](#authors)
- [License](#license)


# Installation & Setup

## Managing the SpliceLauncher Submodule

A modified version of SpliceLauncher is integrated into SpliceVariantRNA as a **Git submodule**.
After cloning the main repository, **make sure to properly retrieve the submodule** to avoid an empty directory.

### Cloning the repository with its submodule
If you **haven’t cloned the repository yet**, use **this command** to retrieve **SpliceVariantRNA and its submodule** in one step:
```bash
git clone --recursive https://github.com/CMarcoBioinfo/SpliceVariantRNA.git
```

### Initializing the submodule (if not cloned with `--recursive`)
If you **have already cloned the repository but did not use `--recursive`**, you must manually initialize the submodule:
```bash
git submodule update --init --recursive
```

### Updating the submodule
If updates have been made in **SpliceLauncher** and you want to retrieve them, run:
```bash
git submodule update --remote --recursive
```
This ensures your version is synchronized with the latest submodule update from the remote repository.
Without this step, your submodule may still be pointing to an older commit instead of the latest version!

---

## Installing Dependencies
To ensure that all required **Python and R libraries** are correctly installed, use the provided requirement files.

### Python dependencies
The file **`requirements.txt`** contains all necessary Python packages.
Using a virtual environment is recommended.
To install them, run:
```bash
pip install -r requirements.txt
```

### R dependencies
The file **`packages.R`** contains all required R libraries.
To install them, run:
```bash
Rscript packages.R
```

This ensures that **all dependencies are installed** and properly configured for the pipeline.
Make sure to verify the installed versions with the following commands:
```bash
python --version
R --version
```

---

## Requirements
The following versions were used for this pipeline, but **other versions may also work**.
However, compatibility is **not guaranteed** for versions significantly different from those listed below.

Before you start, make sure you have installed:
- **Python 3.12**
- **R (4.1.2)**
- **Perl (v. 5.34.0)**
- **awk (mawk 1.3.4)**
- **Snakemake (8.25.5)**
- **pigz (2.6)**
- **fastp (0.20.1)**
- **FastQC (v. 0.11.9)**
- **MultiQC (v. 1.24.1)**
- **STAR (2.7.11b)**
- **samtools (1.13)**
- **bedtools (v2.30.0)**
- **ggsashimi (v1.1.5)**

### Checking versions
Run the following commands to check if everything is installed correctly:
```bash
python --version
R --version
snakemake --version
perl -v
awk -W version
pigz -V
fastp --version
fastqc --version
multiqc --version
STAR --version
samtools --version
bedtools --version
ggsashimi --version
```

### Configuration File ([config_example.yaml](config_example.yaml))
Before running the pipeline, make sure to **properly configure the [config_example.yaml](config_example.yaml)
file** according to your environment.
This file allows you to **define paths, parameters, and module activations** for the workflow.

## Input Files & Data Structure
Before running the pipeline, ensure that your **input files are available** and that their paths are correctly defined in the configuration file [config_example.yaml](config_example.yaml).

### Required Input Files
- **Sample metadata ([samples_metadata_example.tsv](/data/1-raw_data/metadata/samples_path/samples_metadata_example.tsv))**: This file must be filled in by the user (including sample identifiers, file paths, sequencing date, and technology platform.)
- **Genome reference (`{your_reference}.fasta`)**: FASTA file of the reference genome used for mapping.
- **Annotation files (`{your_annotation}.gff3`)**: This file corresponds to the GFF3 annotation of the reference genome and is crucial for splicing analysis.
- **Raw RNA-seq files (`fastq/`)** → FASTQ files containing sequencing data (`.fastq` or `.fastq.gz`).
- **MANE annotation ([MANE.to.SpliceLauncher.txt](data/1-raw_data/annotations/MANE.to.SpliceLauncher.txt/MANE.to.SpliceLauncher.txt))**: **A default MANE file is provided**, but users can use their own if needed.


### Input Directory Structure
**The pipeline does NOT impose a strict directory structure** for input files. However, it's best to follow this structure.
Paths to each file must be specified in [config_example.yaml](config_example.yaml), which allows users to **organize their data freely**.

Example of a common directory structure:
```plaintext
data/
└── 1-raw_data/
    ├── annotations/
    │   ├── MANE.to.SpliceLauncher.txt
    │   ├── {your_annotation}.gff3
    │   ├── ({your_annotation}.gtf)
    ├── fastq/
    │   ├── sample_001_R1.fastq.gz
    │   ├── sample_001_R2.fastq.gz
    │   ├── sample_002_R1.fastq.gz
    │   ├── sample_002_R2.fastq.gz
    │   ├── sample_003_R1.fastq.gz
    │   ├── sample_003_R2.fastq.gz
    ├── metadata/
    │   ├── genes_of_interest/
    │   │   └── (genes_of_interest_example.txt)
    │   ├── others/
    │   │   └── palette.txt
    │   └── samples_path/
    │       └── samples_metadata_example.tsv
    └── references/
        └── {your_reference}.fasta
```

### Sample Metadata ([samples_metadata_example.tsv](data/1-raw_data/metadata/samples_path/samples_metadata_example.tsv))
The **[samples_metadata_example.tsv](data/1-raw_data/metadata/samples_path/samples_metadata_example.tsv)** file is essential for the pipeline, as it provides metadata about each sequencing sample.
Users **must fill in this file** before running the pipeline to ensure proper handling of input data.

#### File Format
The file follows a **tab-separated format (`.tsv`)**, with the following **columns**:

| Column Name      | Description |
|-----------------|-------------|
| `id`            | **Sample identifier** (unique name for each sample) |
| `path_read1`    | **Path to the forward (R1) read file** |
| `path_read2`    | **Path to the reverse (R2) read file** |
| `sequencing_date` | **Date of sequencing** in YYYYMMDD format or others |
| `technology`   | **Sequencing technology** (e.g., Illumina, Nanopore) |

#### Example
Here is an example of a correctly formatted file:

```plaintext
id  path_read1  path_read2  sequencing_date technology
SAMPLE_001  /data/fastq/sample_001_R1.fastq.gz  /data/fastq/sample_001_R2.fastq.gz  20240510    Illumina
SAMPLE_002  /data/fastq/sample_002_R1.fastq.gz  /data/fastq/sample_002_R2.fastq.gz  20240510    Nanopore
SAMPLE_003  /data/fastq/sample_003_R1.fastq.gz  /data/fastq/sample_003_R2.fastq.gz  20240510    Illumina
```

### Genes of Interest ([genes_of_interest_example.txt](data/1-raw_data/metadata/genes_of_interest/genes_of_interest_example.txt))
**This file lists specific genes that should be retained in the final results.**
By specifying these genes, the pipeline ensures that only relevant splicing events related to them are highlighted.
If no gene list is provided, the pipeline analyzes all detected splicing events without filtering.

#### Example of a genes of interest file
```plaintext
KIF5A
LMNB1
MAPT
MATR3
NEK1
OPTN
```

**Users can modify this file to include other genes of interest** based on their study goals.
If no file is provided, the pipeline will analyze all detected splicing events without filtering.

### Color Palette ([palette.txt](data/1-raw_data/metadata/others/palette.txt))
**This file controls the color for Sashimi plots** and highlights the studied sample compared to others randomly selected in the visualization.
**Users can modify this file to apply custom colors** and enhance visibility in the plots.
If no file is provided, the default settings will be used for visualization.

## Running the SpliceVariantRNA

Before running the pipeline, ensure that you have correctly followed the previous steps.
By default, all output files are written under the data/ directory in structured subfolders.

### Dry Run (Simulation Mode)
Before starting the pipeline, we recommend that you check that everything is properly configured:
```bash
snakemake --configfile {your_config.yaml} --cores 1 --dry-run
```
- `--configfile {your_config.yaml}` : Specifies the configuration file.
- `--dry-run` : Runs the pipeline in simulation mode.

### Running the pipeline
Once everything is set up, launch the pipeline with:
```bash
snakemake --configfile {your_config.yaml} -p -c {N}
``` 
- `-p` : Print out the shell commands that will be executed.
- `-c {N}` : Use at most **N CPU cores/jobs** in parallel. Replace `N` with the number of available CPU cores.

### Generating Dependency Graphs
Snakemake allows you to visualize the structure and dependencies of the workflow using different graph types. These graphs help you understand execution order and rule dependencies.

```bash
snakemake --configfile {your_config.yaml} --dag | dot -Tpng > dag.png
snakemake --configfile {your_config.yaml} --rulegraph | dot -Tpng > rulegraph.png
snakemake --configfile {your_config.yaml} --filegraph | dot -Tpng > filegraph.png
``` 
- `--dag` : Generates a graph showing the dependencies between different rules in the pipeline.
- `--rulegraph` : Displays only the relationships between the workflow rules, without showing files.
- `--filegraph` : Visualizes dependencies between input and output files used in the pipeline.

# Processed Data Structure  

Once the pipeline is executed, the generated files are organized into three major categories: 

## Reference data for mapping (`references/`) 
This directory contains **reference files generated by `SpliceLauncherDB.r`**, a script from SpliceLauncher used in the pipeline for sequence alignment and splicing analysis.

```bash
data/
└── 2-processed_data/
    └── references/
        └── {your_genome}/
            ├── {your_genome}
            ├── BEDannotation.bed
            ├── exonInfo.tab
            ├── geneInfo.tab
            ├── SJDBannotation.sjdb
            ├── SpliceLauncherAnnot.txt
            ├── transcriptInfo.tab
            └── other reference files...
```

---

## Filtered FASTQ files **if trimming is enabled** (`fastq_trimmed/`)  
This directory contains **trimmed FASTQ files**. 
- **Only created if trimming is enabled (`TRIMMING: 1`)**  
- **Reads shorter than `{N}bp` are removed** (configurable length threshold).  

```bash
data/
└── 2-processed_data/
    └── fastq_trimmed/
        ├── sample_001.{N}bp.1.fastq.gz  # Reads after filtering
        ├── sample_001.{N}bp.2.fastq.gz
        ├── sample_002.{N}bp.1.fastq.gz
        ├── sample_002.{N}bp.2.fastq.gz
        ├── sample_003.{N}bp.1.fastq.gz
        ├── sample_003.{N}bp.2.fastq.gz
        └── other filtered FASTQ files...
```
> ⚠ **If trimming is disabled (`TRIMMING: 0`), this directory is not created**, and files remain in `1-raw_data/`.

---

## Read alignments (`BAM/`) 
This directory contains **BAM alignment files generated by STAR**, along with detected splicing junctions analyzed by SpliceLauncher.  
- **Organized by reference genome (`Homo_sapiens.GRCh37.dna.primary_assembly.chr`)**  
- **Divided into two subdirectories:**  
  - **`mapping/`** : Stores BAM files aligned after STAR processing.  
  - **`SpliceLauncher/`** : Contains detected splicing junctions in BED format.  

```bash
data/
└── 2-processed_data/
    └── BAM/
        └── {your_reference}/
            ├── mapping/
            │   ├── sample_001[.{N}bp].sorted.bam  # STAR-aligned reads
            │   ├── sample_001[.{N}bp].sorted.bam.bai  # BAM index file
            │   ├── sample_002[.{N}bp].sorted.bam 
            │   ├── sample_002[.{N}bp].sorted.bam.bai  
            │   ├── sample_003[.{N}bp].sorted.bam
            │   ├── sample_003[.{N}bp].sorted.bam.bai
            │   └── log_star/  # STAR log files
            └── SpliceLauncher/
                ├── sample_001[.{N}bp]_juncs.bed  # Detected splice junctions
                ├── sample_002[.{N}bp]_juncs.bed
                ├── sample_003[.{N}bp]_juncs.bed
                └── getClosestExons/  # Closest exon information
                        ├── sample_001[.{N}bp].count
                        ├── sample_002[.{N}bp].count
                        └── sample_003[.{N}bp].count
```

> **Note**: In all filenames, parts enclosed in brackets (e.g., `[.{N}bp]`) indicate optional suffixes that appear only when trimming is enabled.


## Quality Control Structure

Once the pipeline execution is complete, the quality control files are stored under the `3-Quality_control/` directory.  
This directory contains quality control for raw FASTQ files, trimmed FASTQ files, and BAM files.
These pieces of information are merged using MultiQC, and for each type of data, a dedicated report is generated.

```bash
data/
└── 3-Quality_control/
    ├── fastqc_raw/         ┐     
    ├── fastp_trimming/     | Quality control directories
    ├── fastqc_trimmed/     |
    ├── BAM/                ┘
    └── multiqc/
        ├── fastq_raw/
        │   ├── {your_prefix}_{unique_id}_data/         ┐
        │   └── {your_prefix}_{unique_id}.html          |
        ├── fast_trimmed/                               |
        │   ├── {your_prefix}_{unique_id}_data/         |
        │   └── {your_prefix}_{unique_id}.html          | MultiQC report
        └── BAM/                                        |
            └── {your_reference}/                       |
                ├── {your_prefix}_{unique_id}_data/     |
                └── {your_prefix}_{unique_id}.html      ┘
```

---

# Result Data Structure  

Once the pipeline execution is complete, the results are stored under the `4-Results/` directory.  
Each analysis is grouped under a **prefix-based subdirectory** (`{your_prefix}/`), allowing users to differentiate between multiple runs.  

## Directory Structure  

```bash
data/
└── 4-Results/
    └── {your_prefix}/
        └── SpliceLauncher/
            ├── merged_file/
            │   └── {prefix}_{unique_id}.txt
            ├── {your_prefix}_{unique_id}_report_{date}.txt
            └── {your_prefix}_{unique_id}_results/
                ├── {your_prefix}_{unique_id}.bed
                ├── {your_prefix}_{unique_id}.[xlsx,txt]
                ├── {your_prefix}_{unique_id}.statistical_junctions.[xlsx,txt]
                ├── {your_prefix}_{unique_id}.non_statistical_junctions.[xlsx,txt]
                └── samples_results/
                    ├── sample_001[.{N}bp]/
                    │   ├── sample_001[.{N}bp].pdf
                    │   ├── sample_001[.{N}bp]_genes.pdf
                    │   ├── sample_001[.{N}bp]_statistical_junctions.pdf
                    │   ├── sample_001[.{N}bp]_statistical_junctions.[xlsx,txt]
                    │   ├── sample_001[.{N}bp]_statistical_junctions.filter.[xlsx,txt]
                    │   ├── sample_001[.{N}bp]_non_statistical_junctions.pdf
                    │   ├── sample_001[.{N}bp]_non_statistical_junctions.[xlsx,txt]
                    │   ├── sample_001[.{N}bp]_non_statistical_junctions.filter.[xlsx,txt]
                    │   ├── sample_001[.{N}bp].recap.[xlsx,txt]
                    │   └── sashimi_plot/
                    │       ├── statistical_junctions/
                    │       │   ├── {chr}_{start}_{stop}_{gene}.pdf
                    │       │   └── other sashimi plot ....
                    │       └── non_statistical_junctions/
                    │           ├── unique
                    │           │   ├── {chr}_{start}_{stop}_{gene}.pdf
                    │           │   └── other sashimi plot ....
                    │           ├── threshold_exceeded
                    │           │   ├── {chr}_{start}_{stop}_{gene}.pdf
                    │           │   └── other sashimi plot ....
                    │           └── no_model
                    │               ├── {chr}_{start}_{stop}_{gene}.pdf
                    │               └── other sashimi plot ....
                    │
                    ├── sample_002[.{N}bp]/
                    │   └── .....
                    └── sample_003[.{N}bp]/
                        └── .....
```

---

## SpliceLauncher Output  
This directory contains results from **SpliceLauncher**, including detected junctions, statistical analyses, and visualizations.
For more details, see the [SpliceLauncher documentation](https://github.com/LBGC-CFB/SpliceLauncher?tab=readme-ov-file).

- **`{prefix}_{unique_id}.txt`**: Merged summary of detected junctions.
- **`{your_prefix}_{unique_id}_report_{date}.txt`**: Summary report of the SpliceLauncher execution.
- **`{your_prefix}_{unique_id}.bed`**: BED file listing detected junctions.
- **`{your_prefix}_{unique_id}.[xlsx,txt]`**: Comprehensive report of all junctions detected by SpliceLauncher.

| Column names | Example | Description |
|------------:|:--------:|:------------:|
| Conca | chr13_32915333_32920963 | The junction id (chr_start_end) |
| chr | chr13 | Chromosome number |
| start | 32915333 | Genomic coordinate of start junction<br/>End if on reverse strand |
| end | 32920963 | Genomic coordinate of end junction<br/>Start if on reverse strand |
| strand | + | Strand of the junction ('+': forward;<br/> '-':reverse) |
| Strand_transcript | forward | Strand of transcript |
| NM | NM_000059 | The transcript id according RefSeq nomenclature |
| Gene | BRCA2 | Gene symbol |
| *Sample* | 2250 | Read count |
| *P_Sample* | 15.25659623 | % of relative expression |
| event_type | SkipEx | The nature of junction:<br/>Physio: Natural junction<br/>SkipEx: Exon skipping<br/>5AS: Donor splice site shift<br/>3AS: Acceptor splice site shift<br/>NoData: Unannotated junction |
| AnnotJuncs | ∆12 | The junction names |
| cStart | c.6841 | Transcriptomic start coordinate of the junction |
| cEnd | c.6938 | Transcriptomic end coordinate of the junction |
| mean_percent | 12.60242 | Average in % of relative expression across samples |
| read_mean | 2683.769231 | Average of read count across samples |
| nbSamp | 11 | Number of time that the junction has been seen in samples |
| DistribAjust | - | The Distribution of junction expression (Gamma/N.binom) |
| Significative | NO | If a sample shown an abnormal expression of the junction |
| filterInterpretation | Unique junction | If a sample had an abnormal expression: Aberrant junction<br/>For unmodelized junction, if max expression >1%: Unique junction<br/>If junction does not fit the statistical model (Gamma/N.binom): No Model<br>If any P_sample value > 1000 or is infinite: Percentage threshold exceeded|

- **`sample_001[.{N}bp].pdf`** : Overview of detected junctions for all genes.
- **`sample_001[.{N}bp]_genes.pdf`** : Overview of detected junctions for genes with a junction.

![sample_001[.{N}bp]_genes.pdf](/doc/images/exemple_sample_001.pdf.png)

###
---

## Statistical Junctions
SpliceLauncher applies **statistical models** (`gamma` or `negative binomial`) to determine **which junctions are significantly expressed**.
Only junctions observed in **at least 5 samples** can undergo statistical testing.

### Filtering Criteria
Statistical junctions are filtered using **`THRESHOLD_SIGNIFICANCE_LEVEL`** in [config_example.yaml](config_example.yaml):  
- **`0`: No filtering (keep all junctions)**
- **`1` → Keep junctions with *, **, or *** significance**
- **`2` → Keep junctions with ** or *** significance**
- **`3` → Keep only junctions with *** significance**
- **`MAX_SAMPLES`**: The **maximum number of samples** in which a junction is classified as significant (default: `2`).

### Files Generated  
- **`{your_prefix}_{unique_id}.statistical_junctions.[xlsx,txt]`**: Full list of statistically significant junctions for sample
- **`sample_001[.{N}bp]_statistical_junctions.filter.[xlsx,txt]`**: Filtered set based on `THRESHOLD_SIGNIFICANCE_LEVEL` and `MAX_SAMPLES`
For **significant junctions**, three extra columns appear:  

| Column name | Example | Description |
|------------:|:--------:|:------------:|
| **nbSignificantSamples** | 2 | Number of samples where the junction is classified as significant |
| **p_value** | 0.0023 | Statistical p-value of the junction per sample |
| **SignificanceLevel** | `**` | Interpretation of the p-value:<br/>`*` (moderate), `**` (high), `***` (very high) |

- **`sample_001[.{N}bp]_statistical_junctions.pdf`** 

![sample_001[.{N}bp]_statistical_junctions.pdf](/doc/images/exemple_sample_001_statistical_or_non_statistical_junction.pdf.png)

#### Sashimi Plots for Statistical Junctions  
Filtered significant junctions (`sample_001[.{N}bp]_statistical_junctions.filter.[xlsx,txt]`) are visualized using **sashimi plots**: 

```bash
└── sashimi_plot/
    └── statistical_junctions/
        └── {chr}_{start}_{stop}_{gene}.pdf  # Sashimi plot of filtered junctions.
```

![exemple_sample_001_sashimi_plot_statistical_junctions.pdf](/doc/images/exemple_sample_001_sashimi_plot.pdf.png)

> ⚠ In red, the sashimi plot for the studied sample; Its NM transcript is also highlighted in red.

---

## Non Statistical Junctions
Non-statistical junctions are splicing events that cannot be evaluated statistically. These include several categories:

- **No Model**: Junctions that do not conform to the statistical model (Gamma/N.binom)
- **Unique Junctions**: Junctions with unmodeled expression profiles but a maximum expression above 1%, potentially observed in only one or a few samples.
- **Threshold Exceeded**: Junctions for which at least one P_sample value exceeds 1000 or is infinite, indicating extreme or atypical expression.

Although these junctions do not reach statistical significance, they may still be biologically or clinically meaningful, particularly in the context of rare disease studies, where statistical power is often limited.

### Filtering Criteria  
Non Statistical junctions are filtered using **two parameters** from `config.yaml`:  

- **`MAX_SAMPLES`**: The **maximum number of samples** in which a junction is classified as unique (default: `2`).
If a junction is found in more samples than `MAX_SAMPLES`, it **is excluded** from the unique junctions list.
If it appears in `MAX_SAMPLES` or fewer samples, it **is retained for further investigation**.

- **`MIN_READS_SAMPLE`**: The **minimum number of reads required in a sample** for a junction to be considered valid (default: `10`).
Samples containing reads at the junction **must pass the `MIN_READS_SAMPLE` threshold** to be counted.
If a sample does not meet the required read count, it **is not considered valid** for the unique junction classification.

### Files Generated
- **`{your_prefix}_{unique_id}.non_statistical_junctions.[xlsx,txt]`**: List of unique junctions that did not meet statistical testing criteria.
- **`sample_001[.{N}bp]_non_statistical_junctions.filter.[xlsx,txt]`**: Filtered set, retaining only junctions that meet `MAX_SAMPLES` and `MIN_READS_SAMPLE` thresholds.

For **Non Statistical junctions**, two extra columns appear: 

| Column name | Example | Description |
|------------:|:--------:|:------------:|
| **nbSampleFilter** | `2` | Number of samples that passed the filter and contain reads at the junction |
| **SampleRead** | `sample_001 reads = 11; sample_002 reads = 23` | Samples that passed filtering, along with their respective read counts |

Some columns seen above are not included here.

- **`sample_001_non_statistical_junctions.pdf`** : Same plot as for significant junctions

#### Sashimi Plots for Non Statistical Junctions
Filtered unique junctions (`sample_001_non_statistical_junctions.filter.[xlsx,txt]`) are visualized with **sashimi plots**:

```bash
└── sashimi_plot/
    └── non_statistical_junctions/
        ├── unique
        │   ├── {chr}_{start}_{stop}_{gene}.pdf # Sashimi plot of filtered junctions.
        │   └── other sashimi plot ....
        ├── threshold_exceeded
        │   ├── {chr}_{start}_{stop}_{gene}.pdf
        │   └── other sashimi plot ....
        └── no_model
            ├── {chr}_{start}_{stop}_{gene}.pdf
            └── other sashimi plot ....  
```

**Sashimi plots for unique junctions are identical to significant junction plots.**

**Why Non Statistical Junctions Matter?**  
- In **rare diseases**, an individual may **be the only one exhibiting a pathogenic mutation** at a given junction.
- Even though **these junctions lack statistical validation**, they could **highlight disease-driving mutations** that would otherwise be overlooked.

## Recap file

The recap file is a sample-focused summary that brings together all important splicing junctions detected for a given individual.
It consolidates the filtered statistical and non-statistical junctions into a single, simplified file, retaining only the most relevant information to facilitate interpretation.

**`sample_001[.{N}bp].recap.[xlsx,txt]`**: Unified report containing key events extracted from:
- **`sample_001[.{N}bp]_statistical_junctions.filter.[xlsx,txt]`**
- **`sample_001[.{N}bp]_non_statistical_junctions.filter.[xlsx,txt]`**

Depending on the output format:

- In .xlsx, the file is organized into 4 dedicated sheets:
    - Statistical Junctions
    - Unique Junctions
    - Threshold Exceeded Junctions
    - No Model Junctions

- In .txt, the content is structured into clearly labeled sections.

#### HGVS Annotation

For each junction retained in the recap file, a proposed HGVS-like annotation is computed using the transcript information (NM), the genomic coordinates, and the junction type. This annotation follows the HGVS recommendations for RNA splicing, with the following logic:

The NM transcript ID is checked and, if needed, updated using the MANE reference.

The junction symbol is interpreted:
- ∆ (but not ▼) → interpreted as a deletion
- ▼ (but not ∆) → interpreted as an insertion

A synthetic HGVS string is then generated:

- For deletions: NM_000059:r.6841_6938del
- For insertions: NM_000059:r.6841_6938insAUGC
- For other complex or ambiguous events: "Complex annotation"

> ⚠ This field is automatically generated to support downstream interpretation, but it is provided as a preliminary suggestion only and may require manual curation in clinical contexts.

## Authors

- **Corentin Marco** - [CMarcoBioinfo](https://github.com/CMarcoBioinfo)
    - You can contact me at: corentin.marco@chu-nimes.fr or corentin.marco.bioinfo@gmail.com

## License

This project is licensed under the MIT license - see the [LICENSE](LICENSE) file for details