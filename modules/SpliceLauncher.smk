# This rule extracts splice junctions from a sorted BAM file and converts them into a BED file. 
# It filters reads based on strand orientation, block count, and mapping quality, using samtools, bedtools, and custom SpliceLauncher scripts.

rule SpliceLauncher_create_bed:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam"),
        csi = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam.csi"),

    output:
        bed = f"{path_bam}{name_genome}/SpliceLauncher/{{group}}/{{reads}}_juncs.bed"

    params:
        samtools = samtools,
        bedtools = bedtools,
        SpliceLauncher = SpliceLauncher,
        perl = perl,
        awk = awk
    
    log:
        stderr = f"{working_directory}/logs/SpliceLauncher/SpliceLauncher_create_bed/{{group}}/{{reads}}.err"

    shell:
        "{params.samtools} view "
        "-b -f 0x40 -F 1024 {input.bam} | "
        "{params.bedtools} bamtobed -bed12 -i stdin | "
        "{params.awk} '{{if($10>1){{print $0}}}}' | "
        "{params.perl} {params.SpliceLauncher}/scripts/bedBlocks2IntronsCoords.pl y - | "
        "{params.awk} '{{if($5==255){{print $0}}}}' > {output.bed} 2> {log.stderr} && "
        "{params.samtools} view "
        "-b -f 0x80 -F 1024 {input.bam} | "
        "{params.bedtools} bamtobed -bed12 -i stdin | "
        "{params.awk} '{{if($10>1){{print $0}}}}' | "
        "{params.perl} {params.SpliceLauncher}/scripts/bedBlocks2IntronsCoords.pl n - | "
        "{params.awk} '{{if($5==255){{print $0}}}}' >> {output.bed} 2>> {log.stderr}"


# This rule processes junction information from a BED file and counts splice junctions by intersecting them with a reference annotation.
# It outputs a summary of counts for each junction, including genomic coordinates, closest exons, and read support.

rule SpliceLauncher_count_junctions:
    input:
        bed_ref = f"{working_directory}/2-processed_data/references/{name_genome}/BEDannotation.bed",
        bed = f"{path_bam}{name_genome}/SpliceLauncher/{{group}}/{{reads}}_juncs.bed"

    output:
        counts = f"{path_bam}{name_genome}/SpliceLauncher/getClosestExons/{{group}}/{{reads}}.count"

    params:
        bedtools = bedtools,
        awk = awk
        
    log:
        stderr = f"{working_directory}/logs/SpliceLauncher/SpliceLauncher_count_junctions/{{group}}/{{reads}}.err"

    shell:
        "sort -k1,1 -k2,2n {input.bed} | "
        "uniq -c | "
        "{params.awk} 'BEGIN{{OFS=\"\\t\"}}{{print $2,$3,$4,$1,$6,$7}}' | "
        "{params.bedtools} intersect -s -wa -wb -a stdin -b {input.bed_ref} | "
        "{params.awk} 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,$3,$6,$10,$4}}' "
        "> {output.counts} 2> {log.stderr}"


genes_of_interest = get_config_value(config["SPLICELAUNCHER"]["ANALYSE"].get("GENES_OF_INTEREST"), "", str)
if not os.path.isfile(genes_of_interest):
    genes_of_interest = f"{working_directory}/1-raw_data/metadata/genes_of_interest/{genes_of_interest}"

if not os.path.isfile(genes_of_interest):
    genes_of_interest = ""

else :
    print_once(f"Fichier {genes_of_interest} ...... OK")
    genes_of_interest = f"-g {genes_of_interest}"


# This rule merges splice junction count files from multiple samples into a single file.
# It uses a Perl script to aggregate data, optionally including genes of interest if specified.
# The output is a merged file, structured for downstream analysis, and stored in a dedicated directory.

rule SpliceLauncher_merge_count:
    input:
        counts = expand(f"{path_bam}{name_genome}/SpliceLauncher/getClosestExons/{{group}}/{{reads}}.count", zip, reads=all_samples, group=groups)

    output:
        merged_files = f"{path_results}/SpliceLauncher/merged_files/{prefix}_{unique_id}.txt"
    
    params:
        SpliceLauncher = SpliceLauncher,
        perl = perl,
        directory = directory(f"{path_results}/SpliceLauncher/merged_files/"),
        genes_of_interest = genes_of_interest
    
    log:
        stderr = f"{working_directory}/logs/SpliceLauncher/SpliceLauncher_merge_count/{prefix}_{unique_id}.err"

    shell:
        "{params.perl} {params.SpliceLauncher}/scripts/joinJuncFiles.pl -c {input.counts} {params.genes_of_interest} > {output.merged_files} 2> {log.stderr}"


sampleNames = result = '|'.join(all_samples)
nbIntervals = get_config_value(config["SPLICELAUNCHER"]["ANALYSE"].get("NB_INTERVALS"), int(10))
threshold = get_config_value(config["SPLICELAUNCHER"]["ANALYSE"].get("JUNCTION_DISPLAY_THRESHOLD"), int(1))
min_cov = get_config_value(config["SPLICELAUNCHER"]["ANALYSE"].get("MIN_COV"), int(5))
minNonStatisticalReads = get_config_value(config["SPLICELAUNCHER"]["POST_ANALYSE"]["NON_STATISTICAL_JUNCTIONS"].get("MIN_READS_SAMPLE"), int(10))

transcriptList = get_config_value(config["SPLICELAUNCHER"]["ANALYSE"].get("TRANSCRIPT_LIST"), "", str)
removeOther = get_config_value(config["SPLICELAUNCHER"]["ANALYSE"].get("REMOVE_OTHER"), "", str)
if transcriptList and os.path.isfile(transcriptList):
    transcriptList = f"--TranscriptList {transcriptList} "
    if removeOther:
        removeOther = "--removeOther "
    else:
        removeOther = ""
else:
    transcriptList = ""
    removeOther = ""

txt = get_config_value(config["SPLICELAUNCHER"]["ANALYSE"].get("TXT"), int(0))
if txt:
    txt = "--text "
    extension = ".txt"

else:
    txt = ""
    extension = ".xlsx"

#Format date and time as unique identifier
current_time = time.localtime()
date = time.strftime("%m-%d-%Y", current_time)
listOutputSpliceLauncher = []
count_report = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_report_{date}.txt"
listOutputSpliceLauncher.append(count_report)
outputSpliceLauncher = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/{prefix}_{unique_id}_outputSpliceLauncher{extension}"
listOutputSpliceLauncher.append(outputSpliceLauncher)

bedOut =  get_config_value(config["SPLICELAUNCHER"]["ANALYSE"].get("BED_OUT"), int(0))
if bedOut:
    bedOut = "--bedOut "
    bed = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/{prefix}_{unique_id}.bed"
    listOutputSpliceLauncher.append(bed)

else:
    bedOut = ""

Graphics =  get_config_value(config["SPLICELAUNCHER"]["ANALYSE"].get("GRAPHICS"), int(0))

if Graphics:
    Graphics = "--Graphics "
    pdf = expand(f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.pdf", reads=all_samples)
    pdf_statistical_genes = expand(f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.statistical_genes.pdf",  reads=all_samples)
    pdf_statistical_junctions = expand(f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.statistical_junctions.pdf", reads=all_samples)
    pdf_non_statistical_junctions = expand(f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.non_statistical_junctions.pdf", reads=all_samples)
    listOutputSpliceLauncher.append(pdf)
    listOutputSpliceLauncher.append(pdf_statistical_genes)
    listOutputSpliceLauncher.append(pdf_statistical_junctions)
    listOutputSpliceLauncher.append(pdf_non_statistical_junctions)


# This rule performs a comprehensive analysis of splice junctions using the SpliceLauncher tool.
# It aggregates count data from multiple samples, combines it with reference annotation, and outputs results in customizable formats (.txt or .xlsx).
# Additionally, it generates optional outputs such as BED files and graphical PDFs if specified in the configuration.
# Key parameters like intervals, thresholds, and minimum coverage can be adjusted for flexibility, while transcripts of interest and related filters are conditionally applied based on configuration.
# Outputs include reports, merged data, BED files, and visualizations stored in an organized directory structure.

rule SpliceLauncher_Analyse:
    input:
        merged_files = f"{path_results}/SpliceLauncher/merged_files/{prefix}_{unique_id}.txt",
        annot = f"{working_directory}/2-processed_data/references/{name_genome}/SpliceLauncherAnnot.txt"

    output:
        listOutputSpliceLauncher

    params:
        SpliceLauncher = SpliceLauncher,
        Rscript = Rscript,
        nbIntervals = nbIntervals,
        threshold = threshold,
        min_cov = min_cov,
        minNonStatisticalReads = minNonStatisticalReads,
        transcriptList = transcriptList,
        removeOther = removeOther,
        txt = txt,
        bedOut = bedOut,
        Graphics = Graphics,
        directory = f"{path_results}/SpliceLauncher/",
        sampleNames = sampleNames

    log:
        stdout = f"{working_directory}/logs/SpliceLauncher/SpliceLauncher_Analyse/{prefix}_{unique_id}.out",
        stderr = f"{working_directory}/logs/SpliceLauncher/SpliceLauncher_Analyse/{prefix}_{unique_id}.err"

    shell:
        "{params.Rscript} {params.SpliceLauncher}/scripts/SpliceLauncherAnalyse.r "
        "--input {input.merged_files} "
        "-O {params.directory} "
        "--RefSeqAnnot {input.annot} "
        "-n {params.nbIntervals} "
        "{params.transcriptList}"
        "{params.removeOther}"
        "{params.txt}"
        "{params.bedOut}"
        "{params.Graphics}"
        "--SampleNames '{params.sampleNames}' "
        "--threshold {params.threshold} "
        "--min_cov {params.min_cov} "
        "--minUniqueReads {minNonStatisticalReads} "
        "1> {log.stdout} 2> {log.stderr}"


length_all_samples = len(all_samples)
outputFilterAnalyse = []


filterFileStatistical = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/{prefix}_{unique_id}_outputSpliceLauncher.statistical_junctions{extension}"
filterFileNonStatistical = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/{prefix}_{unique_id}_outputSpliceLauncher.non_statistical_junctions{extension}"
filterFilesSampleStatistical = expand(f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.statistical_junctions{extension}", reads=all_samples)
filterFilesSampleNonStatistical = expand(f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.non_statistical_junctions{extension}", reads=all_samples)
outputFilterAnalyse.append(filterFileStatistical)
outputFilterAnalyse.append(filterFileNonStatistical)
outputFilterAnalyse.append(filterFilesSampleStatistical)
outputFilterAnalyse.append(filterFilesSampleNonStatistical)


filterFilesSampleStatisticalFilter = expand(f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.statistical_junctions.filter{extension}", reads=all_samples)
filterFilesSampleNonStatisticalFilter = expand(f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.non_statistical_junctions.filter{extension}", reads=all_samples)
outputFilterAnalyse.append(filterFilesSampleStatisticalFilter)
outputFilterAnalyse.append(filterFilesSampleNonStatisticalFilter)

maxStatisticalSamples = get_config_value(config["SPLICELAUNCHER"]["POST_ANALYSE"]["STATISTICAL_JUNCTIONS"].get("MAX_SAMPLES"), int(-1))
maxNonStatisticalSamples = get_config_value(config["SPLICELAUNCHER"]["POST_ANALYSE"]["NON_STATISTICAL_JUNCTIONS"].get("MAX_SAMPLES"), int(-1))
thresholdSignificanceLevel = get_config_value(config["SPLICELAUNCHER"]["POST_ANALYSE"]["STATISTICAL_JUNCTIONS"].get("THRESHOLD_SIGNIFICANCE_LEVEL"),int(0))

# This rule filters splice junction results to identify statistical and non_statistical junctions across multiple samples.
# It uses the SpliceLauncher tool to process the aggregated data, applying user-defined thresholds and sample criteria.
# Outputs include filtered lists of statistical and non_statistical junctions for each sample, as well as merged results,
# optionally formatted as text or Excel files and organized into dedicated directories.

rule SpliceLauncher_filter_analyse:
    input:
        outputSpliceLauncher = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/{prefix}_{unique_id}_outputSpliceLauncher{extension}",

    output:
        outputFilterAnalyse

    params:
        Rscript = Rscript,
        length_all_samples = length_all_samples,
        txt = txt,
        Graphics = Graphics,
        sampleNames = sampleNames,
        maxStatisticalSamples = maxStatisticalSamples,
        maxNonStatisticalSamples = maxNonStatisticalSamples,
        minNonStatisticalReads = minNonStatisticalReads,
        thresholdSignificanceLevel = thresholdSignificanceLevel,
        directory = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/",
        script = os.path.join(SCR, "SpliceLauncher_filter_analyse.r") 
    
    log:
        stdout = f"{working_directory}/logs/SpliceLauncher/SpliceLauncher_filter_analyse/{prefix}_{unique_id}.out",
        stderr = f"{working_directory}/logs/SpliceLauncher/SpliceLauncher_filter_analyse/{prefix}_{unique_id}.err"

    shell:
        "{params.Rscript} {params.script} "
        "--input {input.outputSpliceLauncher} "
        "--outputStatistical {output[0]} "
        "--outputNonStatistical {output[1]} "
        "--length {params.length_all_samples} "
        "--maxStatisticalSamples {params.maxStatisticalSamples} "
        "--minNonStatisticalReads {params.minNonStatisticalReads} "
        "--maxNonStatisticalSamples {params.maxNonStatisticalSamples} "
        "--thresholdSignificanceLevel {params.thresholdSignificanceLevel} "
        "--directory {params.directory} "
        "{params.txt}"
        "{params.Graphics}"
        "--SampleNames '{params.sampleNames}' "
        "1> {log.stdout} 2> {log.stderr}"


gtf = get_config_value(config["ANNOTATION"].get("GTF"), "", str)
if not os.path.isfile(gtf):
    gtf = f"{working_directory}/1-raw_data/annotations/{gtf}"
if not os.path.isfile(gtf):
    gtf = -1

color = get_config_value(config["SPLICELAUNCHER"]["SASHIMI_PLOT"].get("COLOR"), "", str)
if not os.path.isfile(color):
    color = f"{working_directory}/1-raw_data/metadata/others/{color}"
if not os.path.isfile(color):
    color = f"{working_directory}/1-raw_data/metadata/others/palette.txt"
if not os.path.isfile(color):
    color = -1

extend_bp = get_config_value(config["SPLICELAUNCHER"]["SASHIMI_PLOT"].get("EXTEND_BP"), int(50))
MinThresholdNbReads = get_config_value(config["SPLICELAUNCHER"]["SASHIMI_PLOT"].get("MIN_THRESHOLD_NB_READS"), -1)
nb_samples = get_config_value(config["SPLICELAUNCHER"]["SASHIMI_PLOT"].get("NUMBER_SAMPLES"), int(4))

# This rule generates sashimi plots for splice junction events using the ggsashimi tool.
# It processes BAM files and filtered event files to create visual plots organized by sample and junction. 
# The rule supports user-defined thresholds, extensions, and additional options like colors and annotations via GTF files.

rule SpliceLauncher_sashimi_plot:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam"),
        event_file = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.{{junction}}.filter{extension}",

    output:
        directory = directory(f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}/sashimi_plot/{{junction}}/")

    params:
        list_bam = expand(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam", zip, reads=all_samples, group=groups),
        python = python,
        ggsashimi = ggsashimi,
        extend_bp = extend_bp,
        MinThresholdNbReads = MinThresholdNbReads,
        nb_samples = nb_samples,
        gtf = gtf,
        color = color,
        script = os.path.join(SCR, "generate_sashimi_plot.py") 


    log:
        stdout = f"{working_directory}/logs/SpliceLauncher/SpliceLauncher_sashimi_plot/{prefix}_{unique_id}/{{group}}/{{reads}}.{{junction}}.out",
        stderr = f"{working_directory}/logs/SpliceLauncher/SpliceLauncher_sashimi_plot/{prefix}_{unique_id}/{{group}}/{{reads}}.{{junction}}.err"
    
    priority:2

    shell:
        "{params.python} {params.script} "
        "-directory {output.directory} "
        "-ggsashimi {params.ggsashimi} "
        "-bam {input.bam} "
        "-event_file {input.event_file} "
        "-extend_bp {params.extend_bp} "
        "-MinThresholdNbReads {params.MinThresholdNbReads} "
        '-list_bam "{params.list_bam}" '
        "-nb_samples {params.nb_samples} "
        "-gtf {params.gtf} "
        "-color {params.color} "
        "1> {log.stdout} 2> {log.stderr}"



## Dynamic construction of inputs/outputs to be moved
move_inputs = []
move_outputs = []

filterFilesSampleStatistical_inputs = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.statistical_junctions{extension}"
filterFilesSampleStatisticalFilter_inputs = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.statistical_junctions.filter{extension}"
filterFilesSampleNonStatistical_inputs = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.non_statistical_junctions{extension}"
filterFilesSampleNonStatisticalFilter_inputs = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.non_statistical_junctions.filter{extension}"
move_inputs.append(filterFilesSampleStatistical_inputs)
move_inputs.append(filterFilesSampleStatisticalFilter_inputs)
move_inputs.append(filterFilesSampleNonStatistical_inputs)
move_inputs.append(filterFilesSampleNonStatisticalFilter_inputs)

filterFilesSampleStatistical = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.statistical_junctions{extension}"
filterFilesSampleStatisticalFilter = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.statistical_junctions.filter{extension}"
filterFilesSampleNonStatistical = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.non_statistical_junctions{extension}"
filterFilesSampleNonStatisticalFilter = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.non_statistical_junctions.filter{extension}"
move_outputs.append(filterFilesSampleStatistical)
move_outputs.append(filterFilesSampleStatisticalFilter)
move_outputs.append(filterFilesSampleNonStatistical)
move_outputs.append(filterFilesSampleNonStatisticalFilter)

if Graphics:
    pdf_inputs = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.pdf"
    pdf_statistical_genes_inputs = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.statistical_genes.pdf"
    pdf_statistical_junctions_inputs = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.statistical_junctions.pdf"
    pdf_non_statistical_junctions_inputs = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.non_statistical_junctions.pdf"
    move_inputs.append(pdf_inputs)
    move_inputs.append(pdf_statistical_genes_inputs)
    move_inputs.append(pdf_statistical_junctions_inputs)
    move_inputs.append(pdf_non_statistical_junctions_inputs)

    pdf = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.pdf"
    pdf_statistical_genes = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.statistical_genes.pdf"
    pdf_statistical_junctions = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.statistical_junctions.pdf"
    pdf_non_statistical_junctions = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.non_statistical_junctions.pdf"
    move_outputs.append(pdf)
    move_outputs.append(pdf_statistical_genes)
    move_outputs.append(pdf_statistical_junctions)
    move_outputs.append(pdf_non_statistical_junctions)

# Moves SpliceLauncher result files (tables, filters, PDFs) from per-read folders
# into their final {group}/{reads} directory structure for organized storage

rule move_SpliceLauncher_outputs:
    input:
        move_inputs

    output:
        move_outputs

    params:
        inputs_directory = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}",
        outputs_directory = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}"

    priority:1

    shell:
        """
        mkdir -p {params.outputs_directory}
        for f in {input}; do
            mv "$f" {params.outputs_directory}/
        done
        """

# This rule generates a sample-level recap table by integrating statistical and non statistical splice junctions.
# It uses an R script to produce annotated outputs including HGVS nomenclature, and sample-specific features.
# The final table summarizes key splicing events in a standardized format.

rule SpliceLauncher_recap:
    input:
        reference = genome,
        mane = mane,
        statisticalFile = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.statistical_junctions.filter{extension}",
        nonStatisticalFile = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.non_statistical_junctions.filter{extension}",

    output:
        RecapFile = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.recap{extension}"

    params:
        Rscript = Rscript,
        sample = f"{{reads}}",
        script = os.path.join(SCR,"recap_file.r")

    log:
        stdout = f"{working_directory}/logs/SpliceLauncher/SpliceLauncher_recap/{prefix}_{unique_id}/{{group}}/{{reads}}.out",
        stderr = f"{working_directory}/logs/SpliceLauncher/SpliceLauncher_recap/{prefix}_{unique_id}/{{group}}/{{reads}}.err"

    shell:
        "{params.Rscript} {params.script} "
        "--statisticalFile {input.statisticalFile} "
        "--nonStatisticalFile {input.nonStatisticalFile} "
        "--reference {input.reference} "
        "--mane {input.mane} "
        "--sample {params.sample} "
        "--output {output.RecapFile} "
        "1> {log.stdout} 2> {log.stderr}"