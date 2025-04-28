# This rule computes detailed BAM statistics using Samtools.
# It provides comprehensive alignment metrics, including insert size distribution, coverage depth, and quality scores.

rule samtools_stats:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam"),
        bai = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam.bai")

    output:
        stats = f"{path_qc}/BAM/{name_genome}/samtools/{{group}}/{{reads}}.stats"

    params:
        samtools = samtools

    threads:
        config["QUALITY_CONTROL"]["THREADS"]

    log:
        stderr = f"{working_directory}/logs/samtools/{{group}}/{{reads}}_stats.err"

    shell:
        "{params.samtools} stats "
        "-@ {threads} {input.bam} > {output.stats} 2> {log.stderr}"


# This rule runs samtools flagstat to summarize alignment quality.
# It provides key statistics such as mapped read percentages, duplicate rate, and unmapped reads.

rule samtools_flagstat:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam"),
        bai = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam.bai")

    output:
        flagstat = f"{path_qc}/BAM/{name_genome}/samtools/{{group}}/{{reads}}.flagstat"

    params:
        samtools = samtools

    threads:
        config["QUALITY_CONTROL"]["THREADS"]

    log:
        stderr = f"{working_directory}/logs/samtools/{{group}}/{{reads}}_flagstat.err"

    shell:
        "{params.samtools} flagstat "
        "-@ {threads} {input.bam} > {output.flagstat} 2> {log.stderr}"


add_inputs = []
if use_RSeQC:
    bam_stat = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.bam_stat.txt", zip, reads=all_samples, group=groups)
    pos_DupRate = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.pos.DupRate.xls", zip, reads=all_samples, group=groups)
    seq_DupRate = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.seq.DupRate.xls", zip, reads=all_samples, group=groups)
    DupRate_plot_pdf = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.DupRate_plot.pdf", zip, reads=all_samples, group=groups)
    DupRate_plot_r = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.DupRate_plot.r", zip, reads=all_samples, group=groups)
    add_inputs.append(bam_stat)
    add_inputs.append(pos_DupRate)
    add_inputs.append(seq_DupRate)
    add_inputs.append(DupRate_plot_pdf)
    add_inputs.append(DupRate_plot_r)

    if use_bed12:
        infer_experiment = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.infer_experiment.txt", zip, reads=all_samples, group=groups)
        read_distribution = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.read_distribution.txt", zip, reads=all_samples, group=groups)
        geneBody_coverage_txt = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.geneBodyCoverage.txt", zip, reads=all_samples, group=groups)
        geneBody_coverage_r = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.geneBodyCoverage.r", zip, reads=all_samples, group=groups)
        geneBody_coverage_pdf = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.geneBodyCoverage.curves.pdf", zip, reads=all_samples, group=groups)
        junction_annotation_xls = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.junction.xls", zip, reads=all_samples, group=groups)
        junction_annotation_bed = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.junction.bed", zip, reads=all_samples, group=groups)
        junction_annotation_inter = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.junction.Interact.bed", zip, reads=all_samples, group=groups)
        junction_annotation_r = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.junction_plot.r", zip, reads=all_samples, group=groups)
        junction_saturation_pdf = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.junctionSaturation_plot.pdf", zip, reads=all_samples, group=groups)
        junction_saturation_r = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.junctionSaturation_plot.r", zip, reads=all_samples, group=groups)
        inner_distance_txt = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.inner_distance.txt", zip, reads=all_samples, group=groups)
        inner_distance_freq = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.inner_distance_freq.txt", zip, reads=all_samples, group=groups)
        inner_distance_pdf = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.inner_distance_plot.pdf", zip, reads=all_samples, group=groups)
        inner_distance_r = expand(f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.inner_distance_plot.r", zip, reads=all_samples, group=groups)
        add_inputs.append(infer_experiment)
        add_inputs.append(read_distribution)
        add_inputs.append(geneBody_coverage_txt)
        add_inputs.append(geneBody_coverage_r)
        add_inputs.append(geneBody_coverage_pdf)
        add_inputs.append(junction_annotation_xls)
        add_inputs.append(junction_annotation_bed)
        add_inputs.append(junction_annotation_inter)
        add_inputs.append(junction_annotation_r)
        add_inputs.append(junction_saturation_pdf)
        add_inputs.append(junction_saturation_r)
        add_inputs.append(inner_distance_txt)
        add_inputs.append(inner_distance_freq)
        add_inputs.append(inner_distance_pdf)
        add_inputs.append(inner_distance_r)

# This rule consolidates all BAM QC reports into a MultiQC summary.
# It combines alignment statistics, flagstat results, and logs from the mapping process.
# The final output is a single interactive HTML report.

rule multiqc_bam:
    input:
        markdup_metric = expand(os.path.abspath(f"{path_bam}{name_genome}/mapping/log_MarkDuplicates/{{group}}/{{reads}}.markdup.metrics.txt"), zip, reads=all_samples, group=groups),
        samtools_stats = expand(f"{path_qc}/BAM/{name_genome}/samtools/{{group}}/{{reads}}.stats", zip, reads=all_samples, group=groups),
        samtools_flagstats = expand(f"{path_qc}/BAM/{name_genome}/samtools/{{group}}/{{reads}}.flagstat", zip, reads=all_samples, group=groups),
        log_final = expand(f"{path_bam}{name_genome}/mapping/log_star/{{group}}/{{reads}}_Log.final.out", zip, reads=all_samples, group=groups),
        log = expand(f"{path_bam}{name_genome}/mapping/log_star/{{group}}/{{reads}}_Log.out", zip, reads=all_samples, group=groups),
        log_progress = expand(f"{path_bam}{name_genome}/mapping/log_star/{{group}}/{{reads}}_Log.progress.out", zip, reads=all_samples, group=groups),
        tab = expand(f"{path_bam}{name_genome}/mapping/log_star/{{group}}/{{reads}}_SJ.out.tab", zip, reads=all_samples, group=groups),
        extra = add_inputs

    output:
        directory_data = directory(f"{path_qc}/multiqc/BAM/{name_genome}/{prefix}_{unique_id}_data/"),
        html = f"{path_qc}/multiqc/BAM/{name_genome}/{prefix}_{unique_id}.html"

    params:
        name = f"{prefix}_{unique_id}",
        multiqc = multiqc,
        path = f"{path_qc}/multiqc/BAM/{name_genome}/"

    threads:
        config["QUALITY_CONTROL"]["THREADS"]

    log:
        stdout = f"{working_directory}/logs/multiqc/BAM/{prefix}_{unique_id}.out",
        stderr = f"{working_directory}/logs/multiqc/BAM/{prefix}_{unique_id}.err"

    shell:
        "{params.multiqc} "
        "{input.markdup_metric} {input.samtools_stats} {input.samtools_flagstats} {input.log_final} {input.log} {input.log_progress} {input.tab} {input.extra} "
        "--filename {params.name} "
        '--title "Quality Control of bam files" '
        "--dirs --dirs-depth 1 "
        "-o {params.path} 1> {log.stdout} 2> {log.stderr}"