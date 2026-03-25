# Construction of sample in group
groups_raw = groups[:]      # liste brute, alignée avec all_samples
groups_uniq = sorted(set(groups))

group_to_samples = {}

for g, s in zip(groups_raw, all_samples):
    group_to_samples.setdefault(g, []).append(s)


# Create a recap ZIP for each group by collecting all sample recap files.
# → Produces an intermediate group-level ZIP archive.

rule zip_group_recap:
    input:
        recap = lambda wc: expand(f"{path_results}/{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.recap{extension}", group=wc.group, reads=group_to_samples[wc.group])
    
    output:
        group_zip = temp(f"{path_results}/{unique_id}_results/zip/{{group}}/{{group}}_recap.zip")
    
    shell:
        "zip -j {output.group_zip} {input.recap}"


# Combine all group recap ZIPs into a single final recap ZIP for the run.
# → Produces the final recap deliverable.

rule zip_recap:
    input:
        zip_group = expand(f"{path_results}/{unique_id}_results/zip/{{group}}/{{group}}_recap.zip", group=groups_uniq)
    
    output:
        zip_recap = f"{path_results}/{unique_id}_results/zip/{unique_id}_recap.zip"
    
    shell:
        "zip -j {output.zip_recap} {input.zip_group}"


# Package all Sashimi plot outputs for a single sample.
# → Produces an intermediate per-sample Sashimi ZIP.
# Note: parent directory must exist before running zip.

rule zip_sample_sashimi:
    input:
        nonstat = f"{path_results}/{unique_id}_results/samples_results/{{group}}/{{reads}}/sashimi_plot/non_statistical_junctions/",
        stat = f"{path_results}/{unique_id}_results/samples_results/{{group}}/{{reads}}/sashimi_plot/statistical_junctions/"
    
    output:
        zip_sample = temp(os.path.abspath(f"{path_results}/{unique_id}_results/zip/{{group}}/{{reads}}_sashimi.zip"))

    params:
        directory = os.path.abspath(f"{path_results}/{unique_id}_results/samples_results/{{group}}/{{reads}}")
    
    shell:
        "cd {params.directory} && "
        "zip -r {output.zip_sample} sashimi_plot"


# Combine all per-sample Sashimi ZIPs for a group.
# → Produces an intermediate group-level Sashimi ZIP.

rule zip_group_sashimi:
    input:
        zip_samples = lambda wc: expand(os.path.abspath(f"{path_results}/{unique_id}_results/zip/{{group}}/{{reads}}_sashimi.zip"), group=wc.group, reads=group_to_samples[wc.group])

    output:
        zip_group = temp(f"{path_results}/{unique_id}_results/zip/{{group}}/{{group}}_sashimi.zip")

    shell:
        "zip -j {output.zip_group} {input.zip_samples}"


# Combine all group-level Sashimi ZIPs into a single final archive.
# → Produces the final Sashimi deliverable for the run.

rule zip_sashimi:
    input:
        zip_groups = expand(f"{path_results}/{unique_id}_results/zip/{{group}}/{{group}}_sashimi.zip", group=groups_uniq)
    
    output:
        zip_sashimi = f"{path_results}/{unique_id}_results/zip/{unique_id}_sashimi.zip"
    
    shell:
        "zip -j {output.zip_sashimi} {input.zip_groups}"


# Package all SpliceLauncher outputs for a single sample:
#   - statistical / non-statistical junctions
#   - filtered junctions
#   - optional PDF graphics (if Graphics=True)
# → Produces an intermediate per-sample SpliceLauncher ZIP.

rule zip_sample_SpliceLauncher:
    input:
        f"{path_results}/{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.statistical_junctions{extension}",
        f"{path_results}/{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.non_statistical_junctions{extension}",
        f"{path_results}/{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.statistical_junctions.filter{extension}",
        f"{path_results}/{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.non_statistical_junctions.filter{extension}",
        *(
            [
                f"{path_results}/{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.pdf",
                f"{path_results}/{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.statistical_genes.pdf",
                f"{path_results}/{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.statistical_junctions.pdf",
                f"{path_results}/{unique_id}_results/samples_results/{{group}}/{{reads}}/{{reads}}.non_statistical_junctions.pdf",
            ]
            if Graphics else []
            )

    output:
        zip_sample = temp(f"{path_results}/{unique_id}_results/zip/{{group}}/{{reads}}_SpliceLauncher.zip")

    shell:
        "zip -j {output.zip_sample} {input}"


# Combine all per-sample SpliceLauncher ZIPs for a group.
# → Produces an intermediate group-level SpliceLauncher ZIP.

rule zip_group_SpliceLauncher:
    input:
        zip_samples = lambda wc: expand(f"{path_results}/{unique_id}_results/zip/{{group}}/{{reads}}_SpliceLauncher.zip", group=wc.group, reads=group_to_samples[wc.group])

    output:
        zip_group = temp(f"{path_results}/{unique_id}_results/zip/{{group}}/{{group}}_SpliceLauncher.zip")

    shell:
        "zip -j {output.zip_group} {input.zip_samples}"


# Combine all group-level SpliceLauncher ZIPs and global files:
#   - count report
#   - main XLSX output
#   - statistical / non-statistical junction summaries
# → Produces the final SpliceLauncher deliverable for the run.

rule zip_SpliceLauncher:
    input:
        count_report = f"{path_results}/{unique_id}_report_{date}.txt",
        outputSpliceLauncher = f"{path_results}/{unique_id}_results/{unique_id}_outputSpliceLauncher{extension}",
        filterFileStatistical = f"{path_results}/{unique_id}_results/{unique_id}_outputSpliceLauncher.statistical_junctions{extension}",
        filterFileNonStatistical = f"{path_results}/{unique_id}_results/{unique_id}_outputSpliceLauncher.non_statistical_junctions{extension}",
        zip_groups = expand(f"{path_results}/{unique_id}_results/zip/{{group}}/{{group}}_SpliceLauncher.zip", group=groups_uniq)
    
    output:
        zip_SpliceLauncher = f"{path_results}/{unique_id}_results/zip/{unique_id}_SpliceLauncher.zip"
    
    shell:
        "zip -j {output.zip_SpliceLauncher} {input.count_report} {input.outputSpliceLauncher} {input.filterFileStatistical} {input.filterFileNonStatistical} {input.zip_groups}"


# Package all MultiQC outputs (FASTQ raw, trimming, BAM QC).
# → Produces the final QC ZIP deliverable for the run.

rule zip_qc:
    input:
        # FASTQ RAW
        raw_data = f"{path_qc}/multiqc/fastq_raw/{prefix}_{unique_id}_data/",
        raw_html = f"{path_qc}/multiqc/fastq_raw/{prefix}_{unique_id}.html",

        # FASTQ TRIMMED (optionnel)
        *(
            [
                f"{path_qc}/multiqc/fastq_trimmed/{prefix}_{unique_id}_data/",
                f"{path_qc}/multiqc/fastq_trimmed/{prefix}_{unique_id}.html",
            ]
            if use_trimming else []
        ),

        # BAM QC (optionnel)
        *(
            [
                f"{path_qc}/multiqc/BAM/{name_genome}/{prefix}_{unique_id}_data/",
                f"{path_qc}/multiqc/BAM/{name_genome}/{prefix}_{unique_id}.html",
            ]
            if use_mapping else []
        )

    output:
        zip_qc = os.path.abspath(f"{path_results}/{unique_id}_results/zip/{unique_id}_qc.zip")
    
    run:
        cmd = f"cd {path_qc}/multiqc && "
        cmd += f"zip -r {output.zip_qc} "
        cmd += f"fastq_raw/{prefix}_{unique_id}.html "
        cmd += f"fastq_raw/{prefix}_{unique_id}_data "

        if use_trimming:
            cmd += f"fastq_trimmed/{prefix}_{unique_id}.html "
            cmd += f"fastq_trimmed/{prefix}_{unique_id}_data "

        if use_mapping:
            cmd += f"BAM/{name_genome}/{prefix}_{unique_id}.html "
            cmd += f"BAM/{name_genome}/{prefix}_{unique_id}_data "

        shell(cmd)