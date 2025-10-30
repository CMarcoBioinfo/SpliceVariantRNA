rule RSeQC_bam_stat:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam"),
        bai = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam.bai")

    output:
        bam_stat = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.bam_stat.txt"

    params:
        RSeQC = RSeQC
    
    log:
        stderr = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}_bam_stat.err"

    shell:
        "{params.RSeQC}/bam_stat.py "
        "-i {input.bam} "
        "> {output.bam_stat} 2> {log.stderr}"


rule RSeQC_read_duplication:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam"),
        bai = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam.bai")

    output:
        pos_DupRate = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.pos.DupRate.xls",
        seq_DupRate = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.seq.DupRate.xls",
        DupRate_plot_pdf = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.DupRate_plot.pdf",
        DupRate_plot_r = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.DupRate_plot.r"

    params:
        RSeQC = RSeQC,
        prefix = lambda wildcards: f"{path_qc}/BAM/{name_genome}/RSeQC/{wildcards.group}/{wildcards.reads}"
    
    log:
        stdout = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}_read_duplication.out",
        stderr = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}_read_duplication.err"

    shell:
        "{params.RSeQC}/read_duplication.py "
        "-i {input.bam} "
        "-o {params.prefix} "
        "1> {log.stdout} 2> {log.stderr}"


rule RSeQC_infer_experiment:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam"),
        bai = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam.bai"),
        bed = bed12
    
    output:
        infer_experiment = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.infer_experiment.txt"

    params:
        RSeQC = RSeQC
    
    log:
        stderr = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}_infer_experiment.err"

    shell:
        "{params.RSeQC}/infer_experiment.py "
        "-i {input.bam} "
        "-r {input.bed} " 
        "> {output.infer_experiment} 2> {log.stderr}"


rule RSeQC_read_distribution:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam"),
        bai = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam.bai"),
        bed = bed12

    output:
        read_distribution = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.read_distribution.txt"

    params:
        RSeQC = RSeQC
    
    log:
        stderr = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}_read_distribution.err"

    shell:
        "{params.RSeQC}/read_distribution.py "
        "-i {input.bam} "
        "-r {input.bed} " 
        "> {output.read_distribution} 2> {log.stderr}"
        

rule RSeQC_geneBody_coverage:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam"),
        bai = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam.bai"),
        bed = bed12

    output:
        geneBody_coverage_txt = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.geneBodyCoverage.txt",
        geneBody_coverage_r = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.geneBodyCoverage.r",
        geneBody_coverage_pdf = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.geneBodyCoverage.curves.pdf"

    params:
        RSeQC = RSeQC,
        prefix = lambda wildcards: f"{path_qc}/BAM/{name_genome}/RSeQC/{wildcards.group}/{wildcards.reads}"
    
    log:
        stdout = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}_geneBody_coverage.out",
        stderr = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}_geneBody_coverage.err"

    shell:
        "{params.RSeQC}/geneBody_coverage.py "
        "-i {input.bam} "
        "-r {input.bed} "
        "-o {params.prefix} "
        "1> {log.stdout} 2> {log.stderr}"


rule RSeQC_junction_annotation:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam"),
        bai = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam.bai"),
        bed = bed12

    output:
        junction_annotation_xls = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.junction.xls",
        junction_annotation_bed = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.junction.bed",
        junction_annotation_inter = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.junction.Interact.bed",
        junction_annotation_r = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.junction_plot.r"

    params:
        RSeQC = RSeQC,
        prefix = lambda wildcards: f"{path_qc}/BAM/{name_genome}/RSeQC/{wildcards.group}/{wildcards.reads}"
    
    log:
        stdout = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}.junction_annotation.out",
        stderr = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}.junction_annotation.err"

    shell:
        "{params.RSeQC}/junction_annotation.py "
        "-i {input.bam} "
        "-r {input.bed} "
        "-o {params.prefix} "
        "1> {log.stdout} 2> {log.stderr}"


rule RSeQC_junction_saturation:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam"),
        bai = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam.bai"),
        bed = bed12

    output:
        junction_saturation_pdf = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.junctionSaturation_plot.pdf",
        junction_saturation_r = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.junctionSaturation_plot.r"

    params:
        RSeQC = RSeQC,
        prefix = lambda wildcards: f"{path_qc}/BAM/{name_genome}/RSeQC/{wildcards.group}/{wildcards.reads}"
    
    log:
        stdout = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}.junction_saturation.out",
        stderr = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}.junction_saturation.err"

    shell:
        "{params.RSeQC}/junction_saturation.py "
        "-i {input.bam} "
        "-r {input.bed} "
        "-o {params.prefix} "
        "1> {log.stdout} 2> {log.stderr}"


rule RSeQC_inner_distance:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam"),
        bai = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{group}}/{{reads}}.markdup.bam.bai"),
        bed = bed12

    output:
        inner_distance_txt = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.inner_distance.txt",
        inner_distance_freq = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.inner_distance_freq.txt",
        inner_distance_pdf = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.inner_distance_plot.pdf",
        inner_distance_r = f"{path_qc}/BAM/{name_genome}/RSeQC/{{group}}/{{reads}}.inner_distance_plot.r"

    params:
        RSeQC = RSeQC,
        prefix = lambda wildcards: f"{path_qc}/BAM/{name_genome}/RSeQC/{wildcards.group}/{wildcards.reads}"
    
    log:
        stdout = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}_inner_distance.out",
        stderr = f"{working_directory}/logs/RSeQC/{{group}}/{{reads}}_inner_distance.err"

    shell:
        "{params.RSeQC}/inner_distance.py "
        "-i {input.bam} "
        "-r {input.bed} "
        "-o {params.prefix} "
        "1> {log.stdout} 2> {log.stderr}"