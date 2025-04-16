# This rule performs quality trimming and adapter detection on paired-end FASTQ files
# using the tool 'fastp'. It generates trimmed FASTQ files as output, along with
# corresponding HTML and JSON reports for quality control analysis.
# Directories for trimmed outputs and quality control reports are created if they do not exist.

rule fastp_trimming:
    input:
        read1 = os.path.abspath(f"{working_directory}/1-raw_data/fastq/{{reads}}.1.fastq.gz"),
        read2 = os.path.abspath(f"{working_directory}/1-raw_data/fastq/{{reads}}.2.fastq.gz")

    output:
        trimmed_read1 = os.path.abspath(f"{path_fastq}/{{reads}}.{{length_fastp}}bp.1.fastq.gz"),
        trimmed_read2 = os.path.abspath(f"{path_fastq}/{{reads}}.{{length_fastp}}bp.2.fastq.gz"),
        html = f"{path_qc}/fastp_trimming/{{reads}}.{{length_fastp}}bp.html",
        json = f"{path_qc}/fastp_trimming/{{reads}}.{{length_fastp}}bp.json"

    params:
        fastp = fastp,
        length = length_fastp,
        directory_trimmed = os.path.abspath(f"{path_fastq}"),
        directory_fastp = os.path.abspath(f"{path_qc}/fastp_trimming/")

    threads:
        config["TRIMMING"]["THREADS"]

    log:
        stdout = f"{working_directory}/logs/fastp_trimming/{{reads}}.{{length_fastp}}bp.out",
        stderr = f"{working_directory}/logs/fastp_trimming/{{reads}}.{{length_fastp}}bp.err"

    run:
        create_directory_if_not_exists(params.directory_trimmed)
        create_directory_if_not_exists(params.directory_fastp)
        shell("{params.fastp} "
        "--thread {threads} "
        "--length_required {params.length} "
        "--detect_adapter_for_pe "
        "--in1 {input.read1} "
        "--in2 {input.read2} "
        "--out1 {output.trimmed_read1} "
        "--out2 {output.trimmed_read2} "
        "--html {output.html} "
        "--json {output.json} "
        "1> {log.stdout} 2> {log.stderr}")