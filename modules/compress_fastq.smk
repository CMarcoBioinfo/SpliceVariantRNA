# This rule compresses FASTQ files.
# A marker file is created to indicate the task's completion, whether or not compression is performed.
# If the file is already compressed (.gz), no further action is taken.
# Otherwise, the file is compressed using pigz (multi-threaded compression).

rule compress_fastq:
    input:
        read = lambda wildcards: get_inputs(wildcards,0, length_fastp),
        preprocess = f"{working_directory}/.tmp/samples/{{id}}.pre"
        
    output:
        process = f"{working_directory}/.tmp/samples/{{id}}.pro"

    params:
        id = lambda wildcards: wildcards.id,
        pigz = pigz

    threads: 
        config["COMPRESS"]["THREADS"]
    
    priority: 5

    log:
        stdout = f"{working_directory}/logs/compress_fastq/{{id}}.out",
        stderr = f"{working_directory}/logs/compress_fastq/{{id}}.err"

    run:
        if (str(input.read).endswith('.gz') ):
            if not os.path.exists(output.process):
                shell("touch {output.process}")
        else : 
            shell("{params.pigz} "
            "--best "
            "--verbose "
            "--processes {threads} {input.read} "
            "1> {log.stdout} 2> {log.stderr} && touch {output.process}")


# This rule processes paired FASTQ files by creating symbolic links to their final destination.
# It ensures the output directory exists and verifies that input files are compressed (.gz).
# If an input file is not compressed, the rule appends the '.gz' extension.
# Output files are linked dynamically to streamline the workflow and maintain organization.

rule processed_fastq:
    input:
        read1 = lambda wildcards: get_inputs(wildcards, 1, length_fastp),
        read2 =  lambda wildcards: get_inputs(wildcards, 2, length_fastp),
        process1 = f"{working_directory}/.tmp/samples/reads1/{{id}}.pro",
        process2 = f"{working_directory}/.tmp/samples/reads2/{{id}}.pro"

    output:
        read1 = os.path.abspath(f"{working_directory}/1-raw_data/fastq/{{id}}.1.fastq.gz"),
        read2 = os.path.abspath(f"{working_directory}/1-raw_data/fastq/{{id}}.2.fastq.gz")
        
    params:
        id = lambda wildcards: wildcards.id,
        directory = f"{working_directory}/1-raw_data/fastq/"

    priority: 5

    run:
        create_directory_if_not_exists(params["directory"])
        input.read1 = f"{input.read1}.gz" if not is_gz(input.read1) else input.read1
        input.read2 = f"{input.read2}.gz" if not is_gz(input.read2) else input.read2
        shell("ln -sfn "
        "{input.read1} {output.read1} && "
        "ln -sfn "
        "{input.read2} {output.read2}")