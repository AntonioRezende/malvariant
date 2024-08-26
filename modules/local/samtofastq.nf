process PICARD_SAMTOFASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.2.0--hdfd78af_0' :
        'biocontainers/picard:3.2.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("*.fastq"), emit: reads
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!task.memory) {
        log.warn '[Picard SamToFastq] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    }
    def avail_mem = task.memory ? (task.memory.mega*0.8).intValue() : 3072
    def input = "--INPUT ${seqs}"
    """
    picard \\
        -Xmx${avail_mem}M \\
        SamToFastq \\
        ${args} \\
        ${input} \\
        --FASTQ ${prefix}_interleaved.fastq \\
        --CLIPPING_ATTRIBUTE XT \\
        --CLIPPING_ACTION 2 \\
        --INTERLEAVE true \\
        -NON_PF true


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SamToFastq --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
