process PICARD_MARKILLUMINAADAPTERS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.2.0--hdfd78af_0' :
        'biocontainers/picard:3.2.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.txt"), emit: stats
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!task.memory) {
        log.warn '[Picard MarkIlluminaAdapters] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    }
    def avail_mem = task.memory ? (task.memory.mega*0.8).intValue() : 3072
    def input = meta.single_end ? "--INPUT ${seqs}" : "--INPUT ${seqs}"
    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkIlluminaAdapters \\
        ${args} \\
        ${input} \\
        --OUTPUT ${prefix}_markilluminaadapters.bam \\
        --METRICS ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkIlluminaAdapters --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
