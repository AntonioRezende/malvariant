process PICARD_UMIDUPLICATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.2.0--hdfd78af_0' :
        'biocontainers/picard:3.2.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*_duplicate_metrics.txt"), emit: dmetrics
    tuple val(meta), path("*_umi_metrics.txt"), emit: umetrics
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
    def input = "--INPUT ${bam}"
    """
    picard \\
        -Xmx${avail_mem}M \\
        UmiAwareMarkDuplicatesWithMateCigar \\
        ${args} \\
        ${input} \\
        --OUTPUT ${prefix}_UMIMARKED.bam \\
        --UMI_METRICS_FILE ${prefix}_umi_metrics.txt \\
        --METRICS_FILE ${prefix}_duplicate_metrics.txt

        


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SamToFastq --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
//--UMI_METRICS_FILE ${prefix}_umi_metrics.txt
//--METRICS_FILE ${prefix}_duplicate_metrics.txt