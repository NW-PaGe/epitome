process CONDENSE {
    tag "${prefix}"
    label "process_medium"
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), path(seqs), path(clusters)

    output:
    tuple val(taxon), val(segment), path("${prefix}.condensed.jsonl.gz"), emit: results
    path "versions.yml",                                                  emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    tool = "epitome_condense.py"
    """
    ${tool} \\
        --taxon "${taxon}" \\
        --segment "${segment}" \\
        --clusters ${clusters} \\
        --dist ${params.dist_threshold} \\
        --window_size ${params.window_size} \\
        --ksize ${params.ksize} \\
        --scaled ${params.scaled} \\
        --fasta ${seqs}
        
    # version info 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}