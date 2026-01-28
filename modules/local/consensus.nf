process CONSENSUS {
    tag "${prefix}"
    label 'process_medium'
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), val(cluster), path(aln)

    output:
    tuple val(taxon), val(segment), path("${prefix}.fa.gz"), emit: fa
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}-${cluster}"
    tool = "epitome_consensus.py"
    """
    # run script
    ${tool} \\
        --prefix ${prefix} \\
        --aln ${aln}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
