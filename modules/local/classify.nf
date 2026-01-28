process CLASSIFY {
    tag "${prefix}"
    label 'process_high'
    stageInMode 'copy'

    input:
    tuple val(taxon), path(jsonl)

    output:
    tuple val(taxon), path("*.jsonl.gz"), emit: jsonl
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    prefix = "${taxon.replaceAll(' ','_')}"
    tool = "epitome_classify.py"
    """
    ${tool} \\
        ${args} \\
        ${jsonl}
        

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
