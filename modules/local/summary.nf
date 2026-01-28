process SUMMARY {
    label 'process_low'
    tag "${prefix}"

    input:
    tuple val(taxon), val(segment), path(qc_json), path(clusters_json), path(meta_csv), val(method)

    output:
    tuple val(taxon), val(segment), path("${prefix}.jsonl.gz"),   emit: json
    tuple val(taxon), val(segment), path("*.jpg"),                emit: plot, optional: true
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}.${method}"
    tool = "epitome_summary.py"
    """
    # run script
    ${tool} \\
        --taxon "${taxon}" \\
        --segment "${segment}" \\
        --method ${method} \\
        --qc ${qc_json} \\
        --clusters ${clusters_json} \\
        --meta ${meta_csv} \\
        --z_threshold ${params.z_threshold}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
