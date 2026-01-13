process NCBI_DATA {
    tag "${taxon}"
    label 'process_medium'

    input:
    tuple val(taxon), val(min_length)

    output:
    tuple val(taxon), path("*.fa.gz"), path("*.json"), emit: data
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when
    

    script:
    args   = task.ext.args ?: ''
    tool = "epitome_ncbi.py"

    edirect_cmd = """
    #---- NCBI EDirect -----#
    # Gather subtype information (Not included in Datasets report)

    cat data/datasets-genome.fa \\
        | grep '>' \\
        | cut -f 1 -d ' ' \\
        | tr -d '>' \\
        > data/accessions.txt
        
    epost -db nuccore -input data/accessions.txt \\
        | efetch -format docsum -mode json \\
        > data/edirect.json
    """

    """
    mkdir data/ || true

    #---- NCBI Datasets Genome -----#
    # Downloads genome assemblies and associated sample data.
    # Sample data is currently lacking subtype information (Downloaded below)

    datasets download virus genome taxon \\
        "${taxon}" \\
        --no-progressbar \\
        --include genome \\
        ${params.complete_only ? '--complete-only' : '' }
    
    unzip ncbi_dataset.zip
    mv ncbi_dataset/data/genomic.fna data/datasets-genome.fa
    mv ncbi_dataset/data/data_report.jsonl data/datasets-genome.jsonl
    
    #---- NCBI Datasets Taxon -----#
    # Used to identify the species-level classification

    datasets summary taxonomy taxon \\
        "${taxon}" \\
        --rank species \\
        > data/datasets-taxonomy.json
    
    ${ params.edirect ?  edirect_cmd : '' }
    
    #---- COMBINE ----#
    ${tool} \\
        --taxon ${taxon.replace(' ', '_')} \\
        --min-length ${min_length} \\
        --datasets_genome_fasta data/datasets-genome.fa \\
        --datasets_genome_json data/datasets-genome.jsonl \\
        --datasets_taxonomy data/datasets-taxonomy.json \\
        ${params.edirect ? '--edirect data/edirect.json' : ''}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi-datasets: \$(datasets --version | sed 's/.*: //g')
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
