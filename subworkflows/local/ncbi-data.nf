//
// Check input samplesheet and get read channels
//

include { NCBI_DATA    } from '../../modules/local/ncbi-data'

workflow NCBI_DATA_SUBWF {
    take:
    ch_input // [ va(taxon), val(segment), path(assembly), ... ]

    main:
    ch_versions = Channel.empty()

    /*
    =============================================================================================================================
        GATHER DATA FROM NCBI
    =============================================================================================================================
    */
    // MODULE: Gather NCBI data for taxon
    NCBI_DATA(
        ch_input.map{ [ it.taxon, it.min_length ] }.unique()
    )
    ch_versions = ch_versions.mix(NCBI_DATA.out.versions.first())

    NCBI_DATA
        .out
        .data
        .map{ taxon, assembly, metadata -> [ taxon: taxon, assembly: assembly, metadata: metadata ] }
        .set{ch_input}

    emit:
    input    = ch_input
    versions = ch_versions
}