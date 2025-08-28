include { LOOKUP_MATCHES } from "../../modules/lookup"

workflow LOOKUP {
    // Prepare connection and retrieve precalculated matched from the InterPro API
    take:
    ch_seqs               // fasta files of protein sequences to analyse
    matches_api_apps      // member db analyses to run that are in the matches API
    db_releases           // map: [db: version, dirpath]           
    interproscan_version  // major.minor interproscan version number
    url                   // str, url to matches api
    chunk_size            // int
    max_retries           // int

    main:
    
    ch_seqs
        .map { index, fasta ->
            tuple(index, fasta, matches_api_apps, url, chunk_size, max_retries)
        }
        .set { lookup_input }

    LOOKUP_MATCHES(lookup_input)
    precalculatedMatches = LOOKUP_MATCHES.out[0]
    noMatchesFasta       = LOOKUP_MATCHES.out[1]
    noLookupFasta        = LOOKUP_MATCHES.out[2]

    emit:
    precalculatedMatches
    noMatchesFasta
}
