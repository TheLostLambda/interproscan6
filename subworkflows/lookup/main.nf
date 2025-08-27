include { PREPARE_LOOKUP; LOOKUP_MATCHES } from "../../modules/lookup"

workflow LOOKUP {
    // Prepare connection and retrieve precalculated matched from the InterPro API
    take:
    ch_seqs               // fasta files of protein sequences to analyse
    apps                  // member db analyses to run
    db_releases           // map: [db: version, dirpath]           
    interproscan_version  // major.minor interproscan version number
    workflow_manifest     // map, from nextflow.conf
    matches_api_url       // str, from cmd-line
    chunk_size            // int
    max_retries           // int

    main:

    PREPARE_LOOKUP(
        matches_api_url,
        db_releases,
        interproscan_version,
        workflow_manifest
    )
    api_info = PREPARE_LOOKUP.out[0]

    // Pass the url from PREPARE_LOOKUP to LOOKUP_MATCHES to avoid an 'invalid method invocation' error
    api_info
        .combine(ch_seqs)
        .map { url, apiApps, err, index, fasta ->
            tuple(index, fasta, apps, apiApps, url, chunk_size, max_retries, err)
        }
        .set { lookup_input }

    LOOKUP_MATCHES(lookup_input)
    precalculatedMatches = LOOKUP_MATCHES.out[0]
    noMatchesFasta = LOOKUP_MATCHES.out[1]
    noApiFasta = LOOKUP_MATCHES.out[2]

    emit:
    precalculatedMatches
    noMatchesFasta
    noApiFasta
}


// Build and emit the sets of applications
