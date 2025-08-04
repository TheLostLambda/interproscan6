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
    // Initialise channels for outputs
    precalculatedMatches = Channel.empty()
    noMatchesFasta       = Channel.empty()
    noLookupFasta        = Channel.empty()

    PREPARE_LOOKUP(
        matches_api_url,
        db_releases,
        interproscan_version,
        workflow_manifest,
        apps
    )
    matchesApiApps = PREPARE_LOOKUP.out[0]

    matchesApiApps
        .filter { it }
        .combine(ch_seqs)
        .map { api_apps, index, fasta ->
            tuple(index, fasta, matches_api_url, apps, api_apps, chunk_size, max_retries)
        }
        .set { lookup_input }

    LOOKUP_MATCHES(lookup_input)

    precalculatedMatches = LOOKUP_MATCHES.out[0]
    noMatchesFasta       = LOOKUP_MATCHES.out[1]
    noLookupFasta        = LOOKUP_MATCHES.out[2]

    // fallback when no API URL is available
    matchesApiApps
        .filter { !it }
        .map { _ ->
            precalculatedMatches = Channel.empty()
            noMatchesFasta       = ch_seqs
            noLookupFasta        = Channel.empty()
        }

    emit:
    precalculatedMatches
    noMatchesFasta
    noLookupFasta
}