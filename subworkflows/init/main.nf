workflow INIT_PIPELINE {
    // Validate pipeline input parameters
    take:
    input
    applications
    apps_config
    datadir
    formats
    outdir
    outprefix
    no_matches_api
    interpro_version
    skip_intepro
    skip_applications
    goterms
    pathways

    main:
    // Check the input
    fasta = InterProScan.resolveFile(input)
    if (!fasta) {
        log.error "No such file: ${input}"
        exit 1
    }

    // Applications validation
    (apps, error) = InterProScan.validateApplications(applications, skip_applications, apps_config)
    if (!apps) {
        log.error error
        exit 1
    }

    if (skip_intepro && (goterms || pathways)) {
        log.error "--skip_intepro is mutually exclusive with --goterms and --pathways"
        exit 1
    }

    // Check valid output file formats were provided
    (formats, error) = InterProScan.validateFormats(formats)
    if (error) {
        log.error error
        exit 1
    }

    apps_with_data = InterProScan.getAppsWithData(apps, apps_config)
    if (apps_with_data.size() > 0) {
        if (datadir == null) {
            log.error "'--datadir <DATA-DIR>' is required for the selected applications."
            exit 1
        }
    
        (datadir, error) = InterProScan.resolveDirectory(datadir, false, false)
        if (datadir == null) {
            log.error error
            exit 1
        }
    } else {
        datadir = null
    }
  
    version = InterProScan.validateInterProVersion(interpro_version)
    if (version == null) {
        log.error "--interpro <VERSION>: invalid format; expecting number or 'latest'"
        exit 1
    }

    (outdir, error) = InterProScan.resolveDirectory(outdir, false, false)
    if (!outdir) {
        log.error error
        exit 1
    }

    if (outprefix == null) {
        outprefix = "${outdir}/${fasta.split('/').last()}"
    } else if (outprefix.contains("/") || outprefix.contains(File.separator)) {
        log.error "--outprefix must not contain slashes or directory names. Use --outdir to control output location."
        exit 1
    } else {
        outprefix = "${outdir}/${outprefix}"
    }

    emit:
    fasta            // str: path to input fasta file
    apps             // list: list of application to
    datadir          // str: path to data directory, or null if not needed
    outprefix        // str: base path for output files
    formats          // set<String>: output file formats
    version          // str: InterPro version (or "latest")
}
