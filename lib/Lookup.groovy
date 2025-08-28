class Lookup {
    static prepareLookup(List<String> apps, String url, String interproscan_version, workflow_manifest) {
        List<String> allMatchesApiApps = []
        List<String> matchesApiApps = []
        List<String> localOnlyApps  = []
        String matchesUrl
        String apiVersion
        String error
        
        Map info = HTTPRequest.fetch("${HTTPRequest.sanitizeURL(url)}/info".toString(), null, 0, true)
        if (info == null) {
            error = "An error occurred while querying the Matches API [/info];" +
                    " analyses will be run locally"
        } else {
            apiVersion = info.api ?: "X.Y.Z"
            def majorVersion = apiVersion.split("\\.")[0]
            if (majorVersion != "0") {
                error = "${workflow_manifest.name} ${workflow_manifest.version}" +
                        " is not compatible with the Matches API at ${url};" +
                        " analyses will be run locally"
            } else {
                if (info.analyses) {
                    allMatchesApiApps.addAll(info.analyses*.name.collect { it.toLowerCase().replaceAll("[-\\s]", "") })
                    apps.each { app ->
                        (allMatchesApiApps.contains(app.toLowerCase().replaceAll("[-\\s]", "")) ? matchesApiApps : localOnlyApps) << app
                    }
                } else {
                    error = "Could not retrieve the list of available analyses from the Matches API; " +
                            " analyses will be run locally"
                }
            }
        }

        // If there was an error in the PREPARE_LOOKUP step, or no applications are available via the API,
        // skip the API querying and write out all sequences to noApi.fasta
        // matchesUrl used as the go-no-go signal in LOOKUP_MATCHES
        if (!error && !matchesApiApps.isEmpty() ) {
            matchesUrl = url
        } else {
            matchesUrl = null
            localOnlyApps  = apps.clone() as List<String>
        }

        return [matchesApiApps, localOnlyApps, matchesUrl, apiVersion, error]
    }
}