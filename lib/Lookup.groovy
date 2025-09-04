class Lookup {
    static prepareLookup(List<String> apps, Boolean no_matches_api, String url, String interproscan_version, workflow_manifest) {
        List<String> allMatchesApiApps = []
        List<String> matchesApiApps = []
        List<String> localOnlyApps  = []
        String apiInterproVersion = null
        String apiVersion
        String error

        if (!no_matches_api && url) { // url check for unit tests
            Map info = HTTPRequest.fetch("${HTTPRequest.sanitizeURL(url)}/info".toString(), null, 0, true)
            if (info == null) {
                error = "An error occurred while querying the Matches API;" +
                        " analyses will be run locally"
            } else {
                apiVersion = info.api ?: "X.Y.Z"
                def majorVersion = apiVersion.split("\\.")[0]
                if (majorVersion != "0") {
                    error = "${workflow_manifest.name} ${workflow_manifest.version}" +
                            " is not compatible with the Matches API at ${url};" +
                            " analyses will be run locally"
                } else {
                    apiInterproVersion = info.release
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
        }

        if (error || matchesApiApps.isEmpty() || no_matches_api) {
            localOnlyApps  = apps.clone() as List<String>
        }

        return [matchesApiApps, localOnlyApps, apiInterproVersion, error]
    }
}
