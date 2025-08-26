import groovy.json.JsonSlurper
import java.net.URL
import groovy.json.JsonOutput

process PREPARE_LOOKUP {
    label    'tiny'
    executor 'local'

    input:
    val _url
    val db_releases
    val interproscan_version  // major.minor iprscan version number
    val workflow_manifest

    output:
    tuple val(matchesApiUrl), val(matchesApiApps), val(error)

    exec:
    String _matchesApiUrl = _url  // reassign to avoid 'variable' already used error when logging
    List<String> _matchesApiApps = []
    String _error = null

    // Get MLS metadata: api (version), release, release_date
    Map info = HTTPRequest.fetch("${HTTPRequest.sanitizeURL(_matchesApiUrl)}/info".toString(), null, 0, true)
    if (info == null) {
        _error = "An error occurred while querying the Matches API [/info]; analyses will be run locally"
        matchesApiUrl = null
    } else {
        def apiVersion = info.api ?: "X.Y.Z"
        def majorVersion = apiVersion.split("\\.")[0]
        if (majorVersion != "0") {
            _error = "${workflow_manifest.name} ${workflow_manifest.version}" +
                    " is not compatible with the Matches API at ${_matchesApiUrl};" +
                    " analyses will be run locally"
            matchesApiUrl = null
        } else if (db_releases) {  // can be null if we don't need data for the selected apps (e.g. mobidblite)
            if (db_releases["interpro"]["version"] != info.release) {
                _error = "The local InterPro version does not match the match API release (Local: ${db_releases['interpro']}, Matches API: ${info.release}).\n" +
                        "Pre-calculated matches will not be retrieved, and analyses will be run locally"
                matchesApiUrl = null
            } else {
                if (info?.analyses) {
                    _matchesApiApps.addAll(info.analyses*.name.collect { it.toLowerCase().replaceAll("[-\\s]", "") })
                } else {
                    _error = "Could not retrieve the list of available analyses from the Matches API; analyses will be run locally"
                    matchesApiUrl = null
                }
            }
        }
    }
    matchesApiUrl = _matchesApiUrl
    matchesApiApps = _matchesApiApps
    error = _error
}

process LOOKUP_MATCHES {
    maxForks 1
    label    'tiny'
    executor 'local'

    input:
    tuple val(index), val(fasta), val(applications), val(api_applications), val(url), val(chunk_size), val(max_retries), val(error)

    output:
    tuple val(index), path("calculatedMatches.json")
    tuple val(index), path("noLookup.fasta"), optional: true

    exec:
    String _error = error  // reassign to avoid 'variable' already used error when logging

    def calculatedMatchesPath = task.workDir.resolve("calculatedMatches.json")
    def calculatedMatches = [:]

    def noLookupFastaPath = task.workDir.resolve("noLookup.fasta")
    def noLookupFasta = new StringBuilder()

    Map<String, String> sequences = FastaFile.parse(fasta.toString())  // [md5: sequence]
    def md5List = sequences.keySet().toList().sort()
    def chunks = md5List.collate(chunk_size)

    String baseUrl = HTTPRequest.sanitizeURL(url.toString())
    boolean success = _error ? false : true

    if (success) {
        for (chunk in chunks) {
            String data = JsonOutput.toJson([md5: chunk])
            def response = HTTPRequest.fetch("${baseUrl}/matches", data, max_retries, true)

            if (response != null) {
                response.results.each {
                    String proteinMd5 = it.md5.toUpperCase() // ensure it matches the local seq Db case
                    if (it.found) {
                        calculatedMatches[proteinMd5] = [:]
                        it.matches.each { matchMap ->
                            String library = matchMap.signature.signatureLibraryRelease.library
                            if (library == "MobiDB Lite") {
                                matchMap.signature.signatureLibraryRelease.library = "MobiDB-lite"
                            }
                            
                            String appName = library.toLowerCase().replaceAll("[-\\s]", "")

                            if (applications.contains(appName)) {
                                matchMap = transformMatch(matchMap, sequences[proteinMd5])
                                calculatedMatches[proteinMd5][matchMap.modelAccession] = matchMap
                            }
                        }
                    } else {
                        def seq = sequences[proteinMd5]
                        noLookupFasta.append(">${proteinMd5}\n")
                        noLookupFasta.append("${seq}\n")
                    }
                }
            } else {
                _error = "An error occurred while querying the Matches API; analyses will be run locally"
                success = false
                break
            }
        }
    }

    if (success) {
        List<String> allApps = applications.clone() as List<String>
        List<String> missingApps = allApps.findAll { !api_applications.contains(it) }
        if (missingApps) {
            log.warn "The following applications are not available via the Matches API and will be run locally:\n${missingApps.join(", ")}\n" +
            "Pre-calculated matches will not be retrieved for these applications, and their analyses will be run locally"
            // add writing out fasta for missing apps
        }

        def jsonMatches = JsonOutput.toJson(calculatedMatches)
        new File(calculatedMatchesPath.toString()).write(JsonOutput.toJson(calculatedMatches))
        if (noLookupFasta.length() != 0) { new File(noLookupFastaPath.toString()).write(noLookupFasta.toString()) }
    } else {
        log.warn _error
        // when the connection fails, write out ALL sequences to "noLookup.fasta". Consider all as unknown/novel
        new File(calculatedMatchesPath.toString()).write(JsonOutput.toJson([:]))
        new File(noLookupFastaPath.toString()).write(new File(fasta.toString()).text)
    }
}

def Map transformMatch(Map match, String seq) {
    // * operator - spread contents of a map or collecion into another map or collection
    return [
        *            : match,
        "treegrafter": ["ancestralNodeID": match["annotationNode"]],
        "locations"  : match["locations"].collect { loc ->
            return [
                *          : loc,
                "hmmBounds": loc["hmmBounds"] ? getReverseHmmBounds(loc["hmmBounds"]) : null,
                "fragments": loc["fragments"].collect { tranformFragment(it) },
                "sites"    : loc["sites"] ?: [],
                "targetAlignment": loc["cigarAlignment"] ? decodeAlignment(loc["cigarAlignment"], seq, loc["start"]) : null
            ]
        },
    ]
}

def getReverseHmmBounds(hmmBounds) {
    return [
        "COMPLETE"            : "[]",
        "N_TERMINAL_COMPLETE" : "[.",
        "C_TERMINAL_COMPLETE" : ".]",
        "INCOMPLETE"          : ".."
    ][hmmBounds]
}

def decodeAlignment(cigarAlignment, sequence, startIndex) {
    def targetAlign = new StringBuilder()
    def index = startIndex - 1 // convert from 1-based numbering to 0-based numbering
    def matcher = (cigarAlignment =~ /(\d+)([MID=X])/)
    matcher.each { match ->
        def len = match[1].toInteger()
        def op = match[2]
        switch(op) {
            case 'M':
            case '=':
            case 'X':
                targetAlign << sequence.substring(index, index + len)
                index += len
                break
            case 'I':
                targetAlign << sequence.substring(index, index + len).toLowerCase()
                index += len
                break
            case 'D':
                targetAlign << '-' * len
                break
            default:
                throw new IllegalArgumentException("Unsupported CIGAR operation: $op")
        }
    }
    return targetAlign.toString()
}

def Map tranformFragment(Map fragment) {
    return [
        "start"   : fragment["start"],
        "end"     : fragment["end"],
        "dcStatus": fragment["type"]
    ]
}
