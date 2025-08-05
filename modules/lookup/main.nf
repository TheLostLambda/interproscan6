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
    val apps

    output:
    val matchesApiApps

    exec:
    def _matchesApiApps   = []    // apps in the matches api
    String _matchesApiUrl = _url  // reassign to avoid 'variable' already used error when logging
    // Get MLS metadata: api (version), release, release_date
    def _info = HTTPRequest.fetch("${HTTPRequest.sanitizeURL(_matchesApiUrl)}/info".toString(), null, 0, true)
    if (_info == null) {
        log.warn "An error occurred while querying the Matches API; analyses will be run locally"
        matchesApiUrl = null
    } else {
        def apiVersion = _info.api ?: "X.Y.Z"
        def majorVersion = apiVersion.split("\\.")[0]
        if (majorVersion != "0") {
            log.warn "${workflow_manifest.name} ${workflow_manifest.version}" +
                    " is not compatible with the Matches API at ${_matchesApiUrl};" +
                    " analyses will be run locally"
            matchesApiUrl = null
        } else if (db_releases) {  // can be null if we don't need data for the selected apps (e.g. mobidblite)
            if (db_releases["interpro"]["version"] != _info.release) {
                log.warn "The local InterPro version does not match the match API release (Local: ${db_releases['interpro']}, Matches API: ${_info.release}).\n" +
                        "Pre-calculated matches will not be retrieved, and analyses will be run locally"
                matchesApiUrl = null
            } else {
                for (analyse in _info.analyses) {
                    name = analyse.name.toString().toLowerCase().replace(' ', '').replace('-', '')
                    _matchesApiApps << name
                }
            }
        }
    }
    List<String> allApps = apps.clone() as List<String>
    List<String> _missing_apps = allApps.findAll { !(_matchesApiApps.contains(it)) }
    if (_missing_apps) {
        log.warn "The following applications are not available in the Matches API: ${_missing_apps.join(", ")}.\n" +
                 "Pre-calculated matches will not be retrieved for these applications, and analyses will be run locally."
    }
    matchesApiApps = _matchesApiApps ? [_matchesApiApps] : null
}

process LOOKUP_MATCHES {
    maxForks 1
    label    'tiny'
    executor 'local'

    input:
    tuple val(index), val(fasta), val(url), val(applications), val(api_apps), val(chunkSize), val(maxRetries)

    output:
    tuple val(index), path("calculatedMatches.json")
    tuple val(index), path("noMatches.fasta"), optional: true
    tuple val(index), path("noLookup.fasta"), optional: true

    exec:
    def calculatedMatchesPath = task.workDir.resolve("calculatedMatches.json")
    def calculatedMatches = [:]

    def noMatchesFastaPath = task.workDir.resolve("noMatches.fasta")
    def noMatchesFasta = new StringBuilder()

    def noLookupFastaPath = task.workDir.resolve("noLookup.fasta")

    // Check for apps who are not listed in the matches API
    // We will need a FASTA file with all sequences if some apps are not in the API
    List<String> allApps = applications.clone() as List<String>
    List<String> _all_api_apps = api_apps.clone() as List<String>
    List<String> _missing_apps = allApps.findAll { !(_all_api_apps.contains(it)) }
    if (_missing_apps) {
        def sourceFasta = new File(fasta.toString())
        def noLookupFasta = new File(noLookupFastaPath.toString())
        noLookupFasta.text = sourceFasta.text
    }

    Map<String, String> sequences = FastaFile.parse(fasta.toString())  // [md5: sequence]
    def md5List = sequences.keySet().toList().sort()
    def chunks = md5List.collate(chunkSize)

    String baseUrl = HTTPRequest.sanitizeURL(url.toString())
    boolean success = true
    for (chunk in chunks) {
        String data = JsonOutput.toJson([md5: chunk])
        def response = HTTPRequest.fetch("${baseUrl}/matches", data, maxRetries, true)

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
                    noMatchesFasta.append(">${proteinMd5}\n")
                    noMatchesFasta.append("${seq}\n")
                }
            }
        } else {
            success = false
            break
        }
    }

    if (success) {
        def jsonMatches = JsonOutput.toJson(calculatedMatches)
        new File(calculatedMatchesPath.toString()).write(JsonOutput.toJson(calculatedMatches))
        if (noMatchesFasta.length() != 0) { new File(noMatchesFastaPath.toString()).write(noMatchesFasta.toString()) }
    } else {
        log.warn "An error occurred while querying the Matches API, analyses will be run locally"
        // when the connection fails, write out all sequences to "noMatches.fasta"
        new File(calculatedMatchesPath.toString()).write(JsonOutput.toJson([:]))
        if (noMatchesFasta.length() != 0) { new File(noMatchesFastaPath.toString()).write(noMatchesFasta.toString()) }
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
