import groovy.json.JsonOutput

import Match

process RUN_TMBED_CPU {
    label 'large', 'tmbed_container'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("tmbed.pred")

    script:
    """
    tmbed predict \
        -f ${fasta} \
        -m /opt/tmbed/models \
        -p tmbed.pred
    """
}

process RUN_TMBED_GPU {
    label 'large', 'tmbed_container', 'use_gpu'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("tmbed.pred")

    script:
    """
    tmbed predict \
        -f ${fasta} \
        -m /opt/tmbed/models \
        -p tmbed.pred \
        --use-gpu
    """
}

process PARSE_TMBED {
    label 'tiny'
    executor 'local'

    input:
    tuple val(meta), val(tmbed_out)

    output:
    tuple val(meta), path("tmbed.json")

    exec:
    SignatureLibraryRelease LIBRARY = new SignatureLibraryRelease("TMbed", "1.0.2")
    def MODEL_TYPES = [ // Sig(acc, name, desc, lib, entry)
        "b": new Signature("TMbeta_out-to-in",  "TMbeta_out-to-in",  null, LIBRARY, null),
        "B": new Signature("TMbeta_in-to-out",  "TMbeta-in-to-out",  null, LIBRARY, null),
        "h": new Signature("TMhelix_out-to-in", "TMhelix-out-to-in", null, LIBRARY, null),
        "H": new Signature("TMhelix_in-to-out", "TMhelix-in-to-out", null, LIBRARY, null),
        "S": new Signature("Signal_peptide",    "Signal peptide",    null, LIBRARY, null)
    ]

    def startNewMatch(String symbol) {
        Signature signature = MODEL_TYPES[symbol]
        Match match = hits[seqMd5].computeIfAbsent(signature.accession) {
            Match newMatch = new Match(signature.accession)
            newMatch.signature = modelSig
            newMatch
        }
        return match
    }

    Map<String, Map<String, Match>> hits = [:]
    Match currentMatch = null
    String seqMd5 = null
    def lineCounter = 0
    /* TMbed output:
    >seqId
    sequence
    prediction */
    new File(tmbed_out.toString()).eachLine { line ->
        line = line.trim()
        if (lineCounter % 3 == 0) {
            seqMd5       = line.substring(1)
            currentMatch = null
            hits.computeIfAbsent(seqMd5) { [:] }
        } else if (lineCounter % 3 == 2) {
            line.eachWithIndex { symbol, position ->
                if (character != ".") { // Found a hit
                    if (!currentMatch) {  // Start a new match
                        currentMatch = startNewMatch(symbol)
                        start = position + 1
                    } else if (currentMatch.accession != MODEL_TYPES[symbol].accession) { // Found a new hit
                        currentMatch.addLocation(new Location(start, position))
                        currentMatch = startNewMatch(symbol)
                        start = position + 1
                    }
                    // else, parsing another symbol from the same currentMatch/hit
                } else {
                    currentMatch.addLocation(new Location(start, position))
                    currentMatch = null
                }
            }
        }
        // Add the final match for this sequence
        if (currentMatch) {
            currentMatch.addLocation(new Location(start, position))
        }
        lineCounter += 1
    }

    def outputFilePath = task.workDir.resolve("tmbed.json")
    def json = JsonOutput.toJson(hits)
    new File(outputFilePath.toString()).write(json)
}

