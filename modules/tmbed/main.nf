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
        -m /t5 \
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
        -m /t5 \
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
    def MODEL_TYPES = [
        "b": ["TMbed-b", new Signature("TMbed-b", "TMbeta_out-to-in", "Transmembrane beta strand (OUT->IN orientation)", LIBRARY, null)],
        "B": ["TMbed-B", new Signature("TMbed-B", "TMbeta-in-to-out", "Transmembrane beta strand (IN->OUT orientation)", LIBRARY, null)],
        "h": ["TMbed-h", new Signature("TMbed-h", "TMhelix-out-to-in", "Transmembrane alpha helix (OUT->IN orientation)", LIBRARY, null)],
        "H": ["TMbed-H", new Signature("TMbed-H", "TMhelix-in-to-out", "Transmembrane alpha helix (IN->OUT orientation)", LIBRARY, null)],
        "S": ["TMbed-S", new Signature("TMbed-S", "Signal_peptide", "Signal peptide", LIBRARY, null)]
    ]

    Map<String, Map<String, Match>> hits = [:]
    Match currentMatch = null
    String seqMd5 = null
    String modelAc = null
    Integer start = null
    Integer position = null
    Boolean seqLine = false

    new File(tmbed_out.toString()).eachLine { line ->
        if (line.startsWith(">")) { // Processing a new protein
            // Store a match for the previous protein
            if (modelAc && currentMatch && start != null) {
                currentMatch.addLocation(new Location(start, position))
            }
            // Start a new protein
            seqMd5 = line.replace(">", "").trim()
            modelAc = null
            currentMatch = null
            start = null
            seqLine = true
            hits.computeIfAbsent(seqMd5) { [:] }
        } else { // Processing a site in the protein
            if (seqLine) {
                seqLine = false
                return
            }
            position = 0
            for (character in line) { // SSSSSSSSSSSS.........BBBBBBBBBB.................bbbb
                position++
                if (character != ".") {  // We have a hit!
                    if (!modelAc) { // Found the start of a new hit
                        (hits, modelAc, start, currentMatch) = startNewMatch(MODEL_TYPES, hits, seqMd5, character, position)
                    } else if (modelAc && modelAc != MODEL_TYPES[character][0]) { // Found a new and different hit
                        currentMatch.addLocation(new Location(start, position - 1))  // Process the previous hit
                        (hits, modelAc, start, currentMatch) = startNewMatch(MODEL_TYPES, hits, seqMd5, character, position)
                    }
                    // ELSE (character != '.' and modelAc == MODEL_TYPES[character][0]) -> processing the same match
                } else {
                    if (modelAc && currentMatch) {
                        currentMatch.addLocation(new Location(start, position - 1)) // Process the previous hit
                        modelAc = null // Not currently parsing a hit
                        currentMatch = null
                        start = null
                    }
                }
            }
        }
    }

    // Add the last match if there is one
    if (modelAc && currentMatch && start != null && seqMd5) {
        currentMatch.addLocation(new Location(start, position))
    }

    def outputFilePath = task.workDir.resolve("tmbed.json")
    def json = JsonOutput.toJson(hits)
    new File(outputFilePath.toString()).write(json)
}

def startNewMatch(Map model_type, Map hits, String seq_id, String tmbed_model, Integer position) {
    def (modelAc, modelSig) = model_type[tmbed_model]
    def match = hits[seq_id].computeIfAbsent(modelAc) {
        Match newMatch = new Match(modelAc)
        newMatch.signature = modelSig
        newMatch
    }
    return [hits, modelAc, position, match]
}
