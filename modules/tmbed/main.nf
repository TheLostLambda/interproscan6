import groovy.json.JsonOutput

import Match

process RUN_TMBED_CPU {
    label 'medium', 'tmbed_container'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("tmbed.pred")

    script:
    """
    tmbed predict
        -f ${fasta} \
        -m /t5 \
        -p tmbed.pred
    """
}

process RUN_TMBED_GPU {
    label 'medium', 'tmbed_container', 'use_gpu'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("tmbed_output")

    script:
    """
    tmbed predict
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

    def startNewMatch(String seqId, String tmbed_model, Integer position) {
        def (modelAc, modelSig) = MODEL_TYPES[tmbed_model]
        hits[seqId].computeIfAbsent(modelAc) {
            Match match = new Match(modelAc)
            match.signature = modelSig
            match
        }
        return [modelAc, position]
    }
    def endMatch(Integer start, Integer position) {
        match.addLocation( new Location(start, position) )
    }

    Map<String, Match> hits = [:]
    Match match
    String seqMd5  = null
    String modelAc = null
    Boolean seqLine = false

    new File(tmbed_out.toString()).eachLine { line ->
        if (line.startsWith(">")) { // Processing a new protein
            // Store a match for the previous protein
            if (modelAc) {
                endMatch(start, position)
            }
            // Start a new protein
            seqMd5 = line.trim(">").trim()
            modelAc = null
            seqLine = true
            hits.computeIfAbsent(seqMd5) { [:] }
        } else { // Processing a site in the protein
            if (seqLine) {
                seqLine = false
                return
            }
            position = 0
            for (character in line) { // SSSSSSSSSSSS.........BBBBBBBBBB.................bbbb
                position ++
                if (character != ".") {  // We have a hit!
                    if (!modelAc) { // Found the start of a new hit
                        (modelAc, start) = startNewMatch(seqMd5, character, position)
                    } else if (modelAc && modelAc != MODEL_TYPES[character][0]) { // Found a new and different hit
                        endMatch(start, position)  // Process the previous hit
                        (modelAc, start) = startNewMatch(seqMd5, character, position)
                    }
                    // ELSE (character != '.' and modelAc == MODEL_TYPES[character][0]) -> processing the same match
                } else {
                    if (modelAc) {
                        endMatch(start, position) // Process the previous hit
                        modelAc = null // Not currently parsing a hit
                    }
                }
            }
        }
        if (modelAc && seqMd5) {  // Add the last match
            endMatch(start, position)
        }
    }

    def outputFilePath = task.workDir.resolve("tmbed.json")
    def json = JsonOutput.toJson(hits)
    new File(outputFilePath.toString()).write(json)
}