import com.fasterxml.jackson.databind.ObjectMapper
import java.time.format.DateTimeFormatter
import java.time.LocalDate

process WRITE_TSV {
    label    'tiny'
    executor 'local'

    input:
    val matchesFiles
    val output_file
    val seqDbPath
    val nucleic

    exec:
    SeqDB db = new SeqDB(seqDbPath.toString())
    def tsvFile = new File(output_file)
    tsvFile.text = "" // clear the file if it already exists

    // Each line contains: seqId md5 seqLength memberDb modelAcc sigDesc start end evalue status date entryAcc entryDesc goterms pathways
    def currentDate = LocalDate.now().format(DateTimeFormatter.ofPattern("dd-MM-yyyy"))
    Set<String> seenNucleicMd5s = new HashSet<>()

    matchesFiles.each { matchFile ->
        Map proteins = new ObjectMapper().readValue(new File(matchFile.toString()), Map)
        if (nucleic) {
            nucleicToProteinMd5 = db.groupProteins(proteins)
            nucleicToProteinMd5.each { String nucleicMd5, Set<String> proteinMd5s ->
                if (!seenNucleicMd5s.contains(nucleicMd5)) {
                    seenNucleicMd5s.add(nucleicMd5)

                    seqData = db.nucleicMd5ToNucleicSeq(nucleicMd5)


                    // TODO
                    seqData.each { row ->
                        String seqId = row.id
                        int seqLength = row.sequence.trim().length()
                        match.locations.each { Location loc ->
                            def line = formatLine(seqId, proteinMd5, seqLength, match, loc, memberDb, sigDesc, currentDate, entryAcc, entryDesc, goterms, pathways)
                            tsvFile.append("${line}\n")
                        }
                    }
                }
            }
        } else {
            proteins.each { String proteinMd5, Map proteinMatches ->
                proteinMatches.each { modelAcc, match ->
                    match = Match.fromMap(match)
                    String memberDb = match.signature.signatureLibraryRelease.library
                    String sigDesc = match.signature.description ?: '-'
                    def goterms = match.signature.entry?.goXRefs
                    def pathways = match.signature.entry?.pathwayXRefs
                    String entryAcc = match.signature.entry?.accession ?: '-'
                    String entryDesc = match.signature.entry?.description ?: '-'

                    seqData = db.proteinMd5ToProteinSeq(proteinMd5)

                    seqData.each { row ->
                        String seqId = row.id
                        int seqLength = row.sequence.trim().length()
                        match.locations.each { Location loc ->
                            def line = formatLine(seqId, proteinMd5, seqLength, match, loc, memberDb, sigDesc, currentDate, entryAcc, entryDesc, goterms, pathways)
                            tsvFile.append("${line}\n")
                        }
                    }
                }
            }
        }
        proteins.each { String proteinMd5, Map matchesMap ->
            matchesMap.each { modelAcc, match ->
                match = Match.fromMap(match)
                String memberDb = match.signature.signatureLibraryRelease.library
                String sigDesc = match.signature.description ?: '-'
                def goterms = match.signature.entry?.goXRefs
                def pathways = match.signature.entry?.pathwayXRefs
                String entryAcc = match.signature.entry?.accession ?: '-'
                String entryDesc = match.signature.entry?.description ?: '-'



                seqData = nucleic ? db.proteinMd5ToNucleicSeq(proteinMd5) : db.proteinMd5ToProteinSeq(proteinMd5)
                seqData.each { row ->  // Protein or Nucleic: [id, desc, sequence]
                    String seqId = nucleic ? "${row.nid}_${row.pid}" : row.id
                    int seqLength = row.sequence.trim().length()
                    match.locations.each { Location loc ->
                        def line = formatLine(seqId, proteinMd5, seqLength, match, loc, memberDb, sigDesc, currentDate, entryAcc, entryDesc, goterms, pathways)
                        tsvFile.append("${line}\n")
                    }
                }
            } // end of matches in matchesNode
        } // end of proteins.each
    } // end of matchesFiles
}

def formatLine(seqId, 
               md5, 
               seqLength, 
               match, 
               loc, 
               memberDb, 
               sigDesc, 
               currentDate, 
               entryAcc, 
               entryDesc, 
               interproGoTerms, 
               interproPathways) {
    int start = loc.start
    int end = loc.end
    def scoringValue = "-"
    def pantherGoTerms = []
    switch (memberDb) {
        case ["CDD", "PRINT"]:
            scoringValue = match.evalue
            break
        case ["SignalP-Prok", "SignalP-Euk"]:
            scoringValue = loc.pvalue
            break
        case ["HAMAP", "PROSITE profiles"]:
            scoringValue = loc.score
            break
        case ["COILS", "MobiDB-lite", "Phobius", "PROSITE patterns", "DeepTMHMM"]:
            scoringValue = "-"
            break
        case "PANTHER":
            pantherGoTerms = match.treegrafter.goXRefs.collect { "${it.id}(PANTHER)" }
            scoringValue = loc.evalue
            break
        default:
            scoringValue = loc.evalue
            break
    }

    goTerms = interproGoTerms.collect { "${it.id}(InterPro)" } + pantherGoTerms
    return [
        seqId, 
        md5,
        seqLength,
        memberDb,
        match.signature.accession,
        sigDesc,
        start,
        end,
        scoringValue,
        "T",
        currentDate,
        entryAcc,
        entryDesc,
        goTerms.join("|") ?: "-",
        interproPathways.collect { "${it.databaseName}:${it.id}" }.join("|") ?: "-"
    ].join("\t")
}