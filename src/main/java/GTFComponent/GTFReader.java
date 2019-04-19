package GTFComponent;

import java.io.*;
import java.util.HashMap;

public class GTFReader {
    private HashMap<String, ChromosomeRecord> chromosomeMap = new HashMap<>();

    /**
     * parse the 8th columns of lines GTF file
     * @param attributes the 8th column's content
     * @return GeneAttribute instance
     */
    public GeneAttribute parseAttributes(String attributes) {
        String[] strArr = attributes.split("\";");
        String geneID = "Unknown";
        String geneName = "Unknown";
        String transcriptID = "Unknown";
        String transcriptName = "Unknown";
        String geneBioType = "Unknown";
        for(int i=0; i<strArr.length; i++) {
            strArr[i] = strArr[i].trim();
            if(strArr[i].startsWith("gene_id")) {
                geneID = strArr[i].substring(strArr[i].indexOf("\"") + 1, strArr[i].length());
            }
            if(strArr[i].startsWith("gene_name")) {
                geneName = strArr[i].substring(strArr[i].indexOf("\"") + 1, strArr[i].length());
            }
            if(strArr[i].startsWith("transcript_id")) {
                transcriptID = strArr[i].substring(strArr[i].indexOf("\"") + 1, strArr[i].length());
            }
            if(strArr[i].startsWith("transcript_name")) {
                transcriptName = strArr[i].substring(strArr[i].indexOf("\"") + 1, strArr[i].length());
            }
            if(strArr[i].startsWith("gene_biotype")) {
                geneBioType = strArr[i].substring(strArr[i].indexOf("\"") + 1, strArr[i].length());
            }
        }
        GeneAttribute geneAttr = new GeneAttribute();
        geneAttr.setGeneID(geneID);
        if(geneName.equals("Unknown") && (!geneID.equals("Unknown")) )
            geneAttr.setGeneName(geneID);
        else
            geneAttr.setGeneName(geneName);
        geneAttr.setTranscriptID(transcriptID);
        if(transcriptName.equals("Unknown") && (!transcriptID.equals("Unknown")) )
            geneAttr.setTranscriptName(transcriptID);
        else
            geneAttr.setTranscriptName(transcriptName);
        geneAttr.setBioType(geneBioType);

        return geneAttr;
    }

    /**
     * get information from GTF file to form chromosomeMap
     * @param gtfFile GTF file path
     */
    public void readFromFile(String gtfFile) {
        BufferedReader br = null;
        try
        {
            this.chromosomeMap.clear();
            br = new BufferedReader(
                    new InputStreamReader(new FileInputStream(new File(gtfFile)))
            );
            String strLine = "";
            String[] strArr;
            String chrNum, recordType;
            GeneAttribute geneAttr = null;
            ChromosomeRecord chrRec;
            GeneRecord geneRec;
            TranscriptRecord transcriptRec;
            ElementRecord elementRec;

            while(strLine != null) {
                strLine = br.readLine();
                if(strLine != null && !strLine.startsWith("#")) {
                    strArr = strLine.split("\t");
                    geneAttr = parseAttributes(strArr[8]);
                    // only simulate reads from protein coding genes
                    if (!geneAttr.getBioType().equals("protein_coding"))
                        continue;
                    chrNum = strArr[0];
                    recordType = strArr[2];

                    chrRec = this.chromosomeMap.getOrDefault(chrNum, new ChromosomeRecord(chrNum));

                    if(recordType.equals("gene")) {
                        geneRec = new GeneRecord();
                        geneRec.setGeneRange(Integer.parseInt(strArr[3]), Integer.parseInt(strArr[4]));
                        geneRec.setGeneId(geneAttr.getGeneID());
                        geneRec.setGeneName(geneAttr.getGeneName());
                        geneRec.setBioType(geneAttr.getBioType());
                        geneRec.setStrand(strArr[6]);
                        chrRec.renewGeneOnChr(geneAttr.getGeneID(), geneRec);
                    } else if(recordType.equals("transcript")) {
                        geneRec = chrRec.getChrGenes().get(geneAttr.getGeneID());
                        transcriptRec = geneRec.getTranscriptIsoform().getOrDefault(geneAttr.getTranscriptID(), new TranscriptRecord());
                        transcriptRec.setTranscriptRange(Integer.parseInt(strArr[3]), Integer.parseInt(strArr[4]));
                        transcriptRec.setTranscriptId(geneAttr.getTranscriptID());
                        transcriptRec.setTranscriptName(geneAttr.getTranscriptName());
                        transcriptRec.setBioType(geneAttr.getBioType());
                        transcriptRec.setStrand(strArr[6]);
                        geneRec.renewTranscriptIsoform(geneAttr.getTranscriptID(), transcriptRec);
                    } else {
                        transcriptRec = chrRec.getChrGenes().get(geneAttr.getGeneID()).getTranscriptIsoform()
                                        .get(geneAttr.getTranscriptID());
                        elementRec = new ElementRecord(strArr[2], strArr[6], Integer.parseInt(strArr[3]),
                                                                     Integer.parseInt(strArr[4]));
                        if(!strArr[7].equals("."))
                            elementRec.setFrame(Integer.parseInt(strArr[7]));
                        else
                            elementRec.setFrame(0);
                        transcriptRec.renewElementList(recordType, elementRec);
                    }

                    this.chromosomeMap.put(chrNum, chrRec);
                }
            }
            br.close();
        } catch(IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    public HashMap<String, ChromosomeRecord> getChromosomeMap()
    {
        return chromosomeMap;
    }
}
