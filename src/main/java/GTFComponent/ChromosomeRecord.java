package GTFComponent;

import java.util.HashMap;

public class ChromosomeRecord {
    private String chrName;
    private HashMap<String, GeneRecord> geneOnChr;

    public ChromosomeRecord(String chrName) {
        this.chrName = chrName;
        this.geneOnChr = new HashMap<>();
    }

    public void renewGeneOnChr(String geneName, GeneRecord geneRecord) {
        this.geneOnChr.put(geneName, geneRecord);
    }

    public String getChrName() {
        return this.chrName;
    }

    public HashMap<String, GeneRecord> getChrGenes() {
        return this.geneOnChr;
    }
}
