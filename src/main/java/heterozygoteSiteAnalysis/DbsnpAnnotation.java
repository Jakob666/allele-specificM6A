package heterozygoteSiteAnalysis;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;

/**
 * annotate WES SNP calling result with dbsnp
 */
public class DbsnpAnnotation {
    private String dbsnpFile;
    private HashMap<String, HashSet<String>> dbsnpRecord = new HashMap<>();

    public DbsnpAnnotation(String dbsnpFile) {
        this.dbsnpFile = dbsnpFile;
    }

    public void parseDbsnpFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.dbsnpFile))));
            String line = "", chrNum, position;
            String[] info;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    position = info[1];

                    HashSet<String> chrMutations = this.dbsnpRecord.getOrDefault(chrNum, new HashSet<>());
                    chrMutations.add(position);
                    this.dbsnpRecord.put(chrNum, chrMutations);
                }
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException ie) {
                    ie.printStackTrace();
                }
            }
        }
    }

    public HashMap<String, HashSet<String>> getDbsnpRecord() {
        return this.dbsnpRecord;
    }
}
