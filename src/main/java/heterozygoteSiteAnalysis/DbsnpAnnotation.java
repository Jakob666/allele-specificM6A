package heterozygoteSiteAnalysis;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

/**
 * annotate WES SNP calling result with dbsnp
 */
public class DbsnpAnnotation {
    private String dbsnpFile;
    private HashMap<String, LinkedList<DIYNode>> dbsnpRecord = new HashMap<>();
    private Logger logger;

    public DbsnpAnnotation(String dbsnpFile, Logger logger) {
        this.dbsnpFile = dbsnpFile;
        this.logger = logger;
    }

    public void parseDbsnpFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.dbsnpFile))));
            String line = "", chrNum, position;
            String[] info;
            int pos;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    position = info[1];
                    pos = Integer.valueOf(position);
                    if (!this.dbsnpRecord.containsKey(chrNum))
                        logger.debug("parsing chromosome " + chrNum);
                    LinkedList<DIYNode> chrMutations = this.dbsnpRecord.getOrDefault(chrNum, new LinkedList<>());
                    if (chrMutations.size() != 0) {
                        DIYNode lastNode = chrMutations.getLast();
                        if (lastNode.records.size() < 1000000)
                            lastNode.records.add(position);
                    } else {
                        DIYNode node = new DIYNode(pos);
                        node.records.add(position);
                        chrMutations.add(node);
                    }
                    this.dbsnpRecord.put(chrNum, chrMutations);
                    info = null;
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

    public class DIYNode {
        Integer startValue;
        HashSet<String> records;
        DIYNode(Integer startValue) {
            this.startValue = startValue;
            // default 1M record in one node
            this.records = new HashSet<>(1000000);
        }
    }

    public HashMap<String, LinkedList<DIYNode>> getDbsnpRecord() {
        return this.dbsnpRecord;
    }

    public static boolean getSearchRes(HashMap<String, LinkedList<DIYNode>> dbsnpRecord, String chrNum, String position) {
        if (!dbsnpRecord.containsKey(chrNum))
            return false;
        LinkedList<DIYNode> nodes = dbsnpRecord.get(chrNum);
        int pos = Integer.valueOf(position);
        for (DIYNode node: nodes) {
            if (node.startValue > pos)
                break;
            if (node.records.contains(position))
                return true;
        }
        return false;
    }
}
