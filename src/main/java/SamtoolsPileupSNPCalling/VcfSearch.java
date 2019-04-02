package SamtoolsPileupSNPCalling;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

public class VcfSearch {

    private File vcfFile;
    private Logger log;

    public VcfSearch(File vcfFilePath, Logger logger) {
        this.vcfFile = vcfFilePath;
        this.log = logger;
    }

    public HashMap<String, ArrayList<Integer>> vcfList() {

        HashMap<String, ArrayList<Integer>> vcfSearchTree = new HashMap<String, ArrayList<Integer>>();
        try {
            BufferedReader bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.vcfFile))
            );

            String line = "";
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;

                    String[] lineInfo = line.split("\t");
                    String chrNum = lineInfo[0];
                    int position = Integer.parseInt(lineInfo[1]);
                    ArrayList<Integer> vcfPos = vcfSearchTree.getOrDefault(chrNum, new ArrayList<Integer>());
                    if (vcfPos.size() != 0) {
                        if (position != vcfPos.get(vcfPos.size() - 1))
                            vcfPos.add(position);
                    }else
                        vcfPos.add(position);
                    vcfSearchTree.put(chrNum, vcfPos);
                }
            }
            bfr.close();
        } catch (IOException ie) {
            this.log.error(ie.getMessage());
            System.exit(2);
        }

        return vcfSearchTree;
    }

    public boolean binarySearch(ArrayList<Integer> arr, int position) {
        int len = arr.size();
        int start = 0;
        int end = len - 1;
        int mid;
        while (start <= end) {
            mid = start + (end - start) / 2;
            int cmp = Integer.compare(position, arr.get(mid));
            if (cmp == 0)
                return true;
            else if (cmp > 0){
                start = mid + 1;
            } else {
                end = mid - 1;
            }
        }
        return false;
    }
}
