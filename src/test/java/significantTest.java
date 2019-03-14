import betaBinomialMetaAnalysis.SignificantTest;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Set;

public class significantTest {
    @Test
    public void test() throws Exception{
        double test_rho = 0.0048841870010322315;
        int[][] reads = readsCount();
//        int[][] peakReads = groupIntoPeak(reads);
        // 假设前5个SNP被一个peak覆盖
        int peakNum = 45;
        int[] peakMajorReads = new int[peakNum];
        int[] peakMinorReads = new int[peakNum];
        for (int i = 0; i < peakNum; i++) {
            peakMajorReads[i] = reads[i][0];
            peakMinorReads[i] = reads[i][1];
        }

        SignificantTest st = new SignificantTest(peakMajorReads, peakMinorReads, test_rho);
        double pValue = st.testSignificant();
        System.out.println(pValue);
    }

    private int[][] readsCount() throws Exception {
        BufferedReader bf = new BufferedReader(
                new InputStreamReader(new FileInputStream(new File("C:\\Users\\hbs\\Desktop\\SNPReader.tsv")))
        );
        int[][] readsCount = new int[5000][2];

        String line = bf.readLine();
        int major, minor;
        int ord = 0;
        String[] lineSplit;

        while (line != null) {
            if (!line.isEmpty()) {
                lineSplit = line.split("\t");
                major = Integer.parseInt(lineSplit[0]);
                minor = Integer.parseInt(lineSplit[1]);

                readsCount[ord][0] = major;
                readsCount[ord][1] = minor;
                ord++;
            }
            line = bf.readLine();
        }
        return readsCount;
    }

    private int[][] groupIntoPeak(int[][] readsCount) {
        RandomDataGenerator rdg = new RandomDataGenerator();
        int sum = 0;
        int randInt;
        HashMap<Integer, Integer> idx = new HashMap<>();
        int[][] peakReads;
        while (sum <= 5000) {
            randInt = rdg.nextInt(1, 10);
            if (sum + randInt <= 5000) {
                idx.put(sum, sum + randInt -1);
            } else {
                idx.put(sum, 4999);
            }
            sum += randInt;
        }

        Set<Integer> ks = idx.keySet();
        peakReads = new int[ks.size()][2];
        int ord = 0;
        for (Integer k : ks) {
            int majorReads = 0;
            int minorReads = 0;
            int v = idx.get(k);
            for (int i = k; i <= v; i++) {
                majorReads += readsCount[i][0];
                minorReads += readsCount[i][1];
            }
            peakReads[ord][0] = majorReads;
            peakReads[ord][1] = minorReads;
            ord += 1;
        }
        return peakReads;
    }
}
