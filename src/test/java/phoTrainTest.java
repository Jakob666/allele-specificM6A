import betaBinomialMetaAnalysis.LogLikelihoodFunc;
import betaBinomialMetaAnalysis.RhoEstimator;
import org.junit.Test;

import java.io.*;

public class phoTrainTest {

    @Test
    public void testEstimate() {
        int[] majorReads, minorReads;
        try {
            int[][] reads = readsCount();
            int sampleNum = reads.length;
            majorReads = new int[sampleNum];
            minorReads = new int[sampleNum];

            for (int i=0; i<sampleNum; i++) {
                majorReads[i] = reads[i][0];
                minorReads[i] = reads[i][1];
            }

            RhoEstimator re = new RhoEstimator(majorReads, minorReads, 0.001, 0.0000001);
            // the parameter pho is quite small, the initial value recommend a tiny value
            double estimateRho = re.gradientAscend(0.05);
            System.out.println(estimateRho);
        } catch (Exception e) {
            e.printStackTrace();
        }
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

    @Test
    public void llfTest() {
        int[] majorReads, minorReads;
        try {
            int[][] reads = readsCount();
            int sampleNum = reads.length;
            majorReads = new int[sampleNum];
            minorReads = new int[sampleNum];

            for (int i=0; i<sampleNum; i++) {
                majorReads[i] = reads[i][0];
                minorReads[i] = reads[i][1];
            }

            LogLikelihoodFunc llf = new LogLikelihoodFunc(majorReads, minorReads);
            double res = llf.logLikelihoodFunc(-155);
//            if (res - Double.NEGATIVE_INFINITY > 0.00001) {
//                System.out.println("hhh");
            System.out.println(res);
//            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
