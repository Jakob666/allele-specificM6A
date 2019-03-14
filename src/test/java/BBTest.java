import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.random.RandomData;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.junit.Test;

import java.io.FileWriter;
import java.io.IOException;

public class BBTest {

    @Test
    public void samplingTest() {
        int totalSamples = 5000;

        int[] majorSNPReads = new int[totalSamples];
        int[] minorSNPReads = new int[totalSamples];
        RandomDataGenerator rdg = new RandomDataGenerator();
        double theta = 0.6;
        double tao = 0.005;

        int totalReads;
        int[] samplingRes;

        for (int i=0; i < totalSamples; i++) {
            totalReads = rdg.nextInt(100, 200);
            samplingRes = sampling(1, totalReads, theta, tao);
            majorSNPReads[i] = samplingRes[0];
            minorSNPReads[i] = totalReads - samplingRes[0];
        }

        try {
            FileWriter fw = new FileWriter("C:\\Users\\hbs\\Desktop\\SNPReader.tsv");
            for (int i = 0; i < totalSamples; i++) {
                String readsCount = majorSNPReads[i] + "\t" + minorSNPReads[i] + "\n";
                fw.write(readsCount);
            }
            fw.flush();

            fw.close();
        } catch (IOException ie) {
            ie.printStackTrace();
        }
    }

    public int[] sampling(int sampleSize, int n, double theta, double tao) {
        int[] randomSample = new int[sampleSize];

        UniformRealDistribution ufd = new UniformRealDistribution();
        double[] sampleProba = ufd.sample(sampleSize);
        for (int i = 0; i < sampleSize; i++) {
            randomSample[i] = sampleBinarySearch(n, theta, tao, sampleProba[i]);
        }
        return randomSample;
    }

    public double betaBinomialDistributionProbability(int n, int y, double theta, double tao) {
        double alpha1 = y+theta/tao;
        double beta1 = n-y+(1-theta)/tao;
        double alpha2 = theta/tao;
        double beta2 = (1-theta)/tao;

        double logBetaBinomial = logCombination(n, y) + logBetaFunction(alpha1, beta1) - logBetaFunction(alpha2, beta2);
        return Math.exp(logBetaBinomial);
    }

    @Test
    public void betaBinomialProbaTest() {
        double res = betaBinomialDistributionProbability(10, 2, 0.5, 2);
        System.out.println(res);
    }


    public double betaBinomialCdf(int lowerBound, int upperBound, int n, double theta, double tao) {
        double cdf = 0.0;
        while (lowerBound <= upperBound) {
            cdf += betaBinomialDistributionProbability(n, lowerBound, theta, tao);
            lowerBound++;
        }
        return cdf;
    }

    @Test
    public void betaBinomialCdfTest() {
        double cdf = betaBinomialCdf(0, 5, 10, 0.5, 2);
        System.out.println(cdf);
    }


    private double logBetaFunction(double alpha, double beta){
        return Beta.logBeta(alpha, beta);
    }


    private double logCombination(int total, int part) {
        return CombinatoricsUtils.binomialCoefficientLog(total, part);
    }


    private int sampleBinarySearch(int n, double theta, double tao, double uniformProba) {
        int tmp;
        int lowerBound = 0;
        int upperBound = n;
        int mid = (lowerBound + n) / 2;
        double cdfProbaNMinus1, cdfProbaN;
        final double epsilon = 0.000001;
        while (mid > lowerBound) {
            cdfProbaNMinus1 = betaBinomialCdf(0, mid-1, n, theta, tao);
            cdfProbaN = betaBinomialCdf(0, mid, n, theta, tao);
            if (cdfProbaN-uniformProba >= epsilon && uniformProba-cdfProbaNMinus1 > epsilon) {
                return mid;
            } else if (cdfProbaN-uniformProba < epsilon) {
                tmp = mid;
                mid =  (upperBound + mid) / 2;
                lowerBound = tmp;
            } else if (uniformProba-cdfProbaNMinus1 < epsilon) {
                tmp = mid;
                mid = (lowerBound + mid) / 2;
                upperBound = tmp;
            }
        }
        return mid;
    }

    @Test
    public void sampleBinarySearchTest() {
        // 理论结果为2
        int res = sampleBinarySearch(10, 0.5, 2, 0.5);
        System.out.println(res);
    }
}
