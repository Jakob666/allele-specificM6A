import org.apache.commons.math3.distribution.NormalDistribution;
import org.junit.Test;

public class NormalDistributionApp {

    @Test
    public void normTest(){
        NormalDistribution nd = new NormalDistribution(0.177, 1.44);
        System.out.println(nd.cumulativeProbability(0));
    }

}
