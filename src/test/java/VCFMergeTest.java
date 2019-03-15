import heterozygoteSiteAnalysis.MergeVCFFile;
import org.junit.Test;

import java.io.IOException;

public class VCFMergeTest {
    private MergeVCFFile mvf;

    @Test
    public void vcfMergeTest() throws IOException {
        this.mvf = new MergeVCFFile("src/test/testResource", "src/test/testResource/vcfMergeOutput.vcf");
        this.mvf.mergeVCF();
    }
}
