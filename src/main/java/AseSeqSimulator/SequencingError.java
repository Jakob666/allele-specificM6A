package AseSeqSimulator;

import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

/**
 * reads sequencing error in PCR procedure. A uniform error model
 */
public class SequencingError {

    private double sequencingError;
    private int readLength;
    private UniformRealDistribution uid;
    private ArrayList<String> bases = new ArrayList<>(Arrays.asList("A", "T", "C", "G"));

    public SequencingError(double sequencingError, int readLength) {
        this.sequencingError = sequencingError;
        this.readLength = readLength;
        this.uid = new UniformRealDistribution(0.0, 1.0);
    }

    /**
     * get sequence error positions on PCR reads
     * @return ArrayList contains sequence error position
     */
    public ArrayList<Integer> simulateSequenceError() {
        ArrayList<Integer> sequenceErrorPositions = new ArrayList<>();
        double errorProb;
        for (int i = 0; i < this.readLength; i++) {
            errorProb = this.uid.sample();
            if (errorProb <= this.sequencingError)
                sequenceErrorPositions.add(i);
        }

        return sequenceErrorPositions;
    }

    /**
     * generate read with PCR sequence error base
     * @param originalReads original reads
     * @return sequencing error base
     */
    public String pcrErrorReads(String originalReads) {
        ArrayList<Integer> pcrErrorPositions = this.simulateSequenceError();
        String originalBase, errorBase, mutateRead = null;
        if (pcrErrorPositions.size() == 0)
            mutateRead = originalReads;
        else {
            for (Integer position: pcrErrorPositions) {
                originalBase = originalReads.substring(position, position+1);
                errorBase = originalBase;
                while (errorBase.equals(originalBase)) {
                    Collections.shuffle(this.bases);
                    errorBase = this.bases.get(0);
                }
                mutateRead = originalReads.substring(0, position) + errorBase + originalReads.substring(position+1);
            }
        }

        return mutateRead;
    }
}
