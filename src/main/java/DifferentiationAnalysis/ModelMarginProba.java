package DifferentiationAnalysis;


import org.apache.commons.math3.special.Gamma;

import java.util.Arrays;

public class ModelMarginProba {
    private boolean sameLOR;
    private double s1GlobalLOR, s2GlobalLOR, maxLogPosteriorProb, tau, gammaShape, gammaScale;
    private double[] s1ObservedLOR, s2ObservedLOR, s1Variance, s2Variance, s1ExpectedLOR, s2ExpectedLOR;

    public ModelMarginProba(double[] s1ObservedLOR, double[] s2ObservedLOR, double[] s1Variance, double[] s2Variance,
                            double[] s1ExpectedLOR, double[] s2ExpectedLOR, double s1GlobalLOR, double s2GlobalLOR,
                            double maxLogPosteriorProb, double tau, double gammaShape, double gammaScale, boolean sameLOR) {
        this.s1ObservedLOR = s1ObservedLOR;
        this.s1Variance = s1Variance;
        this.s2ObservedLOR = s2ObservedLOR;
        this.s2Variance = s2Variance;
        this.s1ExpectedLOR = s1ExpectedLOR;
        this.s2ExpectedLOR = s2ExpectedLOR;
        this.s1GlobalLOR = s1GlobalLOR;
        this.s2GlobalLOR = s2GlobalLOR;
        this.maxLogPosteriorProb = maxLogPosteriorProb;
        this.tau = tau;
        this.gammaShape = gammaShape;
        this.gammaScale = gammaScale;
        this.sameLOR = sameLOR;
    }

    public double calcMarginProb() {
        double paramNum = this.s1ObservedLOR.length + this.s2ObservedLOR.length + 1;
        paramNum = this.sameLOR? paramNum+1: paramNum+2;
        double cons = paramNum / 2 * Math.log(2*Math.PI);
        double hessianMatrix = 0;
        for(double v: this.s1Variance)
            hessianMatrix += secondDerivativeTheta(v);
        for(double v: this.s2Variance)
            hessianMatrix += secondDerivativeTheta(v);
        hessianMatrix += this.secondDerivativeU();
        hessianMatrix += this.secondDerivativeTau();

        return cons + hessianMatrix + this.maxLogPosteriorProb;
    }

    // ∂^2 f / ∂ theta^2
    private double secondDerivativeTheta(double var) {
        return -1 / var - 1 / Math.pow(this.tau, 2);
    }

    // ∂^2 f / ∂ u^2
    private double secondDerivativeU() {
        double val = 0, s1GlobalLORVar = 0, s2GlobalLORVar = 0;
        for (int idx=0; idx<this.s1ObservedLOR.length; idx++) {
            val -= 1 / Math.pow(this.tau, 4) / (1 / Math.pow(this.tau, 2) + 1 / this.s1Variance[idx]);
            s1GlobalLORVar += 1 / (Math.pow(this.tau, 2) + this.s1Variance[idx]);
        }
        for (int idx=0; idx<this.s2ObservedLOR.length; idx++) {
            val -= 1 / Math.pow(this.tau, 4) / (1 / Math.pow(this.tau, 2) + 1 / this.s2Variance[idx]);
            s2GlobalLORVar += 1 / (Math.pow(this.tau, 2) + this.s2Variance[idx]);
        }
        val -= s1GlobalLORVar;
        val -= s2GlobalLORVar;

        return val;
    }

    // ∂^2 f / ∂ t^2
    private double secondDerivativeTau() {
        double val = 0, lor, var, expLor, prob, totalProba = 0;
        for (int idx=0; idx<this.s1ObservedLOR.length; idx++) {
            lor = this.s1ObservedLOR[idx];
            var = this.s1Variance[idx];
            expLor = this.s1ExpectedLOR[idx];
            prob = this.thetaProbSecDerTau(this.s1GlobalLOR, expLor, lor, var);
            // NaN caused by tiny value
            if (Double.compare(prob, Double.NaN) == 0)
                prob = 0;
            totalProba += prob;
            val += prob;
        }
        val += uProbSecDerTau(this.s1GlobalLOR, this.s1ObservedLOR, this.s1Variance);
        val += secDerTau(this.gammaShape, this.gammaScale, this.s1ObservedLOR, this.s1Variance);
        if (Double.compare(totalProba, 0)!=0) {
            for (int idx=0; idx<this.s2ObservedLOR.length; idx++) {
                lor = this.s2ObservedLOR[idx];
                var = this.s2Variance[idx];
                expLor = this.s2ExpectedLOR[idx];
                prob = this.thetaProbSecDerTau(this.s2GlobalLOR, expLor, lor, var);
                if (Double.compare(prob, Double.NaN) == 0)
                    prob = 0;
                val += prob;
            }
        }
        val += uProbSecDerTau(this.s2GlobalLOR, this.s2ObservedLOR, this.s2Variance);
        val += secDerTau(this.gammaShape, this.gammaScale, this.s2ObservedLOR, this.s2Variance);

        return val;
    }

    private double thetaProbSecDerTau(double u, double expectedLOR, double observedLOR, double variance) {
        double v1 = 1 / Math.pow(this.tau, 2) + 1 / variance;
        double v2 = u / Math.pow(this.tau, 2) + observedLOR / variance;
        double v3 = v2 / v1;
        double v4 = Math.sqrt(1 / v1);
        double v5 = Math.pow(expectedLOR - v3, 2);
        double cons = Math.sqrt(2 * Math.PI);
        double v7 = cons*v4;
        double v8 = 2*u/Math.pow(this.tau, 3)/v1-2*v2/Math.pow(this.tau, 3)/Math.pow(v1, 2);

        return 1/Math.pow(this.tau, 3)*cons*Math.pow(v4, 3)*(-1*v4/cons/Math.pow(this.tau, 3)+(-1*v1*v8*(expectedLOR-v3)+v5/Math.pow(this.tau, 3))/v7) +
               v7*(v1*v8*(expectedLOR-v3)-v5/Math.pow(this.tau, 3))*(v4/cons/Math.pow(this.tau, 3)+(-1*v1*v8*(expectedLOR-v3)+v5/Math.pow(this.tau, 3))/v7) +
               v7*(3*v4/cons/Math.pow(this.tau, 4)-Math.pow(v4, 3)/cons/Math.pow(this.tau, 6)+1/v7*(-1*v1*Math.pow(v8, 2)-v1*(8*u/Math.pow(this.tau, 6)/Math.pow(v1, 2)-6*u/Math.pow(this.tau, 4)/v1-8*v2/Math.pow(this.tau, 6)/Math.pow(v1, 3)+6*v2/Math.pow(this.tau, 4)/Math.pow(v1, 2))*(expectedLOR-v3)+4*v8*(expectedLOR-v3)/Math.pow(this.tau, 3)-3*v5/Math.pow(this.tau, 4))-1/Math.pow(this.tau, 3)*Math.sqrt(2/Math.PI)*v4*(-1*v1*v8*(expectedLOR-v3)+v5/Math.pow(this.tau, 3))+Math.pow(-1*v1*v8*(expectedLOR-v3)+v5/Math.pow(this.tau, 3), 2)/v7);
    }

    private double uProbSecDerTau(double u, double[] observedLOR, double[] variance) {
        double uExp, uVar, uExpDTau, uVarDTau, uExpSecDTau, uVarSecDTau;
        double var, v1 = 0, v2 = 0, v3 = 0;
        double[] kList = new double[observedLOR.length], pList = new double[observedLOR.length],
                 qList = new double[observedLOR.length], sList = new double[observedLOR.length];
        for (int idx=0; idx<observedLOR.length; idx++) {
            var = variance[idx]+Math.pow(this.tau, 2);
            kList[idx] = 1 / var;
            pList[idx] = observedLOR[idx] / var;
            qList[idx] = observedLOR[idx] / Math.pow(var, 2);
            sList[idx] = 1 / Math.pow(var, 2);
            v1 += (2*observedLOR[idx]*Math.pow(var, 2)-2*this.tau*2*var*2*this.tau)/Math.pow(var, 4);
            v2 += (2*observedLOR[idx]*var-2*this.tau*observedLOR[idx]*2*this.tau)/Math.pow(var, 2);
            v3 += (2*Math.pow(var, 2)-2*this.tau*2*var*2*this.tau)/Math.pow(var, 4);
        }
        uExp = Arrays.stream(pList).sum()/ Arrays.stream(kList).sum();
        uVar = 1 / Arrays.stream(kList).sum();
        uExpDTau = -2*this.tau*Arrays.stream(qList).sum()/Arrays.stream(kList).sum() + 2*this.tau* Arrays.stream(pList).sum();
        uVarDTau = 2*this.tau* Arrays.stream(sList).sum()/Math.pow(Arrays.stream(kList).sum(), 2);
        uExpSecDTau = (-1*v1* Arrays.stream(kList).sum()-2*this.tau* Arrays.stream(qList).sum()*2*this.tau* Arrays.stream(sList).sum())/Math.pow(Arrays.stream(kList).sum(), 2)+v2;
        uVarSecDTau = (v3*Math.pow(Arrays.stream(kList).sum(), 2)+2*this.tau* Arrays.stream(sList).sum()*2* Arrays.stream(kList).sum()*2*this.tau* Arrays.stream(qList).sum())/Math.pow(Arrays.stream(kList).sum(), 4);

        double val = (-1/uVar*uExpDTau-(u-uExp)/Math.pow(uVar, 2)*uVarDTau)*uExpDTau+(u-uExp)/uVar*uExpSecDTau;
        double v4 = Math.pow(u-uExp, 2)/2/uVar;
        double b = Math.exp(v4)*Math.sqrt(2*Math.PI)*(Math.exp(-1*v4)*Math.pow(u-uExp, 2)/2/Math.sqrt(2*Math.PI)/Math.pow(uVar, 2.5)-Math.exp(-1*v4)/2/Math.sqrt(2*Math.PI)/Math.pow(uVar, 1.5))*Math.sqrt(uVar);
        double bDUExp = Math.exp(v4)*Math.sqrt(2*Math.PI)*(u-uExp)*(Math.exp(-1*v4)*Math.pow(u-uExp, 2)/2/Math.sqrt(2*Math.PI)/Math.pow(uVar, 2.5)-Math.exp(-1*v4)/2/Math.sqrt(2*Math.PI)/Math.pow(uVar, 1.5))/Math.sqrt(uVar)+
                        Math.exp(v4)*Math.sqrt(2*Math.PI)*(Math.exp(-1*v4)*Math.pow(u-uExp, 3)/2/Math.sqrt(2*Math.PI)/Math.pow(uVar, 3.5)-3*Math.exp(-1*v4)*(u-uExp)/2/Math.sqrt(2*Math.PI)/Math.pow(uVar, 2.5))*Math.sqrt(uVar);
        double bDUVar = Math.exp(v4)*Math.sqrt(Math.PI/2)*(u- uExp)*(Math.exp(-1*v4)*Math.pow(u-uExp, 2)/2/Math.sqrt(2*Math.PI)/Math.pow(uVar, 2.5)-Math.exp(-1*v4)/2/Math.sqrt(2*Math.PI)/Math.pow(uVar, 1.5))/Math.pow(uVar, 1.5) +
                        Math.exp(v4)*Math.sqrt(Math.PI/2)*(Math.exp(-1*v4)*Math.pow(u-uExp, 2)/2/Math.sqrt(2*Math.PI)/Math.pow(uVar, 2.5)-Math.exp(-1*v4)/2/Math.sqrt(2*Math.PI)/Math.pow(uVar, 1.5))/Math.sqrt(uVar) +
                        Math.exp(v4)*Math.sqrt(2*Math.PI)*(Math.exp(-1*v4)*Math.pow(u-uExp, 4)/4/Math.sqrt(2*Math.PI)/Math.pow(uVar, 4.5)-3*Math.exp(-1*v4)*Math.pow(u-uExp, 2)/2/Math.sqrt(2*Math.PI)/Math.pow(uVar, 3.5)+3*Math.exp(-1*v4)/4/Math.sqrt(2*Math.PI)/Math.pow(uVar, 2.5))*Math.sqrt(uVar);

        val += (bDUExp*uExpDTau+bDUVar*uVarDTau)*uVarDTau+b*uVarSecDTau;

        return val;
    }

    private double secDerTau(double gammaShape, double gammaScale, double[] observedLOR, double[] variance) {
        double uExp, uVar, uExpDTau, uVarDTau, uExpSecDTau, uVarSecDTau;
        double var, v1 = 0, v2 = 0, v3 = 0;
        double[] kList = new double[observedLOR.length], pList = new double[observedLOR.length],
                qList = new double[observedLOR.length], sList = new double[observedLOR.length];
        for (int idx=0; idx<observedLOR.length; idx++) {
            var = variance[idx]+Math.pow(this.tau, 2);
            kList[idx] = 1 / var;
            pList[idx] = observedLOR[idx] / var;
            qList[idx] = observedLOR[idx] / Math.pow(var, 2);
            sList[idx] = 1 / Math.pow(var, 2);
            v1 += (2*observedLOR[idx]*Math.pow(var, 2)-2*this.tau*2*var*2*this.tau)/Math.pow(var, 4);
            v2 += (2*observedLOR[idx]*var-2*this.tau*observedLOR[idx]*2*this.tau)/Math.pow(var, 2);
            v3 += (2*Math.pow(var, 2)-2*this.tau*2*var*2*this.tau)/Math.pow(var, 4);
        }
        uExp = Arrays.stream(pList).sum()/ Arrays.stream(kList).sum();
        uVar = 1 / Arrays.stream(kList).sum();
        uExpDTau = -2*this.tau*Arrays.stream(qList).sum()/Arrays.stream(kList).sum() + 2*this.tau* Arrays.stream(pList).sum();
        uVarDTau = 2*this.tau* Arrays.stream(sList).sum()/Math.pow(Arrays.stream(kList).sum(), 2);
        uExpSecDTau = (-1*v1* Arrays.stream(kList).sum()-2*this.tau* Arrays.stream(qList).sum()*2*this.tau* Arrays.stream(sList).sum())/Math.pow(Arrays.stream(kList).sum(), 2)+v2;
        uVarSecDTau = (v3*Math.pow(Arrays.stream(kList).sum(), 2)+2*this.tau* Arrays.stream(sList).sum()*2* Arrays.stream(kList).sum()*2*this.tau* Arrays.stream(qList).sum())/Math.pow(Arrays.stream(kList).sum(), 4);

        double val = 0.5*(-1*uVarDTau*uVarDTau/Math.pow(uVar, 2) + 1/uVar*uVarSecDTau);
        double v4 = 0, v5 = 0, v6 = 0, v7 = 0, v8 = 0;
        for (int idx=0; idx<observedLOR.length; idx++) {
            var = variance[idx]+Math.pow(this.tau, 2);
            v4 += (2*var-4*Math.pow(this.tau, 2))/Math.pow(var, 2);
            v5 += (-1*Math.pow(uExpDTau, 2)+(observedLOR[idx]-uExp)*uExpSecDTau)/2/Math.pow(var, 2);
            v6 += (-8*(observedLOR[idx]-uExp)*uExpDTau*this.tau+Math.pow(observedLOR[idx]-uExp, 2)*4)/4/Math.pow(var, 2);
            v7 += 2*(observedLOR[idx]-uExp)*uExpDTau*this.tau/Math.pow(var, 3);
            v8 += 4*Math.pow(observedLOR[idx]-uExp, 2)*Math.pow(this.tau, 2)/Math.pow(var, 3);
        }

        val -= 0.5*(v4 + v5 - v6 - v7 - v8);

        val += (1-gammaShape)*Math.log(this.tau) - 1/this.tau/gammaScale - Gamma.logGamma(gammaShape) - gammaShape*Math.log(gammaScale);

        return val;
    }
}
