package hypermodulesrun;

import jsc.contingencytables.ContingencyTable2x2;
import jsc.tests.H1;
import jsc.distributions.Hypergeometric;

public class FishersExactTest {

        private int populationSize;
        private int totalSuccesses;
        private int sampleSize;
        private int sampleSuccesses;
        private ContingencyTable2x2 ct;
        private H1 alternative = H1.LESS_THAN;
        private static final int OVERREP = 1;
        private static final int UNDERREP = -1;
        private static final int NOREP = 0;
        private int representation;

        /**
         * A two-tailed Fisher's exact test. Use the getResult() method to get the
         * probability. You will need to have the cern.jet.random package (colt.jar)
         * on your classpath to use this.
         * 
         * @param populationSize
         *            The size of the population.
         * @param totalSuccesses
         *            The total number of "successes" that meet the given criteria.
         * @param sampleSize
         *            The size of the sample.
         * @param sampleSuccesses
         *            The total number of "successes" within the sample that meet
         *            the given criteria.
         */
        public FishersExactTest(int populationSize, int totalSuccesses, int sampleSize, int sampleSuccesses) {
                try {
                        setVariables(populationSize, totalSuccesses, sampleSize, sampleSuccesses);
                } catch (SillyInputException e) {
                        System.out.println(e.toString());
                        e.printStackTrace();
                }
        }

        public FishersExactTest(int populationSize, int totalSuccesses, int sampleSize, int sampleSuccesses, H1 alternative) {
                try {
                        this.alternative = alternative;
                        setVariables(populationSize, totalSuccesses, sampleSize, sampleSuccesses);
                } catch (SillyInputException e) {
                        System.out.println(e.toString());
                        e.printStackTrace();
                }
        }

        /**
         * Reset the variables for the two-tailed Fisher's exact test. Use the
         * getResult() method to get the probability.
         * 
         * @param populationSize
         *            The size of the population
         * @param totalSuccesses
         *            The total number of "successes" that meet the given criteria
         * @param sampleSize
         *            The size of the sample
         * @param sampleSuccesses
         *            The total number of "successes" within the sample that meet
         *            the given criteria
         * @throws SillyInputException
         *             An exception thrown on the input of impossible counts
         */
        public void setVariables(int populationSize, int totalSuccesses, int sampleSize, int sampleSuccesses) throws SillyInputException {
                this.populationSize = populationSize;
                this.totalSuccesses = totalSuccesses;
                this.sampleSize = sampleSize;
                this.sampleSuccesses = sampleSuccesses;
                testInput();
                double expectedProportion = (totalSuccesses + 0.0) / populationSize;
                double expectedSamSuc = sampleSize * expectedProportion;
                if (expectedSamSuc < sampleSuccesses) {
                        representation = OVERREP;
                } else if (expectedSamSuc > sampleSuccesses) {
                        representation = UNDERREP;
                } else {
                        representation = NOREP;
                }

                int[][] data = new int[2][2];
                data[0][0] = sampleSuccesses;
                data[0][1] = sampleSize - sampleSuccesses;
                data[1][0] = totalSuccesses - sampleSuccesses;
                data[1][1] = populationSize - data[0][0] - data[0][1] - data[1][0];
                ct = new ContingencyTable2x2(data);

        }

        public int getRepresentation() {
                return representation;
        }

        /**
         * Get the results of the test
         * 
         * @return The results of the two-tailed Fishers exact test.
         */
        public double getResult() {
        	/*
        		System.out.println(populationSize);
        		System.out.println(totalSuccesses);
        		System.out.println(sampleSize);
        		System.out.println(sampleSuccesses);
        	*/
        	try{
                if (sampleSuccesses != 0 && sampleSuccesses != sampleSize) {
                    jsc.contingencytables.FishersExactTest fet = new jsc.contingencytables.FishersExactTest(ct, alternative);
                    return fet.getSP();
            } else {
                    if (alternative == H1.LESS_THAN) {
                            Hypergeometric hg = new Hypergeometric(sampleSize, populationSize, totalSuccesses);
                            return hg.pdf(sampleSuccesses);
                    } else if(alternative == H1.GREATER_THAN){
                            return 1.0;
                    } else {
                            Hypergeometric hg = new Hypergeometric(sampleSize, populationSize, totalSuccesses);
                            double p0 = hg.pdf(sampleSuccesses);
                            double sum = p0;
                            for(int i=1;i<=sampleSize;i++){
                                    double pI = hg.pdf(i);
                                    if(pI<=p0){
                                            sum+=pI;
                                    }
                            }
                            return sum;
                    }
            }
        		
        	}
        	catch (Exception e){
        		return 1.0;
        	}


        }

        private boolean testNan(double d) {
                if (Double.isNaN(d) || Double.isInfinite(d)) {
                        return true;
                } else {
                        return false;
                }
        }

        private void testInput() throws SillyInputException {
                if (populationSize < totalSuccesses) {
                        throwException();
                } else if (populationSize < sampleSize) {
                        throwException();
                } else if (sampleSize < sampleSuccesses) {
                        throwException();
                } else if (totalSuccesses < sampleSuccesses) {
                        throwException();
                }
        }

        private void throwException() throws SillyInputException {
                throw new SillyInputException(populationSize, totalSuccesses, sampleSize, sampleSuccesses);
        }
}
