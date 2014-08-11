package hypermodulesrun;

import jsc.contingencytables.ContingencyTable2x2;
import jsc.contingencytables.FishersExactTest;
import jsc.tests.H1;
import jsc.distributions.Hypergeometric;

public class MyFET{

        private int populationSize;
        private int totalSuccesses;
        private int sampleSize;
        private int sampleSuccesses;
        private ContingencyTable2x2 ct;
        private H1 alternative = H1.GREATER_THAN;
        private static final int OVERREP = 1;
        private static final int UNDERREP = -1;
        private static final int NOREP = 0;
        private int representation;
        private int[][] data;

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
        public MyFET(int populationSize, int totalSuccesses, int sampleSize, int sampleSuccesses) {
                try {
            		//TODO: erase this.
                	/*
                	this.populationSize = 40;
                	this.totalSuccesses = 22;
                	this.sampleSize = 27;
                	this.sampleSuccesses = 16;
                	*/
                        setVariables(populationSize,totalSuccesses,sampleSize,sampleSuccesses);
                } catch (SillyInputException e) {
                        System.out.println(e.toString());
                        e.printStackTrace();
                }
        }

        public MyFET(int populationSize, int totalSuccesses, int sampleSize, int sampleSuccesses, H1 alternative) {
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

                data = new int[2][2];
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
         * @return The results of the one-tailed Fishers exact test//two-tailed Fishers exact test.
         */
        public double getResult() {
        	try{
        		double pval = 0;
        		if (ct.getFrequency(0,0)!=0 && ct.getFrequency(1,1)!=0 &&
        				ct.getFrequency(1,0)!=0 && ct.getFrequency(0,1)!=0) {
        			FishersExactTest fet = new jsc.contingencytables.FishersExactTest(ct, H1.GREATER_THAN);
        			pval = fet.getSP();
        		} else {
        			Hypergeometric hg = new Hypergeometric(sampleSize, populationSize, totalSuccesses);
        			for (int ss=sampleSuccesses; ss<=Math.min(totalSuccesses, sampleSize); ss++) {
        				pval += hg.pdf(ss);
        			}
            	}
        		
        		//System.out.println(data[0][0] + " : " + data[0][1] + " : " + data[1][0] + " : " + data[1][1] + " : " + pval);
        		return pval;

        	}
        	catch (Exception e){
        		return 1.0;
        	}
        	finally{
        		
        	}


        }
        
        
        public Double getLogOdds(){
			double lodds1 = (ct.getFrequency(1,1)*ct.getFrequency(0,0));
			double lodds2 = (ct.getFrequency(0,1)*ct.getFrequency(1,0));
			
			Double lodds = lodds1/lodds2;
			
			if (!lodds.isNaN() && !lodds.isInfinite()){
				return Math.log(lodds);
			}
			if (lodds == Double.POSITIVE_INFINITY){
				return 1000.0;
			}
			if (lodds == 0){
				return -1000.0;
			}
			
			return lodds;
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
