package hypermodulesrun;

public class SillyInputException extends Exception{
    private final int populationSize;
    private final int totalSuccesses;
    private final int sampleSize;
    private final int sampleSuccesses;
    
    public SillyInputException(int pSize, int pSuc, int sSize, int sSuc){
            this.populationSize = pSize;
            this.totalSuccesses = pSuc;
            this.sampleSize = sSize;
            this.sampleSuccesses = sSuc;
            
    }
    @Override
    public String toString() {
            return "\nThe population size: "+populationSize+" must be >= the total successes: "+totalSuccesses+".\n\t" +
                            "The population size: "+populationSize+" must be >= the sample size: "+sampleSize+".\n" +
                            "The sample size: "+sampleSize+" must be >= the sample successes: "+sampleSuccesses+".\n" +
                            "The total successes: "+totalSuccesses+" must be >= sample successes: "+sampleSuccesses+".\n\t"+
                                            "Or this exception will be thrown.\n";
    
    }
}

