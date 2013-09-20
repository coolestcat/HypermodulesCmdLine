package hypermodulesrun;


import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
/**
 * 
 * False Discovery Rate P-Value adjustment class - basically, if we get a p-value for a given module, 
 * we then randomize the associations of genes with patients (without changing the topology of the network);
 * if many p-values that we get from the random networks are lower than the original p-value, this means our original p-value 
 * is not "special"; in other words, the original association of genes to patients doesn't carry much weight, as all networks with 
 * the given topology would have a low p-value. However, if the FDR p-value is low, then not many of the randomly shuffled 
 * networks give a low p-value, so our original result is significant.
 * @author alvinleung
 *
 */
public class FDRAdjust {
	
	private HashMap<String, Double> this_true;
	private Multimap<String, Double> this_rand;
	private static final int RESULT_SCALE = 100;
	
	
	public FDRAdjust(HashMap<String, Double> this_true, Multimap<String, Double> this_rand){
		this.this_true = this_true;
		this.this_rand = this_rand;
	}
	
	public HashMap<String, Double> fdrAdjust(){
		
		if (this_rand.isEmpty()){
			return this_true;
		}
		
		HashMap<String, Double> adjustedResults = new HashMap<String, Double>();
		
		Multimap<Double, String> reversedResults = ArrayListMultimap.create();
		//Multimap<Double, String> reversedRandomResults = ArrayListMultimap.create();
		
		
		
		for (String reverse : this_true.keySet()){
			reversedResults.put(this_true.get(reverse), reverse);
		}
		
		
		/*
		for (String reverse : this_rand.keySet()){
			for (Double doubleAdd : this_rand.get(reverse)){
				reversedRandomResults.put(doubleAdd, reverse);
			}
		}
		 */
		
		ArrayList<Double> resultDoubles = new ArrayList<Double>();
		
		
		int k=0;
		for (String result : this_true.keySet()){
			resultDoubles.add(Double.valueOf(this_true.get(result)));
			k++;
		}
		
		Collections.sort(resultDoubles);
		
		double[] randomResultDoubles = new double[this_rand.size()];
		k=0;
		for(String result : this_rand.keySet()){
			for (Double getDouble : this_rand.get(result)){
				randomResultDoubles[k]=getDouble;
				k++;
			}
		}

		double[] toAdjust = new double[resultDoubles.size()];
		
		for (int i=0; i<toAdjust.length; i++){
			toAdjust[i]=lengthWhich(resultDoubles.get(i), randomResultDoubles)/(double) randomResultDoubles.length;
		}
		
		/*
		for (int i=0; i<toAdjust.length; i++){
			System.out.println(toAdjust[i]);
		}
		*/
		
		/*
		System.out.println(randomResultDoubles.length);
		System.out.println(lengthWhich(resultDoubles.get(0), randomResultDoubles));
		
		System.out.println("toAdjust for Java: ");
		for (int i=0; i<toAdjust.length; i++){
			System.out.println(toAdjust[i]);
		}
		*/
		
		//perform FDR adjustment
		double[] adjustedDoubles = new double[resultDoubles.size()];
        BigDecimal min = new BigDecimal("" + 1);
        BigDecimal mkprk;
        int m = toAdjust.length;
        for (int i = m; i > 0; i--) {
            mkprk = (new BigDecimal("" + m).multiply(new BigDecimal(toAdjust[i - 1]))).divide(new BigDecimal("" + i), RESULT_SCALE, BigDecimal.ROUND_HALF_UP);
            if (mkprk.compareTo(min) < 0) {
                min = mkprk;
            }
            adjustedDoubles[i - 1] = Double.parseDouble(min.toString());
        }

		for (int i=0; i<adjustedDoubles.length; i++){
			Double d = adjustedDoubles[i];
			//d = (double)Math.round(d * 10000) / 10000;
			for (String s : reversedResults.get(resultDoubles.get(i))){
				adjustedResults.put(s, d);
			}
		}
		
		/*
		for (int i=0; i<adjustedDoubles.length; i++){
			System.out.println(adjustedDoubles[i]);
		}
		*/
		
		System.out.println(adjustedResults.size());
		return adjustedResults;
		
	}
	
	public double lengthWhich(double real, double[] random){
		double count = 0;
		for (int i=0; i<random.length; i++){
			if (random[i]<=real){
				count++;
			}
		}
		
		return count;
		
	}
	

}



