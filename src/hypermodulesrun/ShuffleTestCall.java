package hypermodulesrun;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.Callable;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyTableUtil;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

public class ShuffleTestCall implements Callable<HashMap<String, Multimap<String, Double>>> {
	
	private ArrayList<String[]> network;
	private ArrayList<String[]> sampleValues;
	private ArrayList<String[]> clinicalValues;
	
	private int shuffleNumber;
	private String stat;
	private String foregroundvariable;
	
	
	public ShuffleTestCall(ArrayList<String[]> network, 
			 			ArrayList<String[]> sampleValues, 
			 			ArrayList<String[]> clinicalValues,
			 			int shuffleNumber,
			 			String stat,
			 			String foregroundvariable){
		
		this.network = network;
		this.sampleValues = sampleValues;
		this.shuffleNumber = shuffleNumber;
		this.clinicalValues = clinicalValues;
		this.stat = stat;
		this.foregroundvariable = foregroundvariable;
	}


	public HashMap<String, Multimap<String, Double>> call() throws Exception {
		HashMap<String, Multimap<String, Double>> rt = new HashMap<String, Multimap<String, Double>>();
		HypermodulesHeuristicAlgorithm ha = new HypermodulesHeuristicAlgorithm(this.stat, this.foregroundvariable, this.sampleValues, this.clinicalValues, this.network);
		ha.initialize();

		
			HashSet<String> allSeeds = new HashSet<String>();
			for (int i=0; i<sampleValues.size(); i++){
				if (!sampleValues.get(i)[1].equals("no_sample") && sampleValues.get(i)[1]!=null){
					allSeeds.add(sampleValues.get(i)[0]);
				}
			}
			
			HashSet<String> filteredSeeds = new HashSet<String>();
			
			for (int i=0; i<network.size(); i++){
				if (allSeeds.contains(network.get(i)[0])){
					filteredSeeds.add(network.get(i)[0]);
				}
				if (allSeeds.contains(network.get(i)[1])){
					filteredSeeds.add(network.get(i)[1]);
				}
			}
			
			for (String runSeed : filteredSeeds){
				Multimap<String, Double> oneResult = testSeed(ha, runSeed);
				rt.put(runSeed, oneResult);
			}
			
			//System.out.println("finished running.");
			
		return rt;
	}
	
	public Multimap<String, Double> testSeed(HypermodulesHeuristicAlgorithm ha, String seedName){
		
		HashSet<String> allPaths = new HashSet<String>();
		allPaths = ha.getAllPaths(seedName);

		Multimap<String, Double> returnMap = ArrayListMultimap.create();
		
    	for (int x = 0; x < this.shuffleNumber; x++){
        	ha.shuffleLabels();
        	ArrayList<String> compress = ha.compressTokens(allPaths, seedName);
        	HashMap<String, Double> shuffledAnswer = ha.mineHublets(compress);
        	for (String s : shuffledAnswer.keySet()){
        		returnMap.put(s, shuffledAnswer.get(s));
        	}
    	}

		return returnMap;
	}
	
	
}
