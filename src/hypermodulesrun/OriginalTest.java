package hypermodulesrun;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyTableUtil;

public class OriginalTest{
	
	private ArrayList<String[]> network;
	private ArrayList<String[]> sampleValues;
	private ArrayList<String[]> clinicalValues;
	
	private int shuffleNumber;
	private String stat;
	
	
	public OriginalTest(ArrayList<String[]> network, 
			 			ArrayList<String[]> sampleValues, 
			 			ArrayList<String[]> clinicalValues,
			 			int shuffleNumber,
			 			String stat){
		
		this.network = network;
		this.sampleValues = sampleValues;
		this.shuffleNumber = shuffleNumber;
		this.clinicalValues = clinicalValues;
		this.stat = stat;
	}
	
	public HashMap<String, HashMap<String, Double>> callTest(){
		HashMap<String, HashMap<String, Double>> rt = new HashMap<String, HashMap<String, Double>>();
		HypermodulesHeuristicAlgorithm ha = new HypermodulesHeuristicAlgorithm(this.stat, this.sampleValues, this.clinicalValues, this.network);
		ha.initialize();

		HashSet<String> allSeeds = new HashSet<String>();
		for (int i=0; i<sampleValues.size(); i++){
			if (!sampleValues.get(i)[1].equals("no_sample") && sampleValues.get(i)[1]!=null){
					allSeeds.add(sampleValues.get(i)[0]);
			}
		}
		
		System.out.println("allSeeds size: " + allSeeds.size());
			
		for (String runSeed : allSeeds){
			HashMap<String, Double> oneResult = testSeed(ha, runSeed);
			rt.put(runSeed, oneResult);		
		}
		System.out.println("Finished running Original Test.");
		
		return rt;
	}
	
	public HashMap<String, Double> testSeed(HypermodulesHeuristicAlgorithm ha, String seedName){
		HashMap<String, Double> returnMap = new HashMap<String, Double>();
		
		HashSet<String> allPaths = new HashSet<String>();
		allPaths = ha.getAllPaths(seedName);
		
		System.out.println("ALL PATHS SIZE: " + allPaths.size());
		
		
		ArrayList<String> compress = ha.compressTokens(allPaths, seedName);
		
		System.out.println("COMPRESSED SIZE: " + compress.size());
		
		HashMap<String, Double> answer = ha.mineHublets(compress);
		returnMap = answer;
		
		System.out.println("FINAL SIZE: " + returnMap.size());
		return returnMap;
	}
	
	public HashMap<String, HashMap<String, Double>> testHighOrLow(HashMap<String, HashMap<String, Double>> ot){
		HashMap<String, HashMap<String, Double>> rt = new HashMap<String, HashMap<String, Double>>();
		HypermodulesHeuristicAlgorithm ha = new HypermodulesHeuristicAlgorithm(this.stat, this.sampleValues, this.clinicalValues, this.network);
		ha.initialize();
		for (String s : ot.keySet()){
			HashMap<String, Double> newMap = new HashMap<String, Double>();
			for (String t : ot.get(s).keySet()){
				//System.out.println(t);
				if (ha.testModuleBoolean(t)==1){
					newMap.put(t, 1.0);
				}
				else if (ha.testModuleBoolean(t)==0){
					newMap.put(t, 0.0);
				}
				else{
					newMap.put(t, 2.0);
				}
			}
			rt.put(s, newMap);
		}
		
		return rt;
	}
	
}
