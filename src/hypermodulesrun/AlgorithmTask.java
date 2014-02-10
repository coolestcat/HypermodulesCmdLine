package hypermodulesrun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import javax.swing.JOptionPane;

import org.cytoscape.util.swing.FileUtil;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

public class AlgorithmTask implements Runnable {

	/**
	 * original pValues from OriginalTest
	 */
	private HashMap<String, HashMap<String, Double>> originalResults;
	/**
	 * classification of "HIGH" or "LOW" survival for all the modules in originalResults
	 */
	private HashMap<String, HashMap<String, Double>> classification;
	/**
	 * all the shuffled data to perform FDR adjustment
	 */
	private HashMap<String, Multimap<String, Double>> shuffling;
	/**
	 * an arraylist to add to in order to manage having multiple threads - to avoid synchronization
	 * issues related to adding to a hashmap simultaneously from different threads
	 */
	private ArrayList<HashMap<String, Multimap<String, Double>>> combinedShuffling;
	/**
	 * FDR p-values
	 */
	private HashMap<String, HashMap<String, Double>> adjustedResults;
	
	private HashMap<String, HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>>> allResults;
	
	private ArrayList<String[]> network;
	private ArrayList<String[]> sampleValues;
	private ArrayList<String[]> filteredSampleValues;
	private ArrayList<String[]> clinicalValues;
	
	private int shuffleNumber;

	private String stat;
	
	private double pValueCutoff;
	
	private int numberCores;
	
	public AlgorithmTask(ArrayList<String[]> network, 
						 ArrayList<String[]> sampleValues, 
						 ArrayList<String[]> clinicalValues,
						 int shuffleNumber, 
						 String stat,
						 double pValueCutoff,
						 int numberCores){
		
		this.network = network;
		this.sampleValues = sampleValues;
		this.shuffleNumber = shuffleNumber;
		this.stat = stat;
		this.pValueCutoff = pValueCutoff;
		this.numberCores = numberCores;
		
		if (stat.toUpperCase().equals("LOGRANK")){
			
		ArrayList<String[]> sortedClinicals = new ArrayList<String[]>();
		
		Multimap<Double, Integer> followupDaysPosition = ArrayListMultimap.create();
		
		for (int k=0; k<clinicalValues.size(); k++){
			boolean b = true;
			try{
				Double d = Double.valueOf(clinicalValues.get(k)[2]);
			}
			catch (NumberFormatException e){
				b = false;
			}
			finally{
				if (b){
					followupDaysPosition.put(Double.valueOf(clinicalValues.get(k)[2]), Integer.valueOf(k));
				}
			}
		}
		
		ArrayList<Double> sortedTime = new ArrayList<Double>();
		ArrayList<Integer> sortedPositions = new ArrayList<Integer>();
		
		for (int k=0; k<clinicalValues.size(); k++){
			sortedTime.add(Double.valueOf(clinicalValues.get(k)[2]));
		}
		
		Collections.sort(sortedTime);
		HashMap<Double, Boolean> alreadyAdded = new HashMap<Double, Boolean>();
		
		for (Double key : followupDaysPosition.keySet()){
			alreadyAdded.put(key,  false);
		}

		for (int k=0; k<sortedTime.size(); k++){
			if (alreadyAdded.get(sortedTime.get(k))==false){
				Collection<Integer> coll = followupDaysPosition.get(sortedTime.get(k));
				for (Integer value : coll){
					sortedPositions.add(value);
				}
			}
			alreadyAdded.put(sortedTime.get(k), true);
		}
		
		for (int k=0; k<sortedPositions.size(); k++){
			String[] thisString = new String[3];
			thisString[0]=clinicalValues.get(sortedPositions.get(k))[0];
			thisString[1]=clinicalValues.get(sortedPositions.get(k))[1];
			thisString[2]=clinicalValues.get(sortedPositions.get(k))[2];
			//thisString[3]=clinicalValues.get(sortedPositions.get(k))[3];
			sortedClinicals.add(thisString);
			
		}
		this.clinicalValues = sortedClinicals;
		}
		else{
			this.clinicalValues = clinicalValues;
		}
	}
	
	public void run() {
		long before = System.nanoTime();
		System.err.println("Running Original Test...");
		OriginalTest ot = new OriginalTest(network, sampleValues, clinicalValues, shuffleNumber, stat);
		this.originalResults = ot.callTest();
		fixOriginalResults();
		
		if (stat.toUpperCase().equals("LOGRANK")){
			this.classification = ot.testHighOrLow(this.originalResults);
		}

		int nCores = numberCores;
		//int nCores = Runtime.getRuntime().availableProcessors();
		System.err.println("Number of Available Processors: " + nCores);
		System.err.println("Computing Background from " + shuffleNumber + " Random Networks");
		combinedShuffling = new ArrayList<HashMap<String, Multimap<String, Double>>>();
		
		int shuffleCount = 0;
		
		ExecutorService executor = Executors.newFixedThreadPool(nCores);
		List<Future<HashMap<String, Multimap<String, Double>>>> list = new ArrayList<Future<HashMap<String, Multimap<String, Double>>>>();
		
		for (int i=0; i<nCores-1; i++){
			//reinitializeVariables();
			ShuffleTestCall st = new ShuffleTestCall(this.network, this.sampleValues, this.clinicalValues, (int) shuffleNumber/nCores, this.stat);
			Future<HashMap<String, Multimap<String, Double>>> submitPool = executor.submit(st);
			list.add(submitPool);
			shuffleCount += (int) shuffleNumber/nCores;
		}
		
		ProgressBarShuffle st = new ProgressBarShuffle(this.network, this.sampleValues, this.clinicalValues, (int) (shuffleNumber - shuffleCount), this.stat);
		Future<HashMap<String, Multimap<String, Double>>> submitPool = executor.submit(st);
		list.add(submitPool);
		
		for (Future<HashMap<String, Multimap<String, Double>>> future : list){
			try{
				combinedShuffling.add(future.get());
			}
			catch (Exception e){
				e.printStackTrace();
			}
		}
		
		
		executor.shutdown();
		System.err.println();
		System.err.println("Finished Running on All Cores");
		
		moveShuffled();
		System.err.println("Finished Moving");
		
		adjustResults();
		System.err.println("Finished Adjusting");
		
		allResults = resultsFormat();
		HashMap<String, String> parameters = new HashMap<String, String>();
		parameters.put("nShuffled", String.valueOf(this.shuffleNumber));
		parameters.put("stat", this.stat);

		long after = System.nanoTime();
		double timeToRun = (double) (after-before)/1000000000;
		System.err.println("Time to run: " + timeToRun + " seconds");
		
		justPrint();
		//loop();
	}

	
	public void loop(){
		int counter = 0;
		printOptions();
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		while (true){
			int i = 5;
			try {
				i = Integer.valueOf(br.readLine());
				if (i<0 || i>2){
					System.out.println("Please enter a number between 0 and 2");
					printOptions();
				}
			} catch (NumberFormatException e) {
				System.out.println("Please enter a number between 0 and 2");
			} catch (IOException e) {
				System.out.println("Please enter a number between 0 and 2");
			}
			
			if (i==0){
				exportMostCorrelated();
				printOptions();
			}
			else if (i==1){
				printToScreen();
				printOptions();
			}
			else if (i==2){
				System.exit(0);
			}
		}
	}
	
	public void printOptions(){
		System.out.println("Options: ");
		System.out.println("0 - Export results to a CSV file");
		System.out.println("1 - Print results to screen");
		System.out.println("2 - Exit Program (WARNING - ALL RESULTS WILL BE ERASED)");
	}
	
	public void justPrint(){
		ArrayList<HashMap<String, Double>> a = extractMostCorrelated(pValueCutoff);

		HashMap<String, Double> mostCorrelated = a.get(0);
		HashMap<String, Double> mostCorrelatedFDR = a.get(1);
		HashMap<String, Double> mostCorrelatedPatn = a.get(2);
		HashMap<String, Double> mostCorrelatedOddsRatio = a.get(3);
		
		System.out.println("# HyperModules Results");
		System.out.println("# Date: " + '\t' + DateFormat.getDateTimeInstance().format(new Date()));
		
		System.out.println("# Shuffle Number: "+ '\t' + shuffleNumber );
		System.out.println("# Statistical Test: " + '\t'+ stat );
		System.out.println();
		
		System.out.println("Module" + '\t' + "Pvalue_test" + '\t' + "Pvalue_background"  + '\t' + "Number_patients" + '\t' + "Log_odds_ratio" );
		for (String s : mostCorrelated.keySet()){
			System.out.println(s + '\t' + mostCorrelated.get(s) + '\t' + mostCorrelatedFDR.get(s) + '\t' + mostCorrelatedPatn.get(s) + '\t' + mostCorrelatedOddsRatio.get(s));
		}
	}
	
	public void printToScreen(){
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		System.out.println("Please enter a p-Value cutoff (between 0 and 1)");
		double i = 0.05;
		try {
			i = Double.valueOf(br.readLine());
			if (i<0 || i>1){
				System.out.println("Please enter a p-Value cutoff between 0 and 1");
			}
		} catch (NumberFormatException e) {
			System.out.println("Please enter a p-Value cutoff between 0 and 1");
			return;
		} catch (IOException e) {
			System.out.println("Please enter a p-Value cutoff between 0 and 1");
			return;
		}
		
		ArrayList<HashMap<String, Double>> a = extractMostCorrelated(i);

		HashMap<String, Double> mostCorrelated = a.get(0);
		HashMap<String, Double> mostCorrelatedFDR = a.get(1);
		

		
		System.out.println("Module" + '-' + "Statistical Test P-Value" + '-' + "FDR P-Value");
		for (String s : mostCorrelated.keySet()){
			System.out.println(s + '-' + mostCorrelated.get(s) + '-' + mostCorrelatedFDR.get(s));
		}
		
	}
	
	public ArrayList<HashMap<String, Double>> extractMostCorrelated(double cutoff){
		ArrayList<HashMap<String, Double>> rt = new ArrayList<HashMap<String, Double>>();
		HashMap<String, Double> mostCorrelated = new HashMap<String, Double>();
		HashMap<String, Double> mostCorrelatedFDR = new HashMap<String, Double>();
		HashMap<String, Double> mostCorrelatedPatn = new HashMap<String, Double>();
		HashMap<String, Double> mostCorrelatedOddsRatio = new HashMap<String, Double>();
		
		for (String s : allResults.keySet()){
			for (ArrayList<HashMap<String, Double>> ahhs : allResults.get(s).keySet()){
				HashMap<String, Double> original = ahhs.get(0);
				HashMap<String, Double> adjusted = ahhs.get(1);
				HashMap<String, Double> patn = null;
				HashMap<String, Double> oddsratio = null;
				if (this.stat.equals("logRank")){
					patn = ahhs.get(3);
					oddsratio = ahhs.get(4);
				}
				else{
					patn = ahhs.get(2);
					oddsratio = ahhs.get(3);
					
				}
				for (String set : original.keySet()){
					if (adjusted.containsKey(set)){
						if (original.get(set)<cutoff && adjusted.get(set)<cutoff){
							mostCorrelated.put(set, original.get(set));
							mostCorrelatedFDR.put(set, adjusted.get(set));
							mostCorrelatedPatn.put(set, patn.get(set));
							mostCorrelatedOddsRatio.put(set, oddsratio.get(set));
						}
					}
				}
			}
		}
		
		rt.add(mostCorrelated);
		rt.add(mostCorrelatedFDR);
		rt.add(mostCorrelatedPatn);
		rt.add(mostCorrelatedOddsRatio);
		
		return rt;
		
	}
	
	private void fixOriginalResults(){
		HashMap<String, HashMap<String, Double>> newOriginal = new HashMap<String, HashMap<String, Double>>();
		for (String s : this.originalResults.keySet()){
			HashMap<String, Double> newHash = new HashMap<String, Double>();
			for (String t : this.originalResults.get(s).keySet()){
				
				String[] st = t.split(":");
				HashSet<String> hs = new HashSet<String>();
				for (int i=0; i<st.length; i++){
					hs.add(st[i]);
				}
				String newString = "";
				for (String u : hs){
					newString = newString + u + ":";
				}
				newString = newString.substring(0, newString.length()-1);
				Double d = this.originalResults.get(s).get(t);
				//d = (double)Math.round(d * 10000) / 10000;
				newHash.put(newString, d);
			}
			newOriginal.put(s,  newHash);
		}
		
		this.originalResults = newOriginal;
		
	}
	
	private void moveShuffled(){
		shuffling = new HashMap<String, Multimap<String, Double>>();
		HashMap<String, Multimap<String, Double>> lastCore = combinedShuffling.get(combinedShuffling.size()-1);
		for (String str : lastCore.keySet()){
			Multimap<String, Double> mhsd = ArrayListMultimap.create();
			for (String hs : lastCore.get(str).keySet()){
				for (Double d : lastCore.get(str).get(hs)){
					mhsd.put(hs, d);
				}
			}
			shuffling.put(str, mhsd);
		}
		
		for (int i=1; i<combinedShuffling.size(); i++){
			HashMap<String, Multimap<String, Double>> thisCore = combinedShuffling.get(i);
			for (String s : thisCore.keySet()){
				for (String hs : thisCore.get(s).keySet()){
					for (Double d : thisCore.get(s).get(hs)){
						shuffling.get(s).put(hs, d);
					}
				}
			}
		}
	}
	
	private void adjustResults(){
		adjustedResults = new HashMap<String, HashMap<String, Double>>();
		for (String s : originalResults.keySet()){
			FDRAdjust fdr = new FDRAdjust(originalResults.get(s), shuffling.get(s));
			HashMap<String, Double> adjusted = fdr.fdrAdjust();
			adjustedResults.put(s, adjusted);
		}
	}
	
	private HashMap<String, HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>>> resultsFormat(){
		HashMap<String, HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>>> allResults = new HashMap<String, HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>>>();
		for (String s : adjustedResults.keySet()){
			ArrayList<HashMap<String, Double>> ah = new ArrayList<HashMap<String, Double>> ();
			ah.add(originalResults.get(s));
			ah.add(adjustedResults.get(s));
			if (stat.toUpperCase().equals("LOGRANK")){
				ah.add(classification.get(s));
			}
			//ah.add(adjustedWithR.get(s));
			HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>> hah = new HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>>();
			hah.put(ah,  null);
			allResults.put(s, hah);
		}

		
		//filter redundant values in clinicalValues and otherValues?
		
		//filter redundant values in sampleValues
		HashSet<String> sam = new HashSet<String>();
		for (int i=0; i<sampleValues.size(); i++){
			String add = sampleValues.get(i)[0] + ";" + sampleValues.get(i)[1];
			//System.out.println(add);
			sam.add(add);
		}
		
		this.filteredSampleValues = new ArrayList<String[]>();
		for (String s : sam){
			String[] t = s.split(";");
			//System.out.println(t[0] + " : " + t[1]);
			filteredSampleValues.add(t);
		}
		
		for (int i=0; i<filteredSampleValues.size(); i++){
			//System.out.println(filteredSampleValues.get(i)[0] + " : " + filteredSampleValues.get(i)[1]);
		}
		//System.out.println("finished filtering sampleValues");
		allResults = filterRedundantResults(allResults);
		return allResults;

	}
	
	

	/**
	 * Filters the "redundant" results (same p-value, same patients, or same genes + more genes)
	 * Adds a column for number of patients in the module
	 * Adds a column for the odds ratio
	 * Makes sure the genes string has the seed at the head (seed:b:c:...)
	 * 
	 * ArrayList<HashMap<String, Double>>: 
	 * .get(0) - original
	 * .get(1) - adjusted
	 * .get(2) - classification
	 * .get(3) (or 2) - numberPatients
	 * .get(4) (or 3) - odds ratio
	 * @param input
	 * @return
	 */
	public HashMap<String, HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>>> filterRedundantResults(HashMap<String, HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>>> input){
		HashMap<String, HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>>> output = new HashMap<String, HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>>>();
		
		HashMap<String, Double> masterList = new HashMap<String, Double>();
		
		for (String s : input.keySet()){
			for (ArrayList<HashMap<String, Double>> ahsd : input.get(s).keySet()){
				HashMap<String, Double> original = ahsd.get(0);
				for (String i : original.keySet()){
					masterList.put(i, original.get(i));
					//System.out.println(i + " : " + original.get(i));
				}
			}
		}
		
		HashMap<String, Double> m2 = new HashMap<String, Double>();
		for (String s : masterList.keySet()){
			m2.put(s, masterList.get(s));
		}
		
		HashSet<String> rejectedList = new HashSet<String>();
		
		int y = 0;
		System.out.println("masterlist size: " + masterList.size());
		//filter the masterList (O (n^2)):
		for (String s : masterList.keySet()){
			double d = masterList.get(s);
			for (String t : m2.keySet()){
				double e = m2.get(t);
				if ((d == e) && !(s.equals(t))){
					//System.out.println(d + " : " + e + " : " + s + " : " + t);
					if (checkConditions(s, t)){//t is redundant to s
						rejectedList.add(t);
					}
				}	
			}
			y++;
			//System.out.println(y);
		}
		
		System.out.println("rejectedlist size:"  + rejectedList.size());
		System.out.println("finished creating rejection list");
		
		for (String s : input.keySet()){
			//output.put(s, input.get(s));
			
			HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>> hah = new HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>>();
			
			HashMap<ArrayList<HashMap<String, Double>>, Multimap<String, Double>> inputhah = input.get(s);
			
			for (ArrayList<HashMap<String, Double>> ahsd : inputhah.keySet()){
				ArrayList<HashMap<String, Double>> newahsd = new ArrayList<HashMap<String, Double>>();
				
				HashMap<String, Double> orig = ahsd.get(0);
				HashMap<String, Double> neworig = new HashMap<String, Double>();
				for (String o : orig.keySet()){
					if (!rejectedList.contains(o)){
						neworig.put(seedAtBeginning(s, o), roundToSignificantFigures(orig.get(o),5));
					}
				}
				
				HashMap<String, Double> adj = ahsd.get(1);
				HashMap<String, Double> newadj = new HashMap<String, Double>();
				for (String o : adj.keySet()){
					if (!rejectedList.contains(o)){
						newadj.put(seedAtBeginning(s, o), roundToSignificantFigures(adj.get(o),5));
					}
				}
				
				if (this.stat.equals("logRank")){
					HashMap<String, Double> clas = ahsd.get(2);
					HashMap<String, Double> newclas = new HashMap<String, Double>();
					for (String o : clas.keySet()){
						if (!rejectedList.contains(o)){
							newclas.put(seedAtBeginning(s, o), roundToSignificantFigures(clas.get(o),5));
						}
					}
					
					HashMap<String, Double> patn = new HashMap<String, Double>();
					for (String x : neworig.keySet()){
						patn.put(x, (double) getNumPatients(x));
					}
					
					HashMap<String, Double> oddsratio = new HashMap<String, Double>();
					for (String x : neworig.keySet()){
						oddsratio.put(x, roundToSignificantFigures(getRatioLogRank(x), 5));
					}
					
					newahsd.add(neworig);
					newahsd.add(newadj);
					newahsd.add(newclas);
					newahsd.add(patn);
					newahsd.add(oddsratio);

				}
				else{
					HashMap<String, Double> patn = new HashMap<String, Double>();
					for (String x : neworig.keySet()){
						patn.put(x, (double) getNumPatients(x));
					}
					
					HashMap<String, Double> oddsratio = new HashMap<String, Double>();
					for (String x : neworig.keySet()){
						oddsratio.put(x, roundToSignificantFigures(getRatioFisher(x), 5));
					}
					
					newahsd.add(neworig);
					newahsd.add(newadj);
					newahsd.add(patn);
					newahsd.add(oddsratio);
				
				}
				
				hah.put(newahsd, inputhah.get(ahsd));
			}
			output.put(s, hah);
			
		}
		
		
		
		
		
		return output;
	}
	
	
	public static double roundToSignificantFigures(double num, int n) {
		if (Double.isNaN(num) || Double.isInfinite(num)){
			return Double.NaN;
		}
		
	    if(num == 0) {
	        return 0;
	    }

	    final double d = Math.ceil(Math.log10(num < 0 ? -num: num));
	    final int power = n - (int) d;

	    final double magnitude = Math.pow(10, power);
	    final long shifted = Math.round(num*magnitude);
	    return shifted/magnitude;
	}
	
	
	public String seedAtBeginning(String s, String t){
		if (t.equals("none")){
			return "none";
		}
		String[] splitted = t.split(":");
		HashSet<String> hs = new HashSet<String>();
		for (int i=0; i<splitted.length; i++){
			hs.add(splitted[i]);
		}
		
		hs.remove(s);
		String returnVal = s;
		
		for (String x : hs){
			returnVal = returnVal + ":" + x;
		}
		
		return returnVal;
	}
	
	
	public boolean checkConditions(String s, String t){
		String[] genes1 = s.split(":");
		String[] genes2 = t.split(":");
		
		HashSet<String> g1 = new HashSet<String>();
		HashSet<String> g2 = new HashSet<String>();
		
		//System.out.println("s:" + s);
		//System.out.println("t:" + t);
		
		
		for (int i=0; i<genes1.length; i++){
			g1.add(genes1[i]);
		}
		
		for (int i=0; i<genes2.length; i++){
			g2.add(genes2[i]);
		}
		
		boolean sSubsetOfT = true;
		for (String x : g1){
			if (!g2.contains(x)){
				sSubsetOfT = false;
				//System.out.println(x + " is not in g2 but is in g1");
			}
		}
		
		if (sSubsetOfT == true){
			return true;
		}
		
		/*
		HashSet<String> p1 = new HashSet<String>();
		HashSet<String> p2 = new HashSet<String>();
		for (int i=0; i<filteredSampleValues.size(); i++){
			if (g1.contains(filteredSampleValues.get(i)[0])){
				p1.add(filteredSampleValues.get(i)[1]);
			}
			if (g2.contains(filteredSampleValues.get(i)[0])){
				p2.add(filteredSampleValues.get(i)[1]);
			}
		}
		
		if (p1.equals(p2)){
			if (genes2.length >= genes1.length){
				return true;
			}
		}
		*/
		
		return false;
	}
	
	
	public int getNumPatients(String genes){
		int result = 0;
		String[] splitted = genes.split(":");
		HashSet<String> checker = new HashSet<String>();
		for (int i=0; i<splitted.length; i++){
			checker.add(splitted[i]);
		}
		
		HashSet<String> patients = new HashSet<String>();
		for (int i=0; i<filteredSampleValues.size(); i++){
			if (checker.contains(filteredSampleValues.get(i)[0])){
				if (filteredSampleValues.get(i)[1]!="no_sample"){
					patients.add(filteredSampleValues.get(i)[1]);
				}

			}
		}
		
		result = patients.size();
		return result;
	}
	

	public double getRatioLogRank(String genes){
		String[] g = genes.split(":");

		
		HashSet<String> gs = new HashSet<String>();
		for (int i=0; i<g.length; i++){
			gs.add(g[i]);
		}
		
		
		HashSet<String> inModulePatients = new HashSet<String>();
		HashSet<String> outOfModulePatients = new HashSet<String>();
		
		for (int i=0; i<filteredSampleValues.size(); i++){
			if (gs.contains(filteredSampleValues.get(i)[0])){
				inModulePatients.add(filteredSampleValues.get(i)[1]);
			}
			else{
				outOfModulePatients.add(filteredSampleValues.get(i)[1]);
			}
		}
		
		int inModuleVar1 = 0;
		int outOfModuleVar1 = 0;
		
		for (int i=0; i<clinicalValues.size(); i++){
			if (inModulePatients.contains(clinicalValues.get(i)[0])){
				if (clinicalValues.get(i)[1].toLowerCase().equals("alive") ||
					clinicalValues.get(i)[1].toLowerCase().equals("yes") ||	
					clinicalValues.get(i)[1].toLowerCase().equals("y") ||	
					clinicalValues.get(i)[1].toLowerCase().equals("0") ||
					clinicalValues.get(i)[1].toLowerCase().equals("living")){
					inModuleVar1 ++;
					
					
				}
			}
			else{
				if (clinicalValues.get(i)[1].toLowerCase().equals("alive") ||
						clinicalValues.get(i)[1].toLowerCase().equals("yes") ||	
						clinicalValues.get(i)[1].toLowerCase().equals("y") ||	
						clinicalValues.get(i)[1].toLowerCase().equals("0") ||
						clinicalValues.get(i)[1].toLowerCase().equals("living")){
						outOfModuleVar1 ++;
						
						
				}
				
			}
			
		}
		
		double p1 = inModuleVar1/ (double) inModulePatients.size();
		double p2 = outOfModuleVar1/ (double) outOfModulePatients.size();
		double rvalue = p1*(1-p2)/(double) (p2*(1-p1));

		
		if (!Double.isNaN(rvalue) && !Double.isInfinite(rvalue)){
			rvalue = Math.log(rvalue);
		}
		return rvalue;
		
		
	}
	

	
	//variable 1 is the first one to appear in otherValues
	public double getRatioFisher(String genes){
		String[] g = genes.split(":");
		
		String v1 = clinicalValues.get(0)[1];
		
		HashSet<String> gs = new HashSet<String>();
		for (int i=0; i<g.length; i++){
			gs.add(g[i]);
		}
		
		
		HashSet<String> inModulePatients = new HashSet<String>();
		HashSet<String> outOfModulePatients = new HashSet<String>();
		
		for (int i=0; i<filteredSampleValues.size(); i++){
			if (gs.contains(filteredSampleValues.get(i)[0])){
				inModulePatients.add(filteredSampleValues.get(i)[1]);
			}
			else{
				outOfModulePatients.add(filteredSampleValues.get(i)[1]);
			}
		}
		
		int inModuleVar1 = 0;
		int outOfModuleVar1 = 0;
		
		for (int i=0; i<clinicalValues.size(); i++){
			if (inModulePatients.contains(clinicalValues.get(i)[0])){
				if (clinicalValues.get(i)[1].equals(v1)){
					inModuleVar1++;
				}
			}
			else{
				if (clinicalValues.get(i)[1].equals(v1)){
					outOfModuleVar1++;
				}
			}
		}
		
		double p1 = inModuleVar1/ (double) inModulePatients.size();
		double p2 = outOfModuleVar1/ (double) outOfModulePatients.size();
		double rvalue = p1*(1-p2)/(double) (p2*(1-p1));
		
		if (!Double.isNaN(rvalue) && !Double.isInfinite(rvalue)){
			rvalue = Math.log(rvalue);
		}

	
		return rvalue;
	}
	
	
	
	public void exportMostCorrelated(){
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		System.out.println("Please enter a p-Value cutoff (between 0 and 1)");
		double i = 0.05;
		try {
			i = Double.valueOf(br.readLine());
			if (i<0 || i>1){
				System.out.println("Please enter a p-Value cutoff between 0 and 1");
				return;
			}
		} catch (NumberFormatException e) {
			System.out.println("Please enter a p-Value cutoff between 0 and 1");
			return;
		} catch (IOException e) {
			System.out.println("Please enter a p-Value cutoff between 0 and 1");
			return;
		}
		
		ArrayList<HashMap<String, Double>> a = extractMostCorrelated(i);

		HashMap<String, Double> mostCorrelated = a.get(0);
		HashMap<String, Double> mostCorrelatedFDR = a.get(1);
		
		final String lineSep = System.getProperty("line.separator");
		String fileName = null;
		FileWriter fout = null;
		
		System.out.println("Please enter a valid path to save the file");
		try {
			fileName = br.readLine();
		} catch (IOException e1) {
			e1.printStackTrace();
			return;
		}

		try {
				//fileName = file.getAbsolutePath();
				if (!fileName.substring(fileName.length()-4,fileName.length()).equals(".csv")){
					fileName = fileName + ".csv";
				}
				fout = new FileWriter(fileName);
				
				fout.write("# HyperModules Results" + lineSep);
				fout.write("# Date: " + ',' + DateFormat.getDateTimeInstance().format(new Date()) + lineSep + lineSep);
				
				fout.write("# Shuffle Number: "+ ',' + shuffleNumber + lineSep);
				fout.write("# Statistical Test: " + ','+ stat + lineSep + lineSep);
				
				fout.write(lineSep);
				
				fout.write("Module" + ',' + "Statistical Test P-Value" + ',' + "FDR P-Value" + lineSep);
				for (String s : mostCorrelated.keySet()){
					fout.write(s + ',' + mostCorrelated.get(s) + ',' + mostCorrelatedFDR.get(s) + lineSep);
				}

			

		} 
		catch (IOException e) {
			JOptionPane.showMessageDialog(null,
										  e.toString(),
										  "Error Writing to \"" + fileName + "\"",
										  JOptionPane.ERROR_MESSAGE);
		} 
		catch (Exception e){
			e.printStackTrace();
		}
		finally {
			if (fout != null) {
				try {
					fout.close();
					System.out.println("file saved.");
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		
	}
	
	
	
}
