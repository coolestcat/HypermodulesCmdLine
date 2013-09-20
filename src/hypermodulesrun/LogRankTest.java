package hypermodulesrun;


import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import JSci.maths.statistics.NormalDistribution;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
/**
 * My implementation of the log rank test for correlation between two sets of survival data with censorship
 * The data is stored in ArrayDeques for faster addition and removal
 * @author alvinleung
 *
 */
public class LogRankTest {

	private HashSet<Double> allTimeHash;
	private ArrayList<Double> uniqueTimes;
	private int hashSize;
	
	/**
	 * Constructor
	 */
	public LogRankTest(double[] followupDays){
		super();
		this.allTimeHash = new HashSet<Double>();
		for (int i=0; i<followupDays.length; i++){
			allTimeHash.add(followupDays[i]);
		}
		this.hashSize = allTimeHash.size();
		
		//sorted list of unique times 
		this.uniqueTimes = new ArrayList<Double>();
		for (Double add : allTimeHash){
			uniqueTimes.add(add);
		}
		Collections.sort(uniqueTimes);
	}

	/**
	 * The main test method
	 * @param time1 - first set of times, already sorted into ascending order
	 * @param time2 - second set of times, already sorted into ascending order
	 * @param censor1 - censorship data corresponding to time1 - 1 means deceased, 0 means censored
	 * @param censor2 - censorship data corresponding to time2 - 1 means deceased, 0 means censored
	 * @return a double containing the test statistic [0], the variance[1], and the pValue[2] of running the log rank test
	 */
	public Double[] logRank(ArrayDeque<Double> sortedTime1, ArrayDeque<Double> sortedTime2, ArrayDeque<Double> sortedCensor1, ArrayDeque<Double> sortedCensor2){

		/*
		ArrayDeque<Double> sortedTime1 = new ArrayDeque<Double>();
		ArrayDeque<Double> sortedTime2 = new ArrayDeque<Double>();
		ArrayDeque<Double> sortedCensor1 = new ArrayDeque<Double>();
		ArrayDeque<Double> sortedCensor2 = new ArrayDeque<Double>();
		
		for (int i=0; i<time1.length; i++){
			//allTimeHash.add(time1[i]);
			sortedTime1.add(time1[i]);
			sortedCensor1.add(censor1[i]);
		}
		
		for (int i=0; i<time2.length; i++){
			//allTimeHash.add(time2[i]);
			sortedTime2.add(time2[i]);
			sortedCensor2.add(censor2[i]);
		}

		 */
		
		/*
		
		//time/censor map with repeated keys allowed (multimap)
		Multimap<Double, Double> timeOne = ArrayListMultimap.create();
		Multimap<Double, Double> timeTwo = ArrayListMultimap.create();
		
		//sort time1 and time 2 into arraylists, WITH repetition
		ArrayList<Double> sortedTime1 = new ArrayList<Double>();
		ArrayList<Double> sortedTime2 = new ArrayList<Double>();
		
		//allTimeHash removes duplicates (unique time values in time1 AND time2)
		HashSet<Double> allTimeHash = new HashSet<Double>();
		//allTimes includes ALL time values in time1 and time2
		
		//add values from time1 and time2 to the lists
		ArrayList<Double> allTimes = new ArrayList<Double>();
		for (int i=0; i<time1.length; i++){
			allTimes.add(time1[i]);
			allTimeHash.add(time1[i]);
			timeOne.put(time1[i], censor1[i]);
			sortedTime1.add(time1[i]);

		}
		
		for (int i=0; i<time2.length; i++){
			allTimeHash.add(time2[i]);
			allTimes.add(time2[i]);
			timeTwo.put(time2[i], censor2[i]);
			sortedTime2.add(time2[i]);
		}
		
		
		//sort time1 and time2
		Collections.sort(sortedTime1);
		Collections.sort(sortedTime2);
		
		
		//match up censors with sorted arraylist
		ArrayList<Double> sortedCensor1 = new ArrayList<Double>();
		ArrayList<Double> sortedCensor2 = new ArrayList<Double>();
		
		//map <time, boolean> whether or not censors corresponding to time value in multimap already added
		HashMap<Double, Boolean> alreadyAdded = new HashMap<Double, Boolean>();
		
		//initialize each value to false
		for (Double key : timeOne.keySet()){
			alreadyAdded.put(key,  false);
		}
		
		//loop through sortedTime, and for each unique value in sortedTime, we get corresponding collection of censor values from TimeOne/TimeTwo multimap and add it to sortedCensor
		for (int i=0; i<sortedTime1.size(); i++){
			if (alreadyAdded.get(sortedTime1.get(i))==false){
				Collection<Double> coll = timeOne.get(sortedTime1.get(i));
				for (Double value : coll){
					sortedCensor1.add(value);
				}
			}
			alreadyAdded.put(sortedTime1.get(i), true);

		}
		
		//same thing for time2
		HashMap<Double, Boolean> alreadyAdded2 = new HashMap<Double, Boolean>();
		
		
		for (Double key : timeTwo.keySet()){
			alreadyAdded2.put(key,  false);
		}
		
		for (int i=0; i<sortedTime2.size(); i++){
			if (alreadyAdded2.get(sortedTime2.get(i))==false){
				Collection<Double> coll = timeTwo.get(sortedTime2.get(i));
				for (Double value : coll){
					sortedCensor2.add(value);
				}
			}
			
			alreadyAdded2.put(sortedTime2.get(i), true);
		}

		 */
		
		//number still at risk (alive) at the beginning of time i in group 1 (before event occurs)
		double[] n1i = new double[hashSize];
		
		//number still at risk (alive) at the beginning of time i in group 2 (before event occurs)
		double[] n2i = new double[hashSize];
		
		//number of deaths in group 1 at time i (DOES NOT INCLUDE CENSORED)
		double[] o1i = new double[hashSize];
		
		//number of deaths in group 2 at time i (DOES NOT INCLUDE CENSORED)
		double[] o2i = new double[hashSize];

		//TODO: delete: Collections.sort(allTimes);
		
		//number of group 1 still alive after time i
		double n1 = sortedTime1.size();
		//number of group 2 still alive after time i
		double n2 = sortedTime2.size();

		//int n1ix = 0;
		//int n2ix = 0;
		
		//loop through sorted list of unique times
		for (int i=0; i<uniqueTimes.size(); i++){
			
			//initialize number of observed deaths
			 double o1 = 0.0;
			 double o2 = 0.0;
			 
			 //set number still alive at the beginning of period i
			 n1i[i] = n1;
			 n2i[i] = n2;
			
			 //while the first element of sortedTime1 matches the current uniqueTime value, 
			 //we increment number of observed deaths if the corresponding censor is 1;
			 //in both cases (0 or 1) we decrease the number still alive (n1) and remove the first element from sortedTime1 and sortedCensor1
			while(!sortedTime1.isEmpty() && uniqueTimes.get(i).equals(sortedTime1.peekFirst())){

				if (sortedCensor1.remove()==1){
					
					o1++;
					 n1=n1-1;
					 sortedTime1.poll();
					
				}
				else{
					sortedTime1.poll();
					n1=n1-1;
				}
			 }
			 
			//same thing for sortedTime2 and sortedCensor2
			 while(!sortedTime2.isEmpty() && uniqueTimes.get(i).equals(sortedTime2.peekFirst())){
				 if (sortedCensor2.remove()==1){
					 o2++;
					 n2=n2-1;
					 sortedTime2.poll();

				 }
				 
				 else{
					 sortedTime2.poll();
					 n2=n2-1;
					 
				 }
			 }

		//set observed values
		o1i[i]=o1;
		o2i[i]=o2;
		}

		
		/*
		for (int i=0; i<n1i.length; i++){
			System.out.println(n1i[i]);
		}
		System.out.println("|");
		for(int i=0; i<o1i.length; i++){
			System.out.println(o1i[i]);
		}
		System.out.println("|");
		for (int i=0; i<n2i.length; i++){
			System.out.println(n2i[i]);
		}
		System.out.println("|");
		for (int i=0; i<o2i.length; i++){
			System.out.println(o2i[i]);
		}
		*/
		
		double ni = 0.0;
		double oi = 0.0;
		double totalVariance = 0.0;
		double statistic = 0.0;
		
		//loop through uniqueTimes, get test statistic and variance
		//ni = number alive in both group 1 and group 2 at the beginning of time i (before event occurs)
		//oi = observed deaths in both groups at time i
		//statistic = sum from 1 to (number of unique times) of (o1i-e1i), e1i=(oi/ni)*n1i
		//variance = sum of oi*(n1i/ni)*(1-n1i/ni)*(ni-oi)/(ni-1)
		
		for (int i=0; i<uniqueTimes.size(); i++){
			ni = n1i[i]+n2i[i];
			oi = o1i[i]+o2i[i];
			
			statistic += (o1i[i]- ((oi/ni)*n1i[i]));
			if ((ni-1)!=0){
				totalVariance += ((oi*(n1i[i]/ni)*(1-(n1i[i]/ni))*(ni-oi))/(ni-1));
			}

		}
		
		//pValue for normal distribution
        double pValue = 2 * (1 - new NormalDistribution().cumulative(
                Math.abs(statistic / Math.pow(totalVariance, 0.5))));
		
        /*
		System.out.println(statistic);
		System.out.println(totalVariance);
		System.out.println(pValue);
		*/
		
		Double[] returnValue = new Double[3];
		returnValue[0]=statistic;
		returnValue[1]=totalVariance;
		returnValue[2]=pValue;
		return returnValue;
		
	}
	
	
	/**
	 * test whether it is high survival or low survival - expected > observed --> low, expected < observed --> high 
	 * @param sortedTime1
	 * @param sortedTime2
	 * @param sortedCensor1
	 * @param sortedCensor2
	 * @return
	 */
	public boolean logRankSurvivalTest(ArrayDeque<Double> sortedTime1, ArrayDeque<Double> sortedTime2, ArrayDeque<Double> sortedCensor1, ArrayDeque<Double> sortedCensor2){
		double[] n1i = new double[hashSize];
		double[] n2i = new double[hashSize];
		double[] o1i = new double[hashSize];
		double[] o2i = new double[hashSize];

		double n1 = sortedTime1.size();
		double n2 = sortedTime2.size();

		for (int i=0; i<uniqueTimes.size(); i++){
			

			 double o1 = 0.0;
			 double o2 = 0.0;

			 n1i[i] = n1;
			 n2i[i] = n2;
			 
			while(!sortedTime1.isEmpty() && uniqueTimes.get(i).equals(sortedTime1.peekFirst())){

				if (sortedCensor1.remove()==1){
					
					o1++;
					 n1=n1-1;
					 sortedTime1.poll();
					
				}
				else{
					sortedTime1.poll();
					n1=n1-1;
				}
			 }
			 
			 while(!sortedTime2.isEmpty() && uniqueTimes.get(i).equals(sortedTime2.peekFirst())){
				 if (sortedCensor2.remove()==1){
					 o2++;
					 n2=n2-1;
					 sortedTime2.poll();

				 }
				 
				 else{
					 sortedTime2.poll();
					 n2=n2-1;
					 
				 }
			 }

		o1i[i]=o1;
		o2i[i]=o2;
		}


		double c = 0;
		for (int i=0; i<o1i.length; i++){
			c += o1i[i];
		}
		
		double d = 0;
		for (int i=0; i<uniqueTimes.size(); i++){
			double oi = o1i[i] + o2i[i];
			double ni = n1i[i] + n2i[i];
			d += ((oi/ni)*n1i[i]);
		}
		
		if (d<c){
			return false;
		}
		else{
			return true;
		}
	}

}

