package hypermodulesrun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

public class MainRun {
	
	public static void main(String[] args){
		if (args.length!=5){
			System.out.println("Please enter correct number of arguments.");
			System.out.println("USAGE: java -jar [*.jar] [PATH_TO_NETWORK_INTERACTIONS] [PATH_TO_SAMPLE_VALUES] [PATH_TO_CLINICAL_VALUES] [SHUFFLE_NUMBER] [STAT_TEST}");
		}
		else{
			
			ArrayList<String[]> network = readCSVfile(args[0]);
			ArrayList<String[]> sampleValues = readCSVfile(args[1]);
			ArrayList<String[]> clinicalValues = readCSVfile(args[2]);
			int shuffleNumber = Integer.valueOf(args[3]);
			String stat = args[4];
			
			AlgorithmTask at = new AlgorithmTask(network, sampleValues, clinicalValues, shuffleNumber, stat );
			at.run();
		}

	}

	 public static ArrayList<String[]> readCSVfile(String path) {
	   	     ArrayList<String[]> Rs = new ArrayList<String[]>();
		     String[] OneRow;
		     
	            try {
	            BufferedReader brd = new BufferedReader (new FileReader(path));

	            String st = brd.readLine();
	            
	            	while (st!=null) {
	            			OneRow = st.split(",|\\s|;|\t");
	            			Rs.add(OneRow);
	            			st = brd.readLine();
	                } 
		            brd.close();
	            } 
	            catch (Exception e) {
	                String errmsg = e.getMessage();
	                System.out.println ("File not found:" +errmsg);
	            }       
	        return Rs;
	        }

	
}

