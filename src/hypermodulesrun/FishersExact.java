package hypermodulesrun;


import java.util.ArrayList;
import com.google.common.math.BigIntegerMath;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
/**
 * 
 * @author alvinleung
 * Implementation of Fisher's Exact Test using the recursive method for enumeration described by Boulton and Wallace in "Occupancy of a Rectangular Array (1973)"
 *
 */
public class FishersExact {
	public int[][] thisMatrix;
	private int[] n;
	private int[] v;
	private int[] d;
	private int[] r;
	private int N;
	private int nRows;
	private int nColumns;
	private int[] columnTotals;
	private ArrayList<Double> pValues;
	
	/**
	 * constructor
	 * @param thisMatrix - input matrix 
	 */
	public FishersExact(int[][] thisMatrix){
		this.thisMatrix = thisMatrix;
	}
	
	//test main method
	/*
	public static void main(String[] args){
		int[][] thisMatrix = {{0,4}, {0, 150}, {2,29}};
		//int[][] thisMatrix = {{0,4}, {0,15}, {0,35}, {3,70}};
		//int[][] thisMatrix = {{5,0}, {1,4}};
		FishersExact fe = new FishersExact(thisMatrix);
		double value = fe.fisher2c();
		System.out.println(value);
	}
	*/

	/**
	 * pretty print of 2 dimensional matrix
	 * @param mat
	 */
	public static void printMatrix(int[][] mat){
		System.out.print("\n");
		for (int[] is : mat) {
			System.out.print("[");
			for (int i = 0; i < is.length; i++) {
				System.out.print(((i == 0) ? "" : ", ") + is[i]);
			}
			System.out.print("]\n");
		}
		System.out.print("\n");
	}
	
	/**
	 * runs fisher's exact test on 2 columns and x rows. 
	 * @return pValue
	 */
	public double fisher2c(){
		double p = 0;
		
		int nRows = thisMatrix.length;
		int nColumns = thisMatrix[0].length;

		int[] rowTotals = new int[nRows];
		int[] columnTotals =  new int[nColumns]; 
		for (int i=0; i<nRows; i++){
			rowTotals[i]=findSumRow(i, nColumns);
		}
		
		for (int i=0; i<nColumns; i++){
			columnTotals[i]=findSumColumn(i, nRows);
		}
		
		int N = 0;
		for (int i=0; i<rowTotals.length; i++){
			N += rowTotals[i];
		}
		
		int[] n = new int[nRows+1];
		int[] d = new int[nRows];
		int[] v = new int[nRows];
		int[] r = new int[nRows];
		
		n[0]=0;
		for (int i=1; i<n.length; i++){
			n[i]=rowTotals[i-1];
		}
		
		d[0]=0;
		v[0]=columnTotals[0];
		r[0]=N;
		
		
		this.n = n;
		this.r = r;
		this.v = v;
		this.d = d;
		this.nRows = nRows;
		this.nColumns = nColumns;
		this.columnTotals = columnTotals; 
		this.N = N;
		
		 
		this.pValues = new ArrayList<Double>();
		int level = 1;
		enumerate(level);
		
		for (int i=1; i<nRows; i++){
			d[i] = thisMatrix[i-1][0];
		}


		Double original = calculateStatistic();

		for (int i=0; i<pValues.size(); i++){
			if (pValues.get(i)<=original){
				p = p + pValues.get(i);
			}
		}

		p = (double)Math.round(p * 100000000) / 100000000;
		return p;
	}
	

	/**
	 * recursive method for enumeration - adapted from the FORTRAN code found in Boulton and Wallace's paper
	 * on matrix enumeration (finding all possible matrices with the row and column totals provided)
	 * @param level
	 */
	private void enumerate(int level){
		if (level==this.nRows){
			this.pValues.add(calculateStatistic());
			return;
		}
		v[level] = v[level-1] - d[level-1];
		r[level] = r[level-1] - n[level];
		int hibound;
		int lobound;
		if (n[level]<v[level]){
			hibound = n[level];
		}
		else{
			hibound = v[level];
		}
		
		if (v[level]-r[level]>0){
			lobound = v[level]-r[level];
		}
		else{
			lobound = 0;
		}
		//System.out.println("v: " + level + ":" + v[level]);
		//System.out.println("r: " + level + ":" + r[level]);
		
		for (int i=lobound; i<=hibound; i++){
			d[level] = i;
			enumerate(level+1);
		}
	}

	/**
	 * implements the formula for fisher's exact test
	 * @return the testStatistic for a given enumeration matrix
	 */
	private Double calculateStatistic(){
		int[][] matrix = new int[nRows][nColumns];
		int sum = 0;
		for (int i=0; i<nRows-1; i++){
			sum = sum + d[i+1];
			matrix[i][0] = d[i+1];
		}
		BigInteger numerator = new BigInteger("1");
		BigInteger denominator = BigIntegerMath.factorial(N);
		matrix[nRows-1][0] = columnTotals[0] - sum;
		
		for (int i=0; i<nRows; i++){
			matrix[i][1] = n[i+1]-matrix[i][0];
			numerator = numerator.multiply(BigIntegerMath.factorial(n[i+1]));
			denominator = denominator.multiply(BigIntegerMath.factorial(matrix[i][0]).multiply(BigIntegerMath.factorial(matrix[i][1])));
		}
		

		
		for (int i=0; i<nColumns; i++){
			numerator = numerator.multiply(BigIntegerMath.factorial(columnTotals[i]));
		}
		BigDecimal newNum = new BigDecimal(numerator);
		BigDecimal newDenom = new BigDecimal(denominator);

		BigDecimal statistic = newNum.divide(newDenom, 15, BigDecimal.ROUND_DOWN);
		
		//printMatrix(matrix);
		
		return statistic.doubleValue();
	}
	
	/**
	 * find the columntotal of a column
	 * @param column
	 * @param nRows
	 * @return
	 */
	private int findSumColumn(int column, int nRows){
		int sum = 0;
		for (int i=0; i<nRows; i++){
			sum += thisMatrix[i][column];
		}
		return sum;
		
	}
	
	/**
	 * find the rowtotal of a row
	 * @param row
	 * @param nColumns
	 * @return
	 */
	private int findSumRow(int row, int nColumns){
		int sum = 0;
		for (int i=0; i<nColumns; i++){
			sum += thisMatrix[row][i];
		}
		
		return sum;
	}
	
	
}
