import java.util.ArrayList;

/************************************************************************************************** 
 * @author Elliot Ensink
 * Basic global alignment tool for two DNA sequences
 *************************************************************************************************/	
public class SeqAlign 
{
	private int idScore, subScore, indelScore;
	private Integer[][] scoreMatrix;
	private String seqA, seqB;
	private ArrayList<coordinate> path;

	public SeqAlign(String A, String B, int id, int sub, int indel)
	{
		idScore = id;//Score for a nucleotide match
		subScore = sub;//Penalty score for a nucleotide mismatch
		indelScore = indel;//Penalty score for an insertion/deletion
		seqA = "-"+A;//DNA sequence A
		seqB = "-"+B;//DNA sequence B
		scoreMatrix = new Integer[seqA.length()][seqB.length()];
		path = new ArrayList<coordinate>();
	}

	private void align()
	{
		setupMatrix();
		path = findPath(0,0);
		printMatrix();
		deconvolutePath();

	}

	private ArrayList<coordinate> findPath(int i, int j)
	{
		coordinate coord = new coordinate(i,j);
		path.add(coord);
		Integer caseA, caseB, caseC;
		caseA = caseB = caseC = null;
		if((j+1)<seqB.length())
		{
			caseA = scoreMatrix[i][j] + indelScore;
		}
		if((i+1)<seqA.length())
		{
			caseB = scoreMatrix[i][j] + indelScore;
		}
		if((j+1)<seqB.length()&&(i+1)<seqA.length())
		{
			if(seqA.charAt(i) == seqB.charAt(j))
				caseC = scoreMatrix[i][j] + idScore;
			else
				caseC = scoreMatrix[i][j] + subScore;
		}
		Integer[] pathOptions = {caseA,caseB,caseC};
		pathOptions = sortPathOptions(pathOptions);

		for(Integer option : pathOptions)
		{
			if(option != null)
			{
				if(option == caseA)
				{
					if(scoreMatrix[i][j+1] == null || caseA >= scoreMatrix[i][j+1])
					{
						scoreMatrix[i][j+1] = caseA;
						findPath(i,j+1);
					}
				}
				else if(option == caseB)
				{
					if(scoreMatrix[i+1][j] == null || caseB >= scoreMatrix[i+1][j])
					{
						scoreMatrix[i+1][j] = caseB;
						findPath(i+1,j);
					}
				}
				else if(option == caseC)
				{
					if(scoreMatrix[i+1][j+1] == null ||caseC >= scoreMatrix[i+1][j+1])
					{
						scoreMatrix[i+1][j+1] = caseC;
						findPath(i+1,j+1);
					}
				}
			}
		}
		return path;
	}

	/**************************************************************************************************
	 * firstRowColSetup
	 * Intializes the first row and column of the score matrix to contain the appropriate scores
	 * according to the rules:
	 * scoreMatrix(i,0) = i*indelscore
	 * scoreMatrix(0,j) = j*indelscore
	 *************************************************************************************************/
	private void setupMatrix()
	{
		for(int i=0;i<scoreMatrix[0].length;i++)
		{	
			for(int j=0;j<scoreMatrix.length;j++)
			{
				if(j == 0)
					scoreMatrix[i][j] = i*indelScore;
				else if(i == 0)
					scoreMatrix[i][j] = j*indelScore;
				else
					scoreMatrix[i][j] = null;
			}
		}

	}
	/**************************************************************************************************
	 * printMatrix
	 * 
	 * prints the matrix of scores generated during the alignment process
	 *************************************************************************************************/
	private void printMatrix()
	{
		int i;
		int j;
		String spacing = "\t";
		System.out.print("\t\t");
		for(j=0;j<seqB.length();j++)
		{
			System.out.print(j+spacing);
		}
		System.out.println();
		System.out.print("\t\t");
		for(j=0;j<seqB.length();j++)
		{
			System.out.print(seqB.charAt(j)+spacing);
		}
		System.out.println();
		for(i=0;i<seqA.length();i++)
		{
			for(j=0;j<seqB.length()+2;j++)
			{
				if(j==0)
					System.out.print(i+spacing);
				else if(j==1)
					System.out.print(seqA.charAt(i)+spacing);
				else
					System.out.print(scoreMatrix[i][j-2]+spacing);
			}
			System.out.println();
		}
	}

	private Integer[] sortPathOptions(Integer[] pathOpt)
	{
		for(int i = 0; i < pathOpt.length; i++)
		{
			for(int j = 0; j < pathOpt.length-1; j++)
			{
				if(pathOpt[j] == null)
				{
					pathOpt[j] = pathOpt[j+1];
					pathOpt[j+1] = null;
				}
				else if(pathOpt[j+1] != null && pathOpt[j] < pathOpt[j+1])
				{
					int temp = pathOpt[j+1];
					pathOpt[j+1] = pathOpt[j];
					pathOpt[j] = temp;
				}
			}
		}
		return pathOpt;
	}

	private void deconvolutePath()
	{
		for(coordinate c : path)
		{
			System.out.println(c);
		}
	}

	private class coordinate
	{
		public int i;
		public int j;
		public coordinate(int i, int j)
		{
			this.i = i;
			this.j = j;
		}
		private boolean equals(coordinate c)
		{
			if(this.i == c.i && this.j == c.j)
				return true;
			else
				return false;

		}
		public String toString()
		{
			return "("+this.i+","+this.j+")";
		}


	}

	public static void main(String[] args)
	{
		SeqAlign seq = new SeqAlign("ATCGT","TGGTG",1,-1,-2);
		seq.align();
	}


}

