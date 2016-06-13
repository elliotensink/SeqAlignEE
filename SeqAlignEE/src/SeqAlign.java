import java.util.ArrayList;

/************************************************************************************************** 
 * @author Elliot Ensink
 * Basic global alignment tool for two DNA sequences
 *************************************************************************************************/	
public class SeqAlign 
{
	private int idScore, subScore, indelScore; 
	private Integer maxAlignScore;
	private Integer[][] scoreMatrix;
	private String seqA, seqB, alignedSeqA, alignedSeqB;
	private ArrayList<coordinate> path;
	private ArrayList<coordinate> bestPath;
	private ArrayList<ArrayList<coordinate>> finalPaths;

	public SeqAlign(String A, String B, int id, int sub, int indel)
	{
		idScore = id;//Score for a nucleotide match
		subScore = sub;//Penalty score for a nucleotide mismatch
		indelScore = indel;//Penalty score for an insertion/deletion
		maxAlignScore = null;//Final alignment score
		seqA = "-"+A;//DNA sequence A
		seqB = "-"+B;//DNA sequence B
		alignedSeqA = "";//DNA sequence A post alignment
		alignedSeqB = "";//DNA sequence B post alignment

		//Matrix of alignment scores
		scoreMatrix = new Integer[seqA.length()][seqB.length()];

		//Contains all path segments found during alignment
		finalPaths = new ArrayList<ArrayList<coordinate>>();

		//Contains "cheapest" path to reach the end node
		bestPath = new ArrayList<coordinate>();

		//Current path through matrix that is updated during align()
		path = new ArrayList<coordinate>();
	}

	/**********************************************************************************************
	 * Aligns two DNA sequences
	 **********************************************************************************************/
	private void align()
	{
		setupMatrix();
		findPath(0,0, null);
		printMatrix();
		coordinate endNode = new coordinate(scoreMatrix.length-1,scoreMatrix[0].length-1);
		deconvolutePath(endNode);
	}
	
	/********************************************************************************************** 
	 * Recursively steps through matrix to find all paths through matrix starting from (0,0) and
	 * ending at (searchMatrix.length,searchMatrix[0].length)
	 * @param i vertical coordinate
	 * @param j horizontal coordinate
	 * @param prev previous node/coordinate visited
	 *********************************************************************************************/
	private void findPath(int i, int j, coordinate prev)
	{
		//printMatrix();
		coordinate coord = new coordinate(i,j);
		PathOpt caseA, caseB, caseC;
		caseA = new PathOpt(null,'A');
		caseB = new PathOpt(null,'B');
		caseC = new PathOpt(null,'C');
		path.add(coord);//Add coordinate to current path
		
		//Calculate scores for the 3 cases (right,down, diagonal down and right)
		if((j+1)<seqB.length())
		{
			caseA.score = scoreMatrix[i][j] + indelScore;
		}
		if((i+1)<seqA.length())
		{
			caseB.score = scoreMatrix[i][j] + indelScore;
		}
		if((j+1)<seqB.length()&&(i+1)<seqA.length())
		{
			if(seqA.charAt(i+1) == seqB.charAt(j+1))
				caseC.score = scoreMatrix[i][j] + idScore;
			else if(seqA.charAt(i+1) == '-' || seqB.charAt(j+1) == '-')
				caseC.score = scoreMatrix[i][j] + indelScore;
			else
				caseC.score = scoreMatrix[i][j] + subScore;
		}
		
		//Sort cases, largest case will be chosen first
		PathOpt[] options = new PathOpt[]{caseA,caseB,caseC};
		options = sortPathOptions(options);
		
		//For each case, make the next move, update matrix, and continue recursively
		for(int index = 0; index < 3; index++)
		{
			PathOpt opt = options[index];
			if(opt != null)
			{
				if(opt.letter == 'A' && opt.score != null)
				{
					if(scoreMatrix[i][j+1] == null || opt.score >= scoreMatrix[i][j+1])
					{
						scoreMatrix[i][j+1] = opt.score;
						findPath(i,j+1,coord);
					}
				}
				else if(opt.letter == 'B' && opt.score != null)
				{
					if(scoreMatrix[i+1][j] == null || opt.score >= scoreMatrix[i+1][j])
					{
						scoreMatrix[i+1][j] = opt.score;
						findPath(i+1,j,coord);
					}
				}
				else if(opt.letter == 'C' && opt.score != null)
				{
					if(scoreMatrix[i+1][j+1] == null || opt.score >= scoreMatrix[i+1][j+1])
					{
						scoreMatrix[i+1][j+1] = opt.score;
						findPath(i+1,j+1,coord);
					}
				}
			}
		}
		//Once end is reached add the path segment to the list of all path segments
		if(maxAlignScore == null || scoreMatrix[i][j] >= maxAlignScore)
		{
			finalPaths.add(new ArrayList<coordinate>(path));
			path.clear();
			path.add(prev);
		}
		return;
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
		for(int i=0;i<scoreMatrix.length;i++)
		{	
			for(int j=0;j<scoreMatrix[0].length;j++)
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
		System.out.println("-----------------------------------------------------");
	}

	/**********************************************************************************************
	 * Sorts the list of possible cases (path options) largest to smallest
	 * @param options List of cases (path options) possible
	 * @return sorted list of path options
	 *********************************************************************************************/
	private PathOpt[] sortPathOptions(PathOpt[] options)
	{
		//Bubble sort
		for(int i = 0; i < options.length; i++)
		{
			for(int j = 0; j < options.length-1; j++)
			{
				if(options[j].score == null)
				{
					PathOpt temp = new PathOpt(options[j].score,options[j].letter);
					options[j] = options[j+1];
					options[j+1] = temp;
				}
				else if(options[j+1].score != null && options[j].score < options[j+1].score)
				{
					PathOpt temp = new PathOpt(options[j+1].score,options[j+1].letter);
					options[j+1] = options[j];
					options[j] = temp;
				}
			}
		}
		return options;
	}
	/**********************************************************************************************
	 * Deconvolutes path segments to find the best path through the matrix
	 * @param end Used to determine the end point of the matrix (bottom right corner)
	 *********************************************************************************************/
	private void deconvolutePath(coordinate end)
	{
		bestPath = separateEndPaths(end);
		alignedSeqA = "";
		alignedSeqB = "";
		for(int i = bestPath.size()-1; i > 0; i--)
		{
			coordinate c1 = bestPath.get(i);
			coordinate c2 = bestPath.get(i-1);
			if(c2.i < c1.i && c2.j < c1.j)
			{
				alignedSeqA = seqA.charAt(c1.i)+alignedSeqA;
				alignedSeqB = seqB.charAt(c1.j)+alignedSeqB;
			}
			else if(c2.i < c1.i)
			{
				alignedSeqA = seqA.charAt(c1.i)+alignedSeqA;
				alignedSeqB = "-"+alignedSeqB;
			}
			else if(c2.j < c1.j)
			{
				alignedSeqA = "-"+alignedSeqA;
				alignedSeqB = seqB.charAt(c1.j)+alignedSeqB;
			}
			else
				System.out.println("Error in deconvoluting");
		}
		System.out.println("Sequence A: "+ alignedSeqA);
		System.out.println("Sequence B: "+ alignedSeqB);
	}

	/**********************************************************************************************
	 * Removes unnecessary path segments (ones that don't lead to the cheapest route)
	 * @param end Used to determine the end point of the matrix (bottom right corner)
	 * @return path - the cheapes path through the matrix
	 *********************************************************************************************/
	private ArrayList<coordinate> separateEndPaths(coordinate end)
	{
		ArrayList<coordinate> path = new ArrayList<coordinate>();
		for(ArrayList<coordinate> p : finalPaths)
		{
			for(coordinate c : p)
			{
				if(c.equals(end))
				{
					if(path.isEmpty())
					{
						path = new ArrayList<coordinate>(p);
					}
					else
					{
						coordinate begin = p.get(0);
						boolean containsNode = false;
						for(coordinate n : path)
						{
							if(n.equals(begin))
							{
								containsNode = true;
								break;
							}
						}
						if(containsNode)
						{
							int i = path.size()-1;
							while(!path.get(i).equals(begin))
							{
								path.remove(i);
								i--;
							}
							path.remove(i);//remove the beginning node of next segment from path
							path.addAll(p);
						}
						//else the path is not possible
					}
				}
			}
		}
		return path;
	}

	//Structure for keeping track of the path through the matrix
	private class coordinate
	{
		public int i;//Vertical
		public int j;//Horizontal
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

	//Structure for organizing the possible path options
	private class PathOpt
	{
		public Integer score;
		public char letter;//A, B, or C
		public PathOpt(Integer score, char letter)
		{
			this.score = score;
			this.letter = letter;
		}
		public String toString()
		{
			return "("+this.letter+","+this.score+")";
		}
	}

	public static void main(String[] args)
	{
		SeqAlign seq = new SeqAlign("ATCGT","TGGTG",1,-1,-2);
		//SeqAlign seq = new SeqAlign("AGACTAGTTAC","AGACTTAC",1,-1,-2);
		seq.align();
	}


}

