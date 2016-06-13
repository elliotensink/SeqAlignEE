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
	findPath(0,0);
	printMatrix();
	deconvolutePath();
	
}

private void findPath(int i, int j)
{
	coordinate coord = new coordinate(i,j);
	for(coordinate c : path)
	{
		if(c.equals(coord))
		{
			path.remove(c);
			break;
		}
	}
	path.add(coord);
	int caseA, caseB, caseC;
	if((j+1)<seqB.length())
	{
		caseA = scoreMatrix[i][j] + indelScore;
		if(scoreMatrix[i][j+1] == null || caseA >= scoreMatrix[i][j+1])
		{
			scoreMatrix[i][j+1] = caseA;
			findPath(i,j+1);
		}
	}
	if((i+1)<seqA.length())
	{
		caseB = scoreMatrix[i][j] + indelScore;
		if(scoreMatrix[i+1][j] == null || caseB >= scoreMatrix[i+1][j])
		{
			scoreMatrix[i+1][j] = caseB;
			findPath(i+1,j);
		}
	}
	if((j+1)<seqB.length()&&(i+1)<seqA.length())
	{
		if(seqA.charAt(i) == seqB.charAt(j))
			caseC = scoreMatrix[i][j] + idScore;
		else
			caseC = scoreMatrix[i][j] + subScore;
		if(scoreMatrix[i+1][j+1] == null ||caseC >= scoreMatrix[i+1][j+1])
		{
			scoreMatrix[i+1][j+1] = caseC;
			findPath(i+1,j+1);
		}
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

private void deconvolutePath()
{
	for(coordinate c : path)
	{
		System.out.println("("+c.i+","+c.j+")");
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
	
	
}

public static void main(String[] args)
{
	SeqAlign seq = new SeqAlign("ATCGT","TGGTG",1,-1,-2);
	seq.align();
}


}

