/************************************************************************************************** 
 * @author Elliot Ensink
 * Basic global alignment tool for two DNA sequences
 *************************************************************************************************/	
public class SeqAlign 
{
private int idScore, subScore, indelScore;
private int[][] scoreMatrix;
private String seqA, seqB;
	
public SeqAlign(String A, String B, int id, int sub, int indel)
{
	idScore = id;//Score for a nucleotide match
	subScore = sub;//Penalty score for a nucleotide mismatch
	indelScore = indel;//Penalty score for an insertion/deletion
	seqA = "-"+A;//DNA sequence A
	seqB = "-"+B;//DNA sequence B
	scoreMatrix = new int[seqA.length()][seqB.length()];
}

private void align()
{
	scoreMatrix[0][0] = 0;
	findPath(0,0);
	printMatrix();
	
}

private void findPath(int i, int j)
{
	boolean a, b, c; 
	a = b = c = false;
	int caseA, caseB, caseC;
	if((j+1)<seqB.length())
	{
		a = true;
		caseA = scoreMatrix[i][j] + indelScore;
		if(caseA >= scoreMatrix[i][j+1])
			scoreMatrix[i][j+1] = caseA;
	}
	if((i+1)<seqA.length())
	{
		b = true;
		caseB = scoreMatrix[i][j] + indelScore;
		if(caseB >= scoreMatrix[i+1][j])
			scoreMatrix[i+1][j] = caseB;
	}
	if(a&&b)
	{
		c = true;
		if(seqA.charAt(i) == seqB.charAt(j))
			caseC = scoreMatrix[i][j] + idScore;
		else
			caseC = scoreMatrix[i][j] + subScore;
		if(caseC >= scoreMatrix[i+1][j+1])
			scoreMatrix[i+1][j+1] = caseC;
	}
	if(a)
		findPath(i,j+1);
	if(c)
		findPath(i+1,j+1);
	if(b)
		findPath(i+1,j);
	else
		return;
	
	
}

//private void findPath()
//{
//	
//}


/**************************************************************************************************
 * firstRowColSetup
 * Intializes the first row and column of the score matrix to contain the appropriate scores
 * according to the rules:
 * scoreMatrix(i,0) = i*indelscore
 * scoreMatrix(0,j) = j*indelscore
 *************************************************************************************************/
//private void firstRowColSetup()
//{
//	for(int i=0;i<scoreMatrix[0].length;i++)
//			scoreMatrix[i][0] = i*indelScore;
//	for(int j=0;j<scoreMatrix.length;j++)
//		scoreMatrix[0][j] = j*indelScore;
//}

private void printMatrix()
{
	int i;
	int j;
	String spacing = "\t";
	//firstRowColSetup();
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

public static void main(String[] args)
{
	SeqAlign seq = new SeqAlign("ACGT","ACGT",1,-1,-2);
	seq.align();
}


}
