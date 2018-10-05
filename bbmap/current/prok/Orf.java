package prok;

import dna.AminoAcid;
import dna.Feature;
import dna.Gene;
import shared.Tools;
import structures.ByteBuilder;

/**
 * @author bbushnell
 * @date Sep 20, 2018
 *
 */
public class Orf extends Feature {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** 
	 * Bases and coordinates are assumed to be the correct strand.
	 * Minus-strand ORFs can be flipped at the end of the constructor.
	 * @param scafName_
	 * @param start_
	 * @param stop_
	 * @param strand_
	 * @param bases
	 */
	public Orf(String scafName_, int start_, int stop_, int strand_, int frame_, byte[] bases, boolean flip) {
		super(scafName_, start_, stop_, strand_);
		frame=frame_;
		scafLen=bases.length;
		startCodon=getCodon(start, bases);
		stopCodon=getCodon(stop-2, bases);
		
		if(flip && strand==Gene.MINUS){flip();}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Init Helpers         ----------------*/
	/*--------------------------------------------------------------*/
	
	//Assumes bases are in the correct strand
	private static int getCodon(int from, byte[] bases){
		int codon=0;
		for(int i=0; i<3; i++){
//			assert(i+from<bases.length) : i+", "+from+", "+bases.length;
			byte b=bases[i+from];
			int x=AminoAcid.baseToNumber[b];
			codon=(codon<<2)|x;
		}
		return codon;
	}
	
//	//Assumes bases are in the correct strand
//	private static float calcKmerScore(int stop, byte[] bases){
//		return 1;//TODO
//	}

	public float calcOrfScore(){
		return calcOrfScore(0);
	}

	public float calcOrfScore(int overlap){
		double e=0.08;
		double e2=e*1.1;
		double a=Math.sqrt(Tools.max(e2, e+startScore));
		double b=Math.sqrt(Tools.max(e2, e+stopScore));
		double c=Tools.max(e2, e+averageKmerScore());
		assert(a!=Double.NaN);
		assert(b!=Double.NaN);
		assert(c!=Double.NaN);
		c=4*Math.pow(c, 2.2);
		double d=(0.1*a*b*c*(Math.pow(length()-overlap, 2.5)-(overlap<1 ? 0 : Math.pow(overlap+50, 2))));
		if(d>0){d=Math.sqrt(d);}
		assert(d!=Double.NaN);
		return (float)d;
	}
	
	public float averageKmerScore(){
		return kmerScore/(length()-GeneModel.innerKmerLength+1);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void flip(){
		flip(scafLen);
	}
	
	public boolean isValidPrev(Orf prev, int maxOverlap){
		if(prev.stop>=stop || prev.stop>=start+maxOverlap || prev.start>=start){return false;}
		if(prev.frame==frame && prev.strand==strand && prev.stop>=start){return false;}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return toGff();
	}
	
	public String toGff(){
		ByteBuilder bb=new ByteBuilder();
		appendGff(bb);
		return bb.toString();
	}
	
	public ByteBuilder appendGff(ByteBuilder bb){
		bb.append(scafName==null ? "." : scafName).append('\t');
		bb.append("BBTools").append('\t');
		bb.append("CDS").append('\t');
		bb.append(start+1).append('\t');
		bb.append(stop+1).append('\t');
		
		if(orfScore<0){bb.append('.').append('\t');}
		else{bb.append(orfScore, 2).append('\t');}
		
		bb.append(strand<0 ? '.' : Gene.strandCodes2[strand]).append('\t');
		
		bb.append('0').append('\t');

		//bb.append('.');
		bb.append("fr").append(frame).append(',');
//		bb.append(startCodon).append(',');
//		bb.append(stopCodon).append(',');
		bb.append("startScr:").append(startScore, 3).append(',');
		bb.append("stopScr:").append(stopScore, 3).append(',');
		bb.append("innerScr:").append(averageKmerScore(), 3).append(',');
		bb.append("len:").append(length()).append(',');
		bb.append("start:").append(AminoAcid.codonToString(startCodon)).append(',');
		bb.append("stop:").append(AminoAcid.codonToString(stopCodon));
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Overrides           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void flip(int scafLen_){//TODO: Move to Feature
		assert(scafLen_==scafLen);
		super.flip(scafLen_);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final int frame;
	public final int scafLen;

	public final int startCodon;
	public final int stopCodon;
	
	public float startScore;
	public float stopScore;
	public float kmerScore;
	
	public float orfScore;
	
	public float pathScore;
	public Orf prev;
	
//	//I need 6 of these.  They should be objects.
//	public float scoreForward=-1;
//	public float scoreReverse=-1;
//	
//	public Orf prevForward=null;
//	public Orf prevReverse=null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
}
