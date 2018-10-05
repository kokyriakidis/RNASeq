package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import dna.AminoAcid;
import shared.Tools;
import stream.Read;

public class GeneCaller {
	
	GeneCaller(int minLen_, int maxOverlapSameStrand_, int maxOverlapOppositeStrand_, 
			float minStartScore_, float minStopScore_, float minInnerScore_, float minOrfScore_){
		this(minLen_, maxOverlapSameStrand_, maxOverlapOppositeStrand_, minStartScore_, minStopScore_, minInnerScore_, minOrfScore_, null);
	}
	
	GeneCaller(int minLen_, int maxOverlapSameStrand_, int maxOverlapOppositeStrand_, 
			float minStartScore_, float minStopScore_, float minInnerScore_, float minOrfScore_, GeneModel pgm_){
		minLen=minLen_;
		maxOverlapSameStrand=maxOverlapSameStrand_;
		maxOverlapOppositeStrand=maxOverlapOppositeStrand_;
		pgm=pgm_;
		
		minStartScore=minStartScore_;
		minStopScore=minStopScore_;
		minInnerScore=minInnerScore_;
		minOrfScore=minOrfScore_;
	}
	
	public ArrayList<Orf> callGenes(Read r){
		return callGenes(r, pgm);
	}
	
	public ArrayList<Orf> callGenes(Read r, GeneModel pgm_){
		pgm=pgm_;
		
		final String name=r.id;
		final byte[] bases=r.bases;
		ArrayList<Orf>[] frameLists=makeOrfs(name, bases, minLen);
		ArrayList<Orf>[] brokenLists=breakOrfs(frameLists, bases);
		
		//This block will just print all orfs.
//		ArrayList<Orf> temp=new ArrayList<Orf>();
//		for(ArrayList<Orf> broken : brokenLists){
//			temp.addAll(broken);
//		}
//		Collections.sort(temp);
//		return temp;
		
		
		ArrayList<Orf> path=findPath(brokenLists, bases);
		return path;
		
//		ArrayList<Orf> selectedOrfs=findPath(brokenLists); //TODO
//		return selectedOrfs;
	}
	
	/** 
	 * Generates lists of all max-length non-overlapping Orfs per frame.
	 * There IS overlap between frames.
	 * All Orfs come out flipped to + orientation. 
	 * */
	public static ArrayList<Orf>[] makeOrfs(String name, byte[] bases, int minlen){
		@SuppressWarnings("unchecked")
		ArrayList<Orf>[] array=new ArrayList[6];
		for(int strand=0; strand<2; strand++){
			for(int frame=0; frame<3; frame++){
				ArrayList<Orf> list=makeOrfsForFrame(name, bases, frame, strand, minlen);
				array[frame+3*strand]=list;
				if(strand==1){
					for(Orf orf : list){
						assert(orf.frame==frame);
						assert(orf.strand==strand);
						orf.flip();
					}
				}
			}
			AminoAcid.reverseComplementBasesInPlace(bases);
		}
		return array;
	}
	
	public ArrayList<Orf> findPath(ArrayList<Orf>[] frameLists, byte[] bases){
		ArrayList<Orf> all=new ArrayList<Orf>();
		for(ArrayList<Orf> list : frameLists){all.addAll(list);}
		if(all.isEmpty()){return all;}
		Collections.sort(all);
		
		for(Orf orf : all){
			orf.pathScore=-9999;
		}
		
		int[] lastPositionScored=new int[6];
		Arrays.fill(lastPositionScored, -1);
//		for(int i=0; i<frameLists.length; i++){
//			ArrayList<Orf> list=frameLists[i];
//			if(list.isEmpty()){
//				Orf orf=new Orf("FAKE_ORF", i, i+5, i/3, i%3, bases, true);
//				orf.orfScore=-1;
//				list.add(orf);
//			}
//			Orf orf=list.get(0);
//			orf.pathScore=orf.orfScore;
//		}
		
		int[] bestIndex=new int[6];
		Orf[] bestInFrame=new Orf[6];
		for(Orf orf : all){
			final int myListNum=3*orf.strand+orf.frame;
			calcPathScore(orf, frameLists, lastPositionScored, bestIndex);
			if(bestInFrame[myListNum]==null || orf.pathScore>=bestInFrame[myListNum].pathScore){
				bestInFrame[myListNum]=orf;
				bestIndex[myListNum]=lastPositionScored[myListNum];
				assert(frameLists[myListNum].get(lastPositionScored[myListNum])==orf);
			}
		}
		
		Orf best=bestInFrame[0];
		for(int i=1; i<bestInFrame.length; i++){
			if(best==null || (bestInFrame[i]!=null && bestInFrame[i].pathScore>best.pathScore)){
				best=bestInFrame[i];
			}
		}
		ArrayList<Orf> bestPath=new ArrayList<Orf>();
		for(Orf orf=best; orf!=null; orf=orf.prev){
			bestPath.add(orf);
		}
		Collections.sort(bestPath);
		return bestPath;
	}
	
	private float calcPathScore(Orf orf, ArrayList<Orf>[] frameLists, int[] lastPositionScored, int[] bestIndex){
		final int myListNum=3*orf.strand+orf.frame;

//		System.err.println("* "+orf);
//		System.err.println("* "+Arrays.toString(lastPositionScored));
//		System.err.println();
		
		for(int strand=0; strand<2; strand++){
			for(int frame=0; frame<3; frame++){
				int listNum=frame+3*strand;
				 ArrayList<Orf> list=frameLists[listNum];
				 int lastPos=lastPositionScored[listNum];
				 int bestPos=bestIndex[listNum];
				 calcPathScore(orf, list, lastPos, bestPos);
			}
		}
		
//		System.err.println(myListNum+", "+Arrays.toString(lastPositionScored)+", "+frameLists[myListNum].size());
		
		lastPositionScored[myListNum]++;
		assert(frameLists[myListNum].get(lastPositionScored[myListNum])==orf) : myListNum+"\n"+orf+"\n"+frameLists[myListNum].get(lastPositionScored[myListNum])+"\n"
			+Arrays.toString(lastPositionScored)+"\n"+frameLists[myListNum].get(lastPositionScored[myListNum]+1);
		assert(orf.prev!=null || orf.stop<100000);
//		assert(orf.pathScore>-10) : orf.pathScore+"\n"+orf+"\n"+orf.prev+"\n";
		return orf.pathScore;
	}
	
	private float calcPathScore(Orf orf, ArrayList<Orf> list, int lastPos, int bestPos){
		if(lastPos<0){
			if(orf.prev==null){
				orf.pathScore=orf.orfScore;
			}
			return orf.pathScore;
		}
		if(list.isEmpty()){return -1;}
		
		
		final boolean sameFrame=(orf.strand==list.get(0).strand && orf.frame==list.get(0).frame);
//		System.err.println("\nExamining   \t"+orf+"\nlastPos="+lastPos+", bestPos="+bestPos+", sameFrame="+sameFrame);
		boolean found=false;
		for(int i=lastPos, min=Tools.max(0, bestPos-80); i>=min || (i>0 && !found); i--){
			Orf prev=list.get(i);
			assert(prev!=orf) : prev;
//			System.err.println("Comparing to \t"+prev);
			int maxOverlap=(orf.strand==prev.strand ? maxOverlapSameStrand : maxOverlapOppositeStrand);
			if(orf.isValidPrev(prev, maxOverlap)){
				int overlap=Tools.max(0, prev.stop-orf.start+1);
				float orfScore=overlap==0 ? orf.orfScore : orf.calcOrfScore(overlap);
				if(prev.strand!=orf.strand){orfScore=(orfScore-20);}
				float pathScore=prev.pathScore+orfScore;
				
//				System.err.println("valid, overlap="+overlap+", score="+orfScore+", pathScore="+pathScore+", prev="+prev.pathScore);
				
				if(overlap<1 && prev.pathScore>0){found=true;}
				if(pathScore>=orf.pathScore){
					orf.pathScore=pathScore;
					orf.prev=prev;
//					System.err.println("Set as best");
				}
			}
			if(found && prev.stop<maxOverlap-2000 && orf.prev!=null){
				System.err.println("Breaking");
				break;
			}
		}
		return orf.pathScore;
	}
	
	/** 
	 * Generates a list of maximal-length Orfs only (non-overlapping).
	 * All Orfs come out in native orientation (unflipped). 
	 * */
	public static ArrayList<Orf> makeOrfsForFrame(String name, byte[] bases, int startFrame, int strand, int minlen){//123 TODO
//		assert(false) : "TODO";
		assert(minlen>=3);
		if(bases==null || bases.length<minlen){return null;}
		ArrayList<Orf> orfs=new ArrayList<Orf>();
//		int mask=63;
		int code=0;
		int start=-2;
		int frame=0;
		int pos=startFrame;
		
		
		for(; pos<bases.length; pos++){
			byte b=bases[pos];
			int x=AminoAcid.baseToNumber[b];
//			code=((code<<2)|x)&mask;
			code=((code<<2)|x);
			frame++;
			if(frame==3){
				frame=0;
				if(start>=0){
					if(GeneModel.isStopCodon(code) || code<0){//NOTE: This adds a stop codon wherever there are Ns.
						int len=pos-start+1;
						if(len>=minlen){
							Orf f=new Orf(name, start, pos, strand, startFrame, bases, true);
							orfs.add(f);
						}
						start=-1;
					}
				}else{
					if(start==-2 || (start<0 && GeneModel.isStartCodon(code))){
						start=pos-2;
					}
				}
				code=0;
			}
		}

		//Add a stop codon at the sequence end.
		if(start>=0){
			pos--;
			while(frame!=3 && frame!=-1){
				pos--;
				frame--;
			}
			int len=pos-start+1;
			if(len>=minlen){
				assert(pos<bases.length) : start+", "+pos+", "+bases.length;
				Orf f=new Orf(name, start, pos, strand, startFrame, bases, true);
				orfs.add(f);
			}
		}
		
		return orfs;
	}
	
	public ArrayList<Orf>[] breakOrfs(ArrayList<Orf>[] frameLists, byte[] bases){

		@SuppressWarnings("unchecked")
		ArrayList<Orf>[] brokenLists=new ArrayList[6];
		for(int strand=0; strand<2; strand++){
			for(int frame=0; frame<3; frame++){
				int fnum=frame+3*strand;
				ArrayList<Orf> longest=frameLists[fnum];
				ArrayList<Orf> broken=new ArrayList<Orf>();
				for(Orf orf : longest){
					assert(orf.frame==frame);
					assert(orf.strand==strand);
					ArrayList<Orf> temp=breakOrf(orf, bases);
					broken.addAll(temp);
				}
				Collections.sort(broken);
				brokenLists[fnum]=broken;
			}
			AminoAcid.reverseComplementBasesInPlace(bases);
		}
		return brokenLists;
	}
	
	public ArrayList<Orf> breakOrf(Orf longest, byte[] bases){
		assert(longest.start<longest.stop);
		final int flipped=longest.flipped();
		if(flipped==1){longest.flip();}//Now the orf is aligned to its native strand
		
//		IntList starts=new IntList();
//		starts.add(longest.start);
		FrameStats stats=pgm.innerKmerStats;
		
		final String name=longest.scafName;
		final int start=longest.start;
		final int stop=longest.stop;
		final int strand=longest.strand;
		final int max=Tools.min(longest.stop-2, longest.stop-minLen+3);
		final int k=pgm.innerKmerLength;
		final int mask=~((-1)<<(2*k));

		final float stopScore=pgm.calcStopScore(longest.stop, bases);
		
		ArrayList<Orf> broken=new ArrayList<Orf>();
		int created=0;
		
		int codon=0;
		int kmer=0;
		int len=0;
		int numKmers=0;
		float currentScore=0;
		for(int pos=start, currentFrame=0; pos<=stop; pos++){
			final byte b=bases[pos];
			final int x=AminoAcid.baseToNumber[b];
			codon=((codon<<2)|x);
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				if(len>=k){
					float prob=stats.probs[currentFrame][kmer];
					float dif=prob-0.99f;//Prob above 1 is more likely than average
					currentScore+=dif;
//					outstream.println("pos="+pos+"\tdif="+String.format("%.4f", dif)+",\tscore="+String.format("%.4f", currentScore)+
//							"\tasStart="+String.format("%.4f", pgm.calcStartScore(pos-2, bases))+"\tasStop="+String.format("%.4f", pgm.calcStopScore(pos, bases))+
//							"\tcodon="+AminoAcid.kmerToString(kmer, 3)+" frame="+(currentFrame));
				}else{
//					outstream.println("pos="+pos+"\tdif="+String.format("%.4f", 0.0)+",\tscore="+String.format("%.4f", 0.0)+
//							"\tasStart="+String.format("%.4f", pgm.calcStartScore(pos-2, bases))+"\tasStop="+String.format("%.4f", pgm.calcStopScore(pos, bases))+
//							"\tcodon="+AminoAcid.kmerToString(kmer, 3)+" frame="+(currentFrame));
				}
			}else{
				len=0;
				kmer=0;
			}
			
			currentFrame++;
//			outstream.println("pos="+pos+", codon="+AminoAcid.kmerToString(kmer, 3)+", frame="+currentFrame+", start="+start+", isStartCodon="+pgm.isStartCodon(codon));
			if(currentFrame==3){
				currentFrame=0;
				if(pos<max && created<breakLimit && (pos==start+2 || pgm.isStartCodon(codon))){
//					outstream.println(x);
					int glen=stop-pos+3;
					assert(glen>=minLen) : "glen="+glen+", minLen="+minLen+", pos="+pos+", max="+max+", start="+start;
					
					int oStart=pos-2;
					float startScore=pgm.calcStartScore(oStart, bases);
					
					if(startScore>=minStartScore || stopScore>=minStopScore){
						Orf orf=new Orf(name, pos-2, stop, strand, longest.frame, bases, false);
						orf.kmerScore=currentScore;
						orf.startScore=startScore;
						orf.stopScore=stopScore;
						
//						orf.orfScore=orf.calcScore();
//						assert(false) : orf.kmerScore+", "+currentScore+", "+orf;

						assert(orf.frame==longest.frame);
						assert(orf.strand==strand);

						if(strand==1){orf.flip();}
						broken.add(orf);
						created++;

//						outstream.println("strand="+strand+", frame="+longest.frame);
//						outstream.println(new String(bases, longest.start, 3)+" longest; start "+longest.start);
//						outstream.println(new String(bases, orf.start, 3)+" broken; start "+orf.start);
//						outstream.println(new String(bases, 337-1, 3)+" ncbi; start "+(337-1));
//						outstream.println(new String(bases, longest.start-3, 9));
//						outstream.println();
//						outstream.println(new String(bases, longest.stop-2, 3)+" longest; stop "+longest.stop);
//						outstream.println(new String(bases, orf.stop-2, 3)+" broken; stop "+orf.stop);
//						outstream.println(new String(bases, (2799-1)-2, 3)+" ncbi; stop "+(2799-1));
//						outstream.println(new String(bases, longest.stop-5, 9));
//						outstream.println();
//
//						outstream.println("\nlongest: "+longest+"\nbroken:  "+orf+"\n"+
//								"glen="+glen+", minLen="+minLen+", pos="+pos+", max="+max+", start="+start+"\n");
//						assert(pos<500) : new String(bases, orf.start-3, orf.length()+6);
					}
				}
				codon=0;
			}
		}
		
		int removed=0;
		for(int i=0; i<broken.size(); i++){//Fix scores because they were generated together, from start to stop, to make this O(N) instead of O(N^2).
			Orf orf=broken.get(i);
			orf.kmerScore=currentScore-orf.kmerScore;
			orf.orfScore=orf.calcOrfScore();
			if(orf.averageKmerScore()<minInnerScore || orf.orfScore<minOrfScore){
				broken.set(i, null);
				removed++;
			}
		}
		if(removed>0){
			Tools.condenseStrict(broken);
		}
		
//		{
//			for(int i=0; i<broken.size(); i++){
//				Orf orf=broken.get(i);
//				if(orf.startScore>-0.3 || orf.stopScore>-0.2){broken.set(i, null);}
//			}
//			Tools.condenseStrict(broken);
//			
//			for(Orf orf : broken){
//				assert(orf.startScore<0);
//				assert(orf.stopScore<0);
//			}
//		}
		
		if(flipped==1){longest.flip();}
		return broken;
	}
	
	public ArrayList<Orf> findPath(ArrayList<Orf>[] lists){
		assert(false) : "TODO";
		ArrayList<Orf> selectedOrfs=null;
		return selectedOrfs;
	}
	
	static boolean[] makeIsCodon(String[] codons){
		boolean[] array=new boolean[64];
		for(String s : codons){
			int x=AminoAcid.toNumber(s);
			array[x]=true;
		}
		return array;
	}

	GeneModel pgm;
	final int minLen;
	final int maxOverlapSameStrand;
	final int maxOverlapOppositeStrand;
	final float minStartScore;
	final float minStopScore;
	final float minInnerScore;
	final float minOrfScore;
	
	public static int breakLimit=12;
	
	private static PrintStream outstream=System.err;
	
}
