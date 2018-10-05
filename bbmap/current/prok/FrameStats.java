package prok;

import dna.AminoAcid;
import shared.Tools;
import structures.ByteBuilder;

public class FrameStats {
	
	public FrameStats(String name_, int k_, int frames_, int leftOffset_){
		name=name_;
		k=k_;
		frames=frames_;
		kMax=1<<(2*k);
		invFrames=1.0f/frames;
		leftOffset=leftOffset_;
		
//		System.err.println(name+", "+frames);
		
		probs=new float[frames][kMax];
		countsTrue=new long[frames][kMax];
		countsFalse=new long[frames][kMax];
		counts=new long[][][] {countsFalse, countsTrue};
//		assert(!name.toLowerCase().contains("start")) : k+", "+frames+", "+leftOffset+", "+GeneModel.startRightOffset()+", "+GeneModel.startFrames();
	}
	
	public void add(int kmer, int frame, int valid){
		counts[valid][frame][kmer]++;
		validSums[valid]++;
	}
	
	public void add(FrameStats fs){
		assert(fs.k==k);
		assert(fs.frames==frames);

		Tools.add(counts, fs.counts);
		Tools.add(validSums, fs.validSums);
		calculate();
	}
	
	public void calculate(){
		average=(float)((validSums[1]+1.0)/(validSums[0]+validSums[1]+1.0));
		invAvg=1.0f/average;
		
		for(int a=0; a<frames; a++){
			for(int b=0; b<kMax; b++){
				long t=countsTrue[a][b];
				long f=countsFalse[a][b];
				probs[a][b]=(float)(t/(t+f+1.0))*invAvg;
			}
		}
	}

	public void parseData(byte[] line) {
		int a=0, b=0;
		final int valid, frame;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		valid=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		frame=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		final long[] row=counts[valid][frame];
		long sum=0;
		for(int kmer=0; kmer<row.length; kmer++){
			while(b<line.length && line[b]!='\t'){b++;}
			assert(b>a) : "Missing field 1: "+new String(line);
			long count=Tools.parseInt(line, a, b);
			b++;
			a=b;
			row[kmer]=count;
			sum+=count;
		}
		validSums[valid]+=sum;
	}
	
	@Override
	public String toString(){
		return appendTo(new ByteBuilder()).toString();
	}
	
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append('#').append(name).nl();
		bb.append("#valid\tframe");
		for(int i=0; i<kMax; i++){bb.tab().append(AminoAcid.kmerToString(i, k));}
		bb.nl();
		for(int a=0; a<2; a++){
			for(int b=0; b<frames; b++){
				bb.append(a);
				bb.tab().append(b);
				for(int c=0; c<kMax; c++){
					bb.tab().append(counts[a][b][c]);
				}
				bb.nl();
			}
		}
		return bb;
	}

	public final String name;
	public final int k;
	public final int frames;
	public final int kMax;
	public final float invFrames;
	public final int leftOffset;
	
	public final float[][] probs;
	public final long[][] countsTrue;
	public final long[][] countsFalse;
	public final long[][][] counts;

	public final long[] validSums=new long[2];
	private float average=-1;
	private float invAvg=-1;
	
}
