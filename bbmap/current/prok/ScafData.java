package prok;

import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import stream.Read;
import structures.IntList;

class ScafData {
	
	ScafData(Read r){
		this(r.id, r.bases, new byte[r.length()]);
	}
	
	ScafData(String name_, byte[] bases_, byte[] frames_){
		name=name_;
		bases=bases_;
		frames=frames_;
		gLines[0]=new ArrayList<GffLine>();
		gLines[1]=new ArrayList<GffLine>();
	}
	
	void clear(){
		Arrays.fill(frames, (byte)0);
		starts.clear();
		stops.clear();
	}
	
	void reverseComplement(){
		AminoAcid.reverseComplementBasesInPlace(bases);
		strand=1^strand;
	}
	
	void add(GffLine gline){
		assert(gline.strand>=0) : gline+"\n"+gline.strand;
		gLines[gline.strand].add(gline);
	}
	
	int strand(){return strand;}

	public int length() {return bases==null ? 0 : bases.length;}
	
	final String name;
	final byte[] bases;
	final byte[] frames;
	final IntList starts=new IntList(8);
	final IntList stops=new IntList(8);
	private int strand=0;
	
	/** gLines[strand] holds the GffLines for that strand */
	@SuppressWarnings("unchecked")
	ArrayList<GffLine>[] gLines=new ArrayList[2];
}
