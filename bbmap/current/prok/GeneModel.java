package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

import dna.AminoAcid;
import dna.Gene;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;
import structures.IntList;

/**
 * This class is designed to analyze paired prokaryotic fna and gff files
 * to calculate the patterns in coding and noncoding frames, start and stop sites.
 * @author bbushnell
 * @date Sep 24, 2018
 *
 */
public class GeneModel {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public GeneModel(){
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(String genomeFname, String gffFname){
		fnames.add(ReadWrite.stripPath(genomeFname));
		FileFormat fnaFile=FileFormat.testInput(genomeFname, FileFormat.FA, null, true, true);
		FileFormat gffFile=FileFormat.testInput(gffFname, FileFormat.GFF, null, true, true);
		
		ArrayList<ScafData> scafList;
		{//Scoped to save memory
			ArrayList<Read> reads=ReadInputStream.toReads(fnaFile, maxReads);
			readsProcessed+=reads.size();
			scafList=new ArrayList<ScafData>(reads.size());
			for(Read r : reads){
				basesProcessed+=r.length();
				scafList.add(new ScafData(r));
			}
		}
		{//Scoped to save memory
			ArrayList<GffLine> gffLines=GffLine.loadGffFile(gffFile, "CDS");
			genesProcessed+=gffLines.size();
			HashMap<String, ScafData> scafMap=makeScafMap(scafList);
			fillScafData(gffLines, scafMap);
		}
		
		countBases(scafList);
		if(PROCESS_PLUS_STRAND){
			processStrand(scafList, Gene.PLUS);
		}
		if(PROCESS_MINUS_STRAND){
			for(ScafData sd : scafList){
				sd.clear();
				sd.reverseComplement();
			}
			processStrand(scafList, Gene.MINUS);
		}
	}
	
	public void add(GeneModel pgm){
		innerKmerStats.add(pgm.innerKmerStats);
		startStats.add(pgm.startStats);
		stopStats.add(pgm.stopStats);
		
		readsProcessed+=pgm.readsProcessed;
		basesProcessed+=pgm.basesProcessed;
		genesProcessed+=pgm.genesProcessed;
		
		fnames.addAll(pgm.fnames);
		taxIds.addAll(pgm.taxIds);
		Tools.add(baseCounts, pgm.baseCounts);
	}

	public static GeneModel loadModel(String fname) {
		ByteFile bf=ByteFile.makeByteFile1(fname, false);

		ArrayList<byte[]> headers=new ArrayList<byte[]>();
		ArrayList<byte[]> data=new ArrayList<byte[]>();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			assert(line.length>0) : line.length+", "+bf.lineNum();
			if(line[0]=='#'){
				headers.add(line);
				parseHeaderStatic(line);
			}else{
				data.add(line);
			}
		}
		GeneModel pgm=new GeneModel();
		for(byte[] line : headers){pgm.parseHeader(line);}
		pgm.parseData(data);
		return pgm;
	}
	
	private void parseData(ArrayList<byte[]> lines){
		FrameStats[] array=new FrameStats[] {innerKmerStats, startStats, stopStats};
		int i=0;
		for(FrameStats fs : array){
//			System.err.println("Loading "+fs.name+", "+fs.k+", "+fs.frames);
			for(int max=i+2*fs.frames; i<max; i++){
//				System.err.println("Processing line "+i);
				fs.parseData(lines.get(i));
			}
			fs.calculate();
		}
	}
	
	public static void parseHeaderStatic(byte[] line){
		
		assert(line[0]=='#');
		if(Tools.startsWith(line, "#k_inner")){
			int x=(int)parseNumber(line);
			assert(x==innerKmerLength);
			setInnerK(x);
		}else if(Tools.startsWith(line, "#k_end")){
			int x=(int)parseNumber(line);
			assert(x==endKmerLength);
			setEndK(x);
		}else if(Tools.startsWith(line, "#start_left_offset")){
			int x=(int)parseNumber(line);
			assert(x==startLeftOffset);
			setStartLeftOffset(x);
		}else if(Tools.startsWith(line, "#start_right_offset")){
			int x=(int)parseNumber(line);
			assert(x==startRightOffset);
			setStartRightOffset(x);
		}else if(Tools.startsWith(line, "#stop_left_offset")){
			int x=(int)parseNumber(line);
			assert(x==stopLeftOffset);
			setStopLeftOffset(x);
		}else if(Tools.startsWith(line, "#stop_right_offset")){
			int x=(int)parseNumber(line);
			assert(x==stopRightOffset);
			setStopRightOffset(x);
		}
	}
	
	public void parseHeader(byte[] line){
		
		assert(line[0]=='#');
		if(Tools.startsWith(line, "#scaffolds")){
			long x=parseNumber(line);
			readsProcessed=x;
		}else if(Tools.startsWith(line, "#bases")){
			long x=parseNumber(line);
			basesProcessed=x;
		}else if(Tools.startsWith(line, "#genes")){
			long x=parseNumber(line);
			genesProcessed=x;
		}else if(Tools.startsWith(line, "#ACGTN")){
			String[] split=new String(line).split("\t");
			for(int i=0; i<baseCounts.length; i++){
				baseCounts[i]=Long.parseLong(split[i+1]);
			}
		}
	}
	
	private static long parseNumber(byte[] line){
		return Tools.parseInt(line, Tools.indexOf(line, '\t')+1, line.length);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	HashMap<String, ScafData> makeScafMap(ArrayList<ScafData> scafList){
		HashMap<String, ScafData> scafMap=new HashMap<String, ScafData>(scafList.size()*3);
		for(ScafData sd : scafList){scafMap.put(sd.name, sd);}
		for(ScafData sd : scafList){
			String name=sd.name;
			int idx=name.indexOf(' ');
			if(idx>=0){
				String prefix=name.substring(0, idx);
				if(scafMap.containsKey(prefix)){
					assert(false) : "Duplicate degenerate name: '"+name+"', '"+prefix+"'";
				}else{
					scafMap.put(prefix, sd);
				}
			}
		}
		return scafMap;
	}
	
	public void fillScafData(ArrayList<GffLine> gffLines, HashMap<String, ScafData> scafMap){
		for(GffLine gline : gffLines){
			ScafData sd=scafMap.get(gline.seqid);
			assert(sd!=null) : "Can't find scaffold for GffLine "+gline.seqid;
			sd.add(gline);
		}
	}
	
	public void processStrand(ArrayList<ScafData> scafList, int strand){
		for(ScafData sd : scafList){
			ArrayList<GffLine> glines=sd.gLines[strand];
			for(GffLine gline : glines){
				assert(gline.strand==strand);
				processGene(gline, sd);
			}
			processFrames(sd.bases, sd.frames, innerKmerStats);
			BitSet startSet=processEnds(sd.bases, startStats, sd.starts, 1);
			BitSet stopSet=processEnds(sd.bases, stopStats, sd.stops, 1);
//			outstream.println("Processed "+sd.starts.size+" valid starts and "+sd.stops.size+" stops.");
			sd.clear();
			findStartCodons(sd.bases, sd.starts, startSet);
			findStopCodons(sd.bases, sd.stops, stopSet);
//			outstream.println("Found "+sd.starts.size+" invalid starts and "+sd.stops.size+" stops.");
			processEnds(sd.bases, startStats, sd.starts, 0);
			processEnds(sd.bases, stopStats, sd.stops, 0);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private void countBases(ArrayList<ScafData> scafList){
		for(ScafData sd : scafList){
			countBases(sd.bases);
		}
	}
	
	private void countBases(byte[] bases){
		for(byte b : bases){
			int x=AminoAcid.baseToNumberACGTother[b];
			baseCounts[x]++;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Finding Codons        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static void findStopCodons(byte[] bases, IntList list, BitSet valid){
		final int k=3;
		final int mask=~((-1)<<(2*k));
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					int point=i;//End of the stop codon
					if(isStopCodon(kmer) && !valid.get(point)){
						list.add(point);
					}
				}
			}else{len=0;}
		}
	}
	
	private static void findStartCodons(byte[] bases, IntList list, BitSet valid){
		final int k=3;
		final int mask=~((-1)<<(2*k));
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					int point=i-k+1;//Start of the start codon
					if(isStartCodon(kmer) && !valid.get(point)){
						list.add(point);
					}
				}
			}else{len=0;}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------     Processing GffLines      ----------------*/
	/*--------------------------------------------------------------*/
	
	private static void processGene(GffLine gline, ScafData sd){
		final int strand=gline.strand;
		assert(strand==sd.strand());
		final byte[] frames=sd.frames;
		int start=gline.start-1, stop=gline.stop-1;
		if(start<0 || stop>=sd.length()){return;}
		assert(start<stop);
		if(strand==Gene.MINUS){
			int x=sd.length()-start-1;
			int y=sd.length()-stop-1;
			start=y;
			stop=x;

//			String a=new String(sd.bases, start, 3);
//			String b=new String(sd.bases, stop-2, 3);
////			assert(false) : start+", "+stop+"\n"+gline+"\n"+new String(sd.bases, start, 3)+", "+new String(sd.bases, stop-2, 3);
//			outstream.println(a+", "+b+", "+start+", "+stop);
		}
		assert(start>=0) : gline.toString()+"\n"+sd.length()+"\n"+sd.name;
		markFrames(start, stop, frames);
		sd.starts.add(start);
		sd.stops.add(stop);
//		assert(gline.start!=337) : gline+"\n"+start+", "+stop;
	}
	
	/** 
	 * Each frame byte has a bit marked for valid coding frames.
	 * For example, if frames[23]=0b100, then base 23 is the 3rd base in a codon.
	 * If frames[23]=0, then base 23 is not coding on this strand.
	 * @param start
	 * @param stop
	 * @param frames
	 */
	private static void markFrames(int start, int stop, byte[] frames){
		assert(start<stop) : start+", "+stop;
		for(int i=start, frameBit=1; i<=stop; i++){
			frames[i]=(byte)(frames[i]|frameBit);
			frameBit<<=1;
			if(frameBit>4){frameBit=1;}
		}
//		assert(false) : Arrays.toString(Arrays.copyOfRange(frames, start, start+20))+"\n"+start; //This is correct
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Counting Kmers        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static void processFrames(byte[] bases, byte[] validFrames, FrameStats stats){
		final int k=stats.k;
		final int mask=~((-1)<<(2*k));
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length && i<k-1; i++){//Fill the first k-1 bases
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){len++;}
			else{len=0;}
		}
		
		for(int i=k-1, j=0; i<bases.length; i++, j++){
			byte b=bases[i];
			int vf=validFrames[j];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					for(int frame=0; frame<stats.frames; frame++){
						int valid=vf&1;
						stats.add(kmer, frame, valid);
						//For CDS start (0-based) of 189, i=192, j=189, vf=1, frame=0 - all as expected.
//						assert(valid==0) : "vf="+vf+", frame="+frame+", len="+len+", kmer="+AminoAcid.kmerToString(kmer, k)+", i="+i+", j="+j;
						vf=(vf>>1);
					}
				}
			}else{len=0;}
		}
	}
	
	private static BitSet processEnds(byte[] bases, FrameStats stats, IntList list, int valid){
		BitSet points=new BitSet(bases.length);
		for(int i=0; i<list.size; i++){
			int point=list.get(i);
			processPoint(bases, stats, list.get(i), valid);
			points.set(point);
		}
		return points;
	}
	
	private static void processPoint(byte[] bases, FrameStats stats, int point, int valid){
		final int k=stats.k;
		final int mask=~((-1)<<(2*k));
		int kmer=0;
		int len=0;

//		outstream.println("k="+k);
//		outstream.println("mask="+mask);
		
		int start=point-stats.leftOffset;
		
		int i=start, frame=0-k+1;
		while(i<0){i++; frame++;}
		for(; i<bases.length && frame<stats.frames; i++, frame++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
//			outstream.println("b="+(char)b+", kmer="+kmer+", len="+(len+1)+", frame="+frame);
			
			if(x>=0){
				len++;
				if(len>=k){
					stats.add(kmer, frame, valid);
				}
			}else{len=0;}
		}
//		assert(false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Scoring            ----------------*/
	/*--------------------------------------------------------------*/
	
	//Assumes bases are in the correct strand
	public float calcStartScore(int start, byte[] bases){
		float f=scorePoint(start, bases, startStats);
		return f;
	}
	
	//Assumes bases are in the correct strand
	public float calcStopScore(int stop, byte[] bases){
		float f=scorePoint(stop, bases, stopStats);
		return f;
	}
	
	public static float scorePoint(int point, byte[] bases, FrameStats stats){
		final int k=stats.k;
		final int mask=~((-1)<<(2*k));
		
		int kmer=0;
		int len=0;
		float score=0;

//		outstream.println("k="+k);
//		outstream.println("mask="+mask);
		
		int start=point-stats.leftOffset;
		for(int i=start, frame=0-k+1; i<bases.length && frame<stats.frames; i++, frame++){
			byte b=(i>=0 ? bases[i] : (byte)'A');
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
//			outstream.println("b="+(char)b+", kmer="+kmer+", len="+(len+1)+", frame="+frame);
			
			if(x>=0){
				len++;
				if(len>=k){
					float prob=stats.probs[frame][kmer];
					float dif=prob-0.99f;
					score+=dif;
					
//					if(stats.name.equals("startStats")){
//						System.err.println("frame="+frame+" kmer="+AminoAcid.kmerToString(kmer, k)+
//								String.format(" prob=%.4f\tdif=%.4f\tscore=%.4f", prob, dif, score)+
//								"\tvalid="+stats.counts[1][frame][kmer]+"\tinvalid="+stats.counts[0][frame][kmer]);
//					}
				}
			}else{len=0;}
		}
		
		//TODO: Add 0.1 to the score, and compress all negative scores into that 0.1?
//		if(stats.name.equals("startStats")){
//			System.err.println();
//			assert(score>0) : score+", "+stats.invFrames+", "+stats.name;
//		}
		
		return score*stats.invFrames;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           toString           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return appendTo(new ByteBuilder()).toString();
	}
	
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append("#BBMap "+Shared.BBMAP_VERSION_STRING+" Prokaryotic Gene Model\n");
		bb.append("#files");
		for(String fname : fnames){
			bb.tab().append(fname);
		}
		bb.nl();
		bb.append("#taxIDs");
		for(int i=0; i<taxIds.size; i++){
			bb.tab().append(taxIds.get(i));
		}
		bb.nl();
		bb.append("#k_inner\t").append(innerKmerLength).nl();
		bb.append("#k_end\t").append(endKmerLength).nl();
		bb.append("#start_left_offset\t").append(startLeftOffset).nl();
		bb.append("#start_right_offset\t").append(startRightOffset).nl();
		bb.append("#stop_left_offset\t").append(stopLeftOffset).nl();
		bb.append("#stop_right_offset\t").append(stopRightOffset).nl();
		bb.append("#scaffolds\t").append(readsProcessed).nl();
		bb.append("#bases\t").append(basesProcessed).nl();
		bb.append("#genes\t").append(genesProcessed).nl();
		bb.append("#GC\t").append(gc(),2).nl();
		bb.append("#ACGTN");
		for(long x : baseCounts){
			bb.tab().append(x);
		}
		bb.nl();

		innerKmerStats.appendTo(bb);
		startStats.appendTo(bb);
		stopStats.appendTo(bb);
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Stats             ----------------*/
	/*--------------------------------------------------------------*/

	public final FrameStats innerKmerStats=new FrameStats("innerKmerStats", innerKmerLength, 3, 0);
	public final FrameStats startStats=new FrameStats("startStats", endKmerLength, startFrames, startLeftOffset);
	public final FrameStats stopStats=new FrameStats("stopStats", endKmerLength, stopFrames, stopLeftOffset);
	
	/*--------------------------------------------------------------*/
	
	public ArrayList<String> fnames=new ArrayList<String>();
	public IntList taxIds=new IntList();
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	long readsProcessed=0;
	long basesProcessed=0;
	long genesProcessed=0;
	long[] baseCounts=new long[5];
	
	public float gc(){
		long a=baseCounts[0];
		long c=baseCounts[1];
		long g=baseCounts[2];
		long t=baseCounts[3];
		return (float)((g+c)/Tools.max(1.0, a+t+g+c));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	public static void setInnerK(int k){
		innerKmerLength=k;
		innerKmax=1<<(2*innerKmerLength);
	}
	
	public static void setEndK(int k){
		endKmerLength=k;
		endKmax=1<<(2*endKmerLength);
	}
	
	public static void setStartLeftOffset(int x){
		startLeftOffset=x;
		startFrames=startLeftOffset+startRightOffset+1;
//		System.err.println("startLeftOffset="+startLeftOffset+", startRightOffset="+startRightOffset+", frames="+startFrames);
	}
	
	public static void setStartRightOffset(int x){
		startRightOffset=x;
		startFrames=startLeftOffset+startRightOffset+1;
//		System.err.println("startLeftOffset="+startLeftOffset+", startRightOffset="+startRightOffset+", frames="+startFrames);
//		assert(false) : endLeftOffset+", "+endRightOffset+", "+endFrames;
	}
	
	public static void setStopLeftOffset(int x){
		stopLeftOffset=x;
		stopFrames=stopLeftOffset+stopRightOffset+1;
//		System.err.println("stopLeftOffset="+stopLeftOffset+", stopRightOffset="+stopRightOffset+", frames="+stopFrames);
	}
	
	public static void setStopRightOffset(int x){
		stopRightOffset=x;
		stopFrames=stopLeftOffset+stopRightOffset+1;
//		System.err.println("stopLeftOffset="+stopLeftOffset+", stopRightOffset="+stopRightOffset+", frames="+stopFrames);
//		assert(false) : endLeftOffset+", "+endRightOffset+", "+endFrames;
	}
	
	public static int innerKmerLength=4;
	public static int innerKmax=1<<(2*innerKmerLength);
	
	public static int endKmerLength=3;
	public static int endKmax=1<<(2*endKmerLength);

	static int startLeftOffset(){return startLeftOffset;}
	static int startRightOffset(){return startRightOffset;}
	static int startFrames(){return startFrames;}
	
	private static int startLeftOffset=21; //21 works well for k=4
	private static int startRightOffset=10; //10 works well for k=4
	private static int startFrames=startLeftOffset+startRightOffset+1;
	
	private static int stopLeftOffset=12;
	private static int stopRightOffset=15;
	private static int stopFrames=stopLeftOffset+stopRightOffset+1;

	public static boolean PROCESS_PLUS_STRAND=true;
	public static boolean PROCESS_MINUS_STRAND=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         More Statics         ----------------*/
	/*--------------------------------------------------------------*/
	
	//E. coli uses 83% AUG (3542/4284), 14% (612) GUG, 3% (103) UUG[7] and one or two others (e.g., an AUU and possibly a CUG).[8][9]
	public static String[] startCodons=new String[] {"ATG", "GTG", "TTG"};
	public static String[] extendedStartCodons=new String[] {"ATG", "GTG", "TTG", "ATT", "CTG", "ATA"};
	public static String[] stopCodons=new String[] {"TAG", "TAA", "TGA"};
	public static boolean[] isStartCodon=GeneCaller.makeIsCodon(startCodons);
	public static boolean[] isStopCodon=GeneCaller.makeIsCodon(stopCodons);
	public static final boolean isStartCodon(int code){
		return code>=0 && code<=63 && isStartCodon[code];
	}
	public static final boolean isStopCodon(int code){
		return code>=0 && code<=63 && isStopCodon[code];
	}
	
	/*--------------------------------------------------------------*/
	
	private static PrintStream outstream=System.err;
	public static boolean verbose=false;
	public static boolean errorState=false;
	
}
