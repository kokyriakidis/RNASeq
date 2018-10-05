package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import dna.AminoAcid;
import dna.Data;
import dna.Gene;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.PreParser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * @author bbushnell
 * @date Sep 24, 2018
 *
 */
public class CallGenes {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		CallGenes x=new CallGenes(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CallGenes(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;

			outGff=parser.out1;
			maxReads=parser.maxReads;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffoutGff=FileFormat.testOutput(outGff, FileFormat.GFF, null, true, overwrite, append, ordered);
		ffoutAmino=FileFormat.testOutput(outAmino, FileFormat.FA, null, true, overwrite, append, ordered);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

//			outstream.println(arg+", "+a+", "+b);
			if(a.equals("in") || a.equals("infna") || a.equals("fnain") || a.equals("fna")){
				assert(b!=null);
				Tools.addFiles(b, fnaList);
			}else if(a.equals("in") || a.equals("infna") || a.equals("fnain") || a.equals("fna")){
				assert(b!=null);
				Tools.addFiles(b, fnaList);
			}else if(b==null && new File(arg).exists() && FileFormat.isFastaFile(arg)){
				fnaList.add(arg);
			}else if(a.equals("pgm") || a.equals("gm") || a.equals("model")){
				assert(b!=null);
				if(b.equalsIgnoreCase("auto") || b.equalsIgnoreCase("default")){
					b=Data.findPath("?model.pgm");
					pgmList.add(b);
				}else{
					Tools.addFiles(b, pgmList);
				}
			}else if(b==null && new File(arg).exists() && FileFormat.isPgmFile(arg)){
				Tools.addFiles(arg, pgmList);
			}else if(a.equals("outamino") || a.equals("aminoout") || a.equals("outa") || a.equals("outaa") || a.equals("amino")){
				outAmino=b;
			}else if(a.equals("k")){
				int k=Integer.parseInt(b);
				GeneModel.setInnerK(k);
				GeneModel.setEndK(k);
			}else if(a.equals("innerk")){
				int k=Integer.parseInt(b);
				GeneModel.setInnerK(k);
			}else if(a.equals("endk")){
				int k=Integer.parseInt(b);
				GeneModel.setEndK(k);
			}else if(a.equals("startleftoffset")){
				int x=Integer.parseInt(b);
				GeneModel.setStartLeftOffset(x);
			}else if(a.equals("startrightoffset")){
				int x=Integer.parseInt(b);
				GeneModel.setStartRightOffset(x);
			}else if(a.equals("stopleftoffset")){
				int x=Integer.parseInt(b);
				GeneModel.setStopLeftOffset(x);
			}else if(a.equals("stoprightoffset")){
				int x=Integer.parseInt(b);
				GeneModel.setStopRightOffset(x);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(a.equals("plus")){
				GeneModel.PROCESS_PLUS_STRAND=Tools.parseBoolean(b);
			}else if(a.equals("minus")){
				GeneModel.PROCESS_MINUS_STRAND=Tools.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Tools.parseBoolean(b);
			}
			
			else if(a.equals("minlen") || a.equals("minlength")){
				minLen=Integer.parseInt(b);
			}else if(a.equals("maxoverlapss") || a.equals("overlapss") || a.equals("overlapsamestrand") || a.equals("moss")){
				maxOverlapSameStrand=Integer.parseInt(b);
			}else if(a.equals("maxoverlapos") || a.equals("overlapos") || a.equals("overlapoppositestrand") || a.equals("moos")){
				maxOverlapOppositeStrand=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minStartScore")){
				minStartScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minStopScore")){
				minStopScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minInnerScore") || a.equalsIgnoreCase("minKmerScore")){
				minKmerScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minOrfScore") || a.equalsIgnoreCase("minScore")){
				minOrfScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("breakLimit")){
				GeneCaller.breakLimit=Integer.parseInt(b);
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		if(pgmList.isEmpty()){
			String b=Data.findPath("?model.pgm");
			pgmList.add(b);
		}
		
		if(Shared.threads()<2){ordered=false;}
		assert(!fnaList.isEmpty()) : "At least 1 fasta file is required.";
		assert(!pgmList.isEmpty()) : "At least 1 pgm file is required.";
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		fnaList=Tools.fixExtension(fnaList);
		pgmList=Tools.fixExtension(pgmList);
		if(fnaList.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, outGff, outAmino)){
			outstream.println((outGff==null)+", "+outGff);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+outGff+", "+outAmino+"\n");
		}
		
		//Ensure input files can be read
		ArrayList<String> foo=new ArrayList<String>();
		foo.addAll(fnaList);
		foo.addAll(pgmList);
		if(!Tools.testInputFiles(false, true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		foo.add(outGff);
		foo.add(outAmino);
		if(!Tools.testForDuplicateFiles(true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Create read streams and process all data */
	void process(Timer t){
		
		GeneModel pgm=PGMTools.loadAndMerge(pgmList);
		
		ByteStreamWriter bsw=makeBSW(ffoutGff);
		if(bsw!=null){
			bsw.forcePrint("##gff-version 3\n");
		}
		
		ConcurrentReadOutputStream ros=makeCros(ffoutAmino);
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Reset counters
		readsIn=genesOut=0;
		basesIn=bytesOut=0;
		
		for(String fname : fnaList){
			//Create a read input stream
			final ConcurrentReadInputStream cris=makeCris(fname);

			//Process the reads in separate threads
			spawnThreads(cris, bsw, ros, pgm);
			
			//Close the input stream
			errorState|=ReadWrite.closeStreams(cris, ros);
		}
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the output stream
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		outstream.println(Tools.linesBytesOut(readsIn, basesIn, genesOut, bytesOut, 8, false));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private ConcurrentReadInputStream makeCris(String fname){
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ffin, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ByteStreamWriter bsw, ConcurrentReadOutputStream ros, GeneModel pgm){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, bsw, ros, pgm, minLen, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for threads to finish
		waitForThreads(alpt);
		
		//Do anything necessary after processing
		
	}
	
	private void waitForThreads(ArrayList<ProcessThread> alpt){
		
		//Wait for completion of all threads
		boolean success=true;
		for(ProcessThread pt : alpt){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			genesOut+=pt.genesOutT;
			bytesOut+=pt.bytesOutT;
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	private ConcurrentReadOutputStream makeCros(FileFormat ff){
		if(ff==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(4, 64, (Shared.threads()*2)/3) : 4);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ByteStreamWriter bsw_, ConcurrentReadOutputStream ros_,
				GeneModel pgm_, final int minLen, final int tid_){
			cris=cris_;
			bsw=bsw_;
			ros=ros_;
			pgm=pgm_;
			tid=tid_;
			caller=new GeneCaller(minLen, maxOverlapSameStrand, maxOverlapOppositeStrand, minStartScore, minStopScore, minKmerScore, minOrfScore, pgm);
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processList(ln);

				//Fetch a new list
				ln=cris.nextList();
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		void processList(ListNum<Read> ln){
			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;

//			System.err.println(reads.size());
			
			ArrayList<Orf> orfList=new ArrayList<Orf>();
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				final Read r2=r1.mate;
				
				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}
				if(r2!=null && !r2.validated()){r2.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r1.length();
				final int initialLength2=r1.mateLength();

				//Increment counters
				readsInT+=r1.pairCount();
				basesInT+=initialLength1+initialLength2;
				
				{
					//Reads are processed in this block.
					{
						ArrayList<Orf> list=processRead(r1);
						if(list!=null){orfList.addAll(list);}
					}
					if(r2!=null){
						ArrayList<Orf> list=processRead(r2);
						if(list!=null){orfList.addAll(list);}
					}
				}
			}

			genesOutT+=orfList.size();
			ByteBuilder bb=new ByteBuilder();
			
			if(bsw!=null){
				if(bsw.ordered){
					for(Orf orf : orfList){
						orf.appendGff(bb);
						bb.nl();
					}
					bsw.add(bb, ln.id);
					bytesOutT+=bb.length();
				}else{
					for(Orf orf : orfList){
						orf.appendGff(bb);
						bb.nl();
						if(bb.length()>=32000){
							bsw.print(bb);
							bytesOutT+=bb.length();
							bb.clear();
						}
					}
					if(bb.length()>0){
						bsw.addJob(bb);
						bytesOutT+=bb.length();
					}
				}
			}

			//Notify the input stream that the list was used
			cris.returnList(ln);
//			if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		ArrayList<Orf> processRead(final Read r){
			ArrayList<Orf> list=caller.callGenes(r, pgm);
			
			if(ros!=null && list!=null && !list.isEmpty()){
				ArrayList<Read> prots=new ArrayList<Read>();
				for(int strand=0; strand<2; strand++){
					for(Orf orf : list){
						if(orf.strand==strand){
							Read aa=translate(orf, r.bases, r.id);
							prots.add(aa);
						}
					}
					r.reverseComplement();
				}
				ros.add(prots, r.numericID);
			}
			
			return list;
		}
		
		private Read translate(Orf orf, byte[] bases, String id){
			if(orf.strand==1){orf.flip();}
			byte[] acids=AminoAcid.toAAs(bases, orf.start, orf.stop);
			if(orf.strand==1){orf.flip();}
			Read r=new Read(acids, null, id+"\t"+(Gene.strandCodes[orf.strand]+"\t"+orf.start+"-"+orf.stop), 0, Read.AAMASK);
			return r;
		}

		/** Number of reads processed by this thread */
		protected long readsInT=0;
		/** Number of bases processed by this thread */
		protected long basesInT=0;
		
		/** Number of genes called by this thread */
		protected long genesOutT=0;
		/** Number of bytes written by this thread */
		protected long bytesOutT=0;
		
		protected ConcurrentReadOutputStream ros;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ByteStreamWriter bsw;
		/** Gene Model for annotation (not really needed) */
		private final GeneModel pgm;
		/** Gene Caller for annotation */
		private final GeneCaller caller;
		/** Thread ID */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	private long readsIn=0;
	private long basesIn=0;
	private long genesOut=0;
	private long bytesOut=0;
	
	//These call 98.8% of E.coli gene stops (according to Prodigal) with a 5.7% FP rate.
	private int minLen=60;
	private int maxOverlapSameStrand=80;
	private int maxOverlapOppositeStrand=120;
	
	private float minStartScore=-0.04f;
	private float minStopScore=-0.04f;
	private float minKmerScore=0.09f;
	private float minOrfScore=30f; //Higher increases SNR dramatically but reduces TP rate
	
	/*--------------------------------------------------------------*/

	private ArrayList<String> fnaList=new ArrayList<String>();
	private ArrayList<String> pgmList=new ArrayList<String>();
	private String outGff=null;
	private String outAmino=null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffoutGff;
	private final FileFormat ffoutAmino;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	private boolean ordered=true;
	
}
