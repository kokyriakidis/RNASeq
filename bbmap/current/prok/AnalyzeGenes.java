package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/**
 * This class is designed to analyze paired prokaryotic fna and gff files
 * to calculate the patterns in coding and noncoding frames, start and stop sites.
 * @author bbushnell
 * @date Sep 27, 2018
 *
 */
public class AnalyzeGenes {
	
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
		AnalyzeGenes x=new AnalyzeGenes(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public AnalyzeGenes(String[] args){
		
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

			out=parser.out1;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout=FileFormat.testOutput(out, FileFormat.PGM, null, true, overwrite, append, false);
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
			if(a.equals("in") || a.equals("infna") || a.equals("fnain") || a.equals("fna") || a.equals("ref")){
				assert(b!=null);
				Tools.addFiles(b, fnaList);
			}else if(a.equals("gff") || a.equals("ingff") || a.equals("gffin")){
				assert(b!=null);
				Tools.addFiles(b, gffList);
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
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(arg.indexOf('=')<0 && new File(arg).exists() && FileFormat.isFastaFile(arg)){
				fnaList.add(arg);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		if(gffList.isEmpty()){
			for(String s : fnaList){
				String prefix=ReadWrite.stripExtension(s);
				String gff=prefix+".gff";
				File f=new File(gff);
				if(!f.exists()){
					String gz=gff+".gz";
					f=new File(gz);
					assert(f.exists() && f.canRead()) : "Can't read file "+gff;
					gff=gz;
				}
				gffList.add(gff);
			}
		}
		assert(gffList.size()==fnaList.size()) : "Number of fna and gff files do not match: "+fnaList.size()+", "+gffList.size();
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		fnaList=Tools.fixExtension(fnaList);
		gffList=Tools.fixExtension(gffList);
		if(fnaList.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		//Ensure input files can be read
		ArrayList<String> foo=new ArrayList<String>();
		foo.addAll(fnaList);
		foo.addAll(gffList);
		if(!Tools.testInputFiles(false, true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		foo.add(out);
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
	
	void process(Timer t){
		
		GeneModel pgm=makeModel();
		
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);
		
		ByteBuilder bb=new ByteBuilder();
		pgm.appendTo(bb);
		bytesOut+=bb.length;
		
		if(bsw!=null){
			bsw.addJob(bb);
			errorState|=bsw.poisonAndWait();
		}
		
		t.stop();
		
		outstream.println(timeReadsBasesGenesProcessed(t, pgm.readsProcessed, pgm.basesProcessed, pgm.genesProcessed, 8));
		
		//outstream.println("Bytes Out:         \t"+bytesOut);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private static String timeReadsBasesGenesProcessed(Timer t, long readsProcessed, long basesProcessed, long genesProcessed, int pad){
		return ("Time:                         \t"+t+"\n"+readsBasesGenesProcessed(t.elapsed, readsProcessed, basesProcessed, genesProcessed, pad));
	}
	
	private static String readsBasesGenesProcessed(long elapsed, long reads, long bases, long genes, int pad){
		double rpnano=reads/(double)elapsed;
		double bpnano=bases/(double)elapsed;
		double gpnano=genes/(double)elapsed;

		String rstring=Tools.padKM(reads, pad);
		String bstring=Tools.padKM(bases, pad);
		String gstring=Tools.padKM(genes, pad);
		StringBuilder sb=new StringBuilder();
		sb.append("Sequences Processed:").append(rstring).append(String.format(Locale.ROOT, " \t%.2fk seqs/sec", rpnano*1000000)).append('\n');
		sb.append("Bases Processed:    ").append(bstring).append(String.format(Locale.ROOT, " \t%.2fm bases/sec", bpnano*1000)).append('\n');
		sb.append("Genes Processed:    ").append(gstring).append(String.format(Locale.ROOT, " \t%.2fk genes/sec", gpnano*1000000));
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private GeneModel makeModel(){
		GeneModel pgmSum=new GeneModel();
		
		for(int i=0; i<fnaList.size(); i++){
			String fna=fnaList.get(i);
			String gff=gffList.get(i);
			GeneModel pgm=new GeneModel();
			pgm.process(fna, gff);
			pgmSum.add(pgm);
		}
		return pgmSum;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<String> fnaList=new ArrayList<String>();
	private ArrayList<String> gffList=new ArrayList<String>();
	private IntList taxList=new IntList();
	private String out=null;
	
	/*--------------------------------------------------------------*/
	
	private long bytesOut=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}

