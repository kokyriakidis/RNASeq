package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

public class PGMTools {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Combines multiple pgm files into a single file */
	public static void main(String[] args){
		
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean overwrite=true;
		boolean allowDupes=false;
		String out=null;
		ArrayList<String> in=new ArrayList<String>(); 
		
		{
			Parser parser=new Parser();
			for(int i=0; i<args.length; i++){
				String arg=args[i];
				String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if(b!=null && b.equalsIgnoreCase("null")){b=null;}
				
				if(a.equals("in")){
					assert(b!=null);
					Tools.addFiles(b, in);
				}else if(b==null && new File(arg).exists()){
					in.add(arg);
				}else if(a.equals("allowdupes") || a.equals("allowduplicates") || a.equals("dupes")){
					allowDupes=Tools.parseBoolean(b);
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
				}

				else if(parser.parse(arg, a, b)){
					//do nothing
				}else{
					outstream.println("Unknown parameter "+args[i]);
					assert(false) : "Unknown parameter "+args[i];
					//				throw new RuntimeException("Unknown parameter "+args[i]);
				}
			}
			overwrite=parser.overwrite;
			out=parser.out1;
		}
		
		checkFileExistence(in, out, overwrite, allowDupes);
		
		ArrayList<GeneModel> models=loadModels(in);
		GeneModel gm=mergeModels(models);
		boolean errorState=writeModel(gm, out, overwrite);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	public static ArrayList<GeneModel> loadModels(ArrayList<String> fnames){
		ArrayList<GeneModel> models=new ArrayList<GeneModel>(fnames.size());
		for(String s : fnames){
			GeneModel pgm=GeneModel.loadModel(s);
			models.add(pgm);
		}
		return models;
	}
	
	public static GeneModel mergeModels(ArrayList<GeneModel> models){
		GeneModel pgmSum=new GeneModel();
		for(GeneModel pgm : models){
			pgmSum.add(pgm);
		}
		return pgmSum;
	}
	
	public static GeneModel loadAndMerge(ArrayList<String> in) {
		ArrayList<GeneModel> models=loadModels(in);
		return mergeModels(models);
	}
	
	public static boolean writeModel(GeneModel pgm, String out, boolean overwrite){
		FileFormat ffout=FileFormat.testOutput(out, FileFormat.PGM, null, true, overwrite, false, false);
		return writeModel(pgm, ffout);
	}
	
	public static boolean writeModel(GeneModel pgm, FileFormat ffout){
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);

		ByteBuilder bb=new ByteBuilder();
		pgm.appendTo(bb);

		boolean errorState=false;
		if(bsw!=null){
			bsw.addJob(bb);
			errorState|=bsw.poisonAndWait();
		}
		return errorState;
	}
	
	/** Ensure files can be read and written */
	private static void checkFileExistence(ArrayList<String> in, String out, boolean overwrite, boolean allowDupes){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, false, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		//Ensure input files can be read
		ArrayList<String> foo=new ArrayList<String>();
		foo.addAll(in);
		if(!Tools.testInputFiles(allowDupes, true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!allowDupes){
			foo.add(out);
			if(!Tools.testForDuplicateFiles(true, foo.toArray(new String[0]))){
				throw new RuntimeException("\nSome file names were specified multiple times.\n");
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private static PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
