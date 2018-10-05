package dna;

public class Feature implements Comparable<Feature>{

	public Feature(String scafName_, int start_, int stop_, int strand_){
		scafName=scafName_;
		start=start_;
		stop=stop_;
		strand=strand_;
	}
	
	public void flip(int scafLen){
		int a=scafLen-start-1;
		int b=scafLen-stop-1;
		start=b;
		stop=a;
		flipped=flipped^1;
	}
	
	public int currentStrand(){
		return strand^flipped;
	}
	
	public int length(){
		return stop-start+1;
	}
	
	@Override
	public int compareTo(Feature f) {
		int x=scafName.compareTo(f.scafName);
		if(x!=0){return x;}
		if(stop!=f.stop){return stop-f.stop;}
		return start-f.start;
	}
	
	public int flipped(){return flipped;}
	
	public final String scafName;
	public final int strand;
	
	public int start;
	public int stop;
	private int flipped=0;
	
}
