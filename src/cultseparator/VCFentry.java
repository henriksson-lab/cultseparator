package cultseparator;

public class VCFentry {

	String ref;
	String alt;
	//store quality?
	
	
	@Override
	public String toString() {
		return ref+">"+alt;
	}
}
