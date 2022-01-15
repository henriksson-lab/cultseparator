package cultseparator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.Map.Entry;

/**
 * 
 * VCF file reader and some utils to check diffs
 * 
 * @author Johan Henriksson
 *
 */
public class VCFfile {

	HashMap<String, TreeMap<Integer, VCFentry>> entries=new HashMap<String, TreeMap<Integer,VCFentry>>();

	
	public VCFfile(File f) throws IOException {
		System.out.println("Reading vcf "+f);
		BufferedReader br=new BufferedReader(new FileReader(f));
		String line;
		while((line=br.readLine())!=null) {
			if(!line.startsWith("#")) {
				// 1	8894165	.	T	C	43.4147	.	DP=8;VDB=0.00911214;SGB=-0.556411;MQSBZ=0;FS=0;MQ0F=0;AC=2;AN=2;DP4=0,0,2,2;MQ=20	GT:PL	1/1:73,12,0
				StringTokenizer stok=new StringTokenizer(line,"\t");

				String chrom=stok.nextToken();
				int pos=Integer.parseInt(stok.nextToken());
				stok.nextToken();
				
				VCFentry e=new VCFentry();
				e.ref=stok.nextToken();
				e.alt=stok.nextToken();
				
				TreeMap<Integer, VCFentry> onechrom=entries.get(chrom);
				if(onechrom==null)
					entries.put(chrom,onechrom=new TreeMap<Integer, VCFentry>());
				
				//Only keep SNVs
				if(e.ref.length()==1 && e.alt.length()==1)
					onechrom.put(pos,e);
			}
		}
		br.close();
	}
	
	
	/**
	 * Subtract common entries
	 */
	public void subtract(VCFfile other) {
		//For every chromosome
		int startsize=size();
		int removed=0;
		for(String chr:entries.keySet()) {
			TreeMap<Integer, VCFentry> thisE=entries.get(chr);
			TreeMap<Integer, VCFentry> otherE=other.entries.get(chr);
			
			if(otherE!=null) {
				//For every entry
				TreeMap<Integer, VCFentry> loopE=new TreeMap<Integer, VCFentry>(thisE);
				for(Entry<Integer,VCFentry> e:loopE.entrySet()) {
					int pos=e.getKey();
					if(otherE.containsKey(pos)) {
						VCFentry v1=e.getValue();
						VCFentry v2=otherE.get(pos);
						if(v1.alt.equals(v2.alt)) {
							thisE.remove(pos);
							otherE.remove(pos);
							removed++;
						}
					}
				}				
			}
		}
		System.out.println("Removed VCF entries: "+removed+"   out of "+startsize);
		
	}
	
	
	public int size() {
		int size=0;
		for(TreeMap<Integer, VCFentry> chr:entries.values()) {
			size+=chr.size();
		}
		return size;
	}
	
	
	/**
	 * Assume the other VCF has the reference letter if no entry. This is not ideal,
	 * but helps counting. Ideally, should rather use the VCF to go back and calculate the occurence
	 * at interesting sites
	 * 
	 */
	public void createOpposite(VCFfile other) {
		
		//For every chromosome
		int startsize=size();
		int added=0;
		for(String chr:entries.keySet()) {
			TreeMap<Integer, VCFentry> thisE=entries.get(chr);
			TreeMap<Integer, VCFentry> otherE=other.entries.get(chr);
			
			if(otherE==null) {
				other.entries.put(chr,otherE=new TreeMap<Integer, VCFentry>());
			}
			
			//For every entry
			TreeMap<Integer, VCFentry> loopE=new TreeMap<Integer, VCFentry>(thisE);
			for(Entry<Integer,VCFentry> e:loopE.entrySet()) {
				int pos=e.getKey();
				if(!otherE.containsKey(pos)) {
					VCFentry v1=e.getValue();
					VCFentry v2=new VCFentry();
					v2.ref=v1.ref;
					v2.alt=v1.ref;
					otherE.put(pos,v2);
					added++;
				}
			}				
		}
		System.out.println("Added VCF entries: "+added+"   out of "+startsize);
		
	}
	
	
	
	
}
