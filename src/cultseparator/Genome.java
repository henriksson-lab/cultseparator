package cultseparator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

public class Genome {
	
	
	TreeMap<String, String> seqs=new TreeMap<String, String>();
	
	public Genome(File f) throws IOException {
		System.out.println("Reading fasta "+f);
		String seqname=null;
		StringBuilder sb=null;
		
		GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(f));
		BufferedReader br = new BufferedReader(new InputStreamReader(gzip));
		
		String line=null;
		while((line=br.readLine())!=null) {
			if(line.startsWith(">")) {
				/*
				if(seqname!=null)
					break;  /// For now
				*/
				
				if(sb!=null) {
					System.out.println("add chrom: "+seqname+"  len "+sb.length());
					seqs.put(seqname,sb.toString());
				}
				seqname=line.substring(1);
				StringTokenizer stok=new StringTokenizer(seqname," ");
				seqname=stok.nextToken();
				sb=new StringBuilder();
			} else {
				sb.append(line);
			}
		}
		
		if(seqname!=null) {
			//StringTokenizer stok=new StringTokenizer(seqname," ");
			//seqname=stok.nextToken();
			System.out.println("add chrom: "+seqname+"  len "+sb.length());
			seqs.put(seqname,sb.toString());
		}
		br.close();
	}

	public String getSeq(String chrom, int start, int len) {
		//System.out.println("chrom "+chrom);
		return seqs.get(chrom).substring(start, start+len);
	}
	

}
