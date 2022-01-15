package cultseparator;

import java.io.File;
import java.io.IOException;
import java.util.Map.Entry;
import java.util.SortedMap;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * 
 * Takes a BAM file and splits it into 2, based on VCF files
 * 
 * @author Johan Henriksson
 *
 */
public class Main {
	
	/**
	 * Entry point
	 * 
	 * 
	 * how to make a subset: samtools view -h big.bam | head -n 1000 | samtools view -bS - > little.bam
	 * 
	 */
	public static void main(String[] args) throws IOException {
		
		File fVCFa=new File("data/pc3_mono.vcf");
		File fVCFb=new File("data/U937_mono.vcf");
		File fFASTA=new File("data/all.fa.gz");
		
		File fIn=new File("data/little_u937.bam");
		File fOutA=new File("data/outA.bam");
		File fOutB=new File("data/outB.bam");

		if(args.length>0) {
			if(args.length==6) {
				
				fVCFa=new File(args[0]);
				fVCFb=new File(args[1]);
				fFASTA=new File(args[2]);
				
				fIn=new File(args[3]);
				fOutA=new File(args[4]);
				fOutB=new File(args[5]);
			} else {
				System.out.println("1.vcf 2.vcf reference.fa.gz input.bam out1.bam out2.bam");
				System.exit(0);
			}
		}
		
		System.out.println("Reading VCFs");
		VCFfile refA=new VCFfile(fVCFa);
		VCFfile refB=new VCFfile(fVCFb);
		refA.subtract(refB);
		refB.subtract(refA);
		refA.createOpposite(refB);
		refB.createOpposite(refA);
		
		System.out.println("read genome");
		Genome gen=new Genome(fFASTA);
		
		System.out.println("splitting bam");
		
		splitBAM(gen, refA, refB, 
				fIn,
				fOutA,fOutB);
		System.out.println("Done");
	}
	
	
	public static void splitBAM(
			Genome gen,
			VCFfile refA, 
			VCFfile refB, 
			File fBAM,
			File fOutA,
			File fOutB
		) throws IOException {
			//Open BAM file
			final SamReader reader = SamReaderFactory.makeDefault().open(fBAM);

			//Create writers
			boolean isPresorted=reader.getFileHeader().getSortOrder()!=SAMFileHeader.SortOrder.coordinate;
			if(fOutA.exists())
				fOutA.delete();
			if(fOutB.exists())
				fOutB.delete();
			SAMFileWriter writerA=new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), isPresorted, fOutA);
			SAMFileWriter writerB=new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), isPresorted, fOutB);
			
			//Loop through all SAM records
			int readRecords=0;
			int ambigous=0;
			for (final SAMRecord samRecord : reader) {
				
				//Update user about progress
				readRecords++;
				if(readRecords%1000 == 0){
					//Calculate progress
					System.out.println("records so far: "+readRecords);
				}

				String chrom=samRecord.getContig();
				if(!refA.entries.containsKey(chrom)) {
					System.out.println("chrom not in vcf, "+chrom);
					continue;
				}
				
				
				String readseq=samRecord.getReadString();

				//example cigar: 77M140N69M757N4M
				//D and N are the same. but N is large
				//for simplicity, 
				//M is both matches and mismatches. alignment blocks are over M
				
				int matchesA=0;
				int matchesB=0;
				int sizeEntriesA=0;
				int sizeEntriesB=0;
				for(AlignmentBlock ab:samRecord.getAlignmentBlocks()) {
					int start=ab.getReadStart()-1;   //only -1 on the first round???
					int len=ab.getLength();
					int refstart=ab.getReferenceStart();
					String readSeq=readseq.substring(start, start+len-1); 
					
					//Fetch diffs
					SortedMap<Integer,VCFentry> entriesA=refA.entries.get(chrom);
					SortedMap<Integer,VCFentry> entriesB=refB.entries.get(chrom);
					sizeEntriesA+=entriesA.size();
					sizeEntriesB+=entriesB.size();
					
					//limits here are inclusive. Note that VCF starts from 1!
					entriesA=entriesA.subMap(refstart, refstart+len-1);  
					entriesB=entriesB.subMap(refstart, refstart+len-1);  
					
					//Are there any differences?
					for(Entry<Integer,VCFentry> e:entriesA.entrySet()) {
						int j=e.getKey()-refstart;
						char readchar=readSeq.charAt(j);
						if(e.getValue().alt.equals(Character.toString(readchar))) {
							matchesA++;
						} 
					}
					for(Entry<Integer,VCFentry> e:entriesB.entrySet()) {
						int j=e.getKey()-refstart;
						char readchar=readSeq.charAt(j);
						if(e.getValue().alt.equals(Character.toString(readchar))) {
							matchesB++;
						} 
					}						
				}
				
				
				/*
				System.out.println("########### "+
						matchesA+"/"+sizeEntriesA+"  "+
						matchesB+"/"+sizeEntriesB);			
				*/

				//Only bother if there are any matches
				if(matchesA>0 || matchesB>0) {
					
					if(sizeEntriesA>0 & sizeEntriesB>0) {
						//Some sort of numeric comparison. Ignores total number of VCFs. Fair?
						if(matchesA > matchesB) {
							writerA.addAlignment(samRecord);
						} else if(matchesA < matchesB) {
							writerB.addAlignment(samRecord);
						} else {
							ambigous++;
						}
					} else {
						System.out.println("this should not happen");
					}
				} else {
					ambigous++;
				}
			}
			writerA.close();
			writerB.close();

			System.out.println("total records "+readRecords);
			System.out.println("ambious: "+ambigous);
	}
}
