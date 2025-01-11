from mpi4py import MPI
import os

dataset = ["PA220KH-lib09-P19-Tumor_S2","PA221MH-lib09-P19-Norm_S1"]
InputFolder = "Input/" ##Input folder contains .fastq.gz files and Reference folder
OutputFolder = "Results/"

#Creating directory
os.system("mkdir "+OutputFolder+"/1_FastQC")
os.system("mkdir "+OutputFolder+"/2_BWA-Alignment")
os.system("mkdir "+OutputFolder+"/3_BAMConversions")
os.system("mkdir "+OutputFolder+"/4_CallSomaticMutation")


def task2( int ):

####Step 1 FASTQC###
	os.system("time fastqc "+InputFolder+"/"+dataset[int]+"_L001_R1_001.fastq.gz --outdir="+OutputFolder+"/FastQC &> "+OutputFolder+"/FastQC/"+dataset[int]+"_R1.log")
	os.system("time fastqc "+InputFolder+"/"+dataset[int]+"_L001_R2_001.fastq.gz --outdir="+OutputFolder+"/FastQC &> "+OutputFolder+"/FastQC/"+dataset[int]+"_R2.log")


##Step 2 BWA-MEM alignment: In this step, mapping the fastq reads to the current hg38 human reference build to produce a file in SAM (Sequence Alignment mapping) format###	
	os.system("time bwa index "+InputFolder+"/Reference/PANCAN_PDAC_100plex_ref.fa")
	os.system("time samtools faidx "+InputFolder+"/Reference/PANCAN_PDAC_100plex_ref.fa")
 
	os.system("time bwa mem -R '@RG\\tID:"+dataset[int]+"\\tSM:"+dataset[int]+"'  -t 12 "+InputFolder+"/Reference/PANCAN_PDAC_100plex_ref.fa "+InputFolder+"/"+dataset[int]+"_L001_R1_001.fastq.gz "+InputFolder+"/"+dataset[int]+"_L001_R2_001.fastq.gz > "+OutputFolder+"/2_BWA-Alignment/"+dataset[int]+"_alignment.sam")
 
###Step 3 Generating BAM (binary alignment mapping) file and then sort the BAM file based on the coordinate order ###

	os.system("time samtools view -@ 8 -S -b "+OutputFolder+"/2_BWA-Alignment/"+dataset[int]+"_alignment.sam > "+OutputFolder+"/3_BAMConversions/"+dataset[int]+"_aligned.bam")
 
	os.system("time samtools sort -@ 8 -O BAM -o "+OutputFolder+"/3_BAMConversions/"+dataset[int]+"_aligned_sorted.bam "+OutputFolder+"/3_BAMConversions/"+dataset[int]+"_aligned.bam")

	os.system("time samtools index "+OutputFolder+"/3_BAMConversions/"+dataset[int]+"_aligned_sorted.bam")

###Step 4 Somatic genetic variants call using gatk Mutect2### 

	os.system("time gatk CreateSequenceDictionary -R "+InputFolder+"/Reference/PANCAN_PDAC_100plex_ref.fa -O "+InputFolder+"/Reference/PANCAN_PDAC_100plex_ref.dict")

	os.system("time gatk Mutect2 -R "+InputFolder+"/Reference/PANCAN_PDAC_100plex_ref.fa -I "+OutputFolder+"/3_BAMConversions/PA221MH-lib09-P19-Norm_S1_aligned_sorted.bam -I "+OutputFolder+"/3_BAMConversions/PA220KH-lib09-P19-Tumor_S2_aligned_sorted.bam  -normal PA221MH-lib09-P19-Norm_S1 -tumor PA220KH-lib09-P19-Tumor_S2 -A Coverage -A MappingQuality -A StrandBiasBySample -A FisherStrand -A ExcessHet -O "+OutputFolder+"/4_CallSomaticMutation/PA22_Mutect2.vcf.gz")

##Step 5 Extract reuqired Informations == Run this on the terminal
	#bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/FS\t%INFO/MBQ\t%INFO/POPAF\t%INFO/TLOD[\t%GT][\t%AD][\t%DP][\t%FAD][\t%SB]\n' "+OutputFolder+"/4_CallSomaticMutation/PA22_Mutect2.vcf.gz > Varaint_bcftools_Query.csv
 
 

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank < 2:
	task2(rank)
