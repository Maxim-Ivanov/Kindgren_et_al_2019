##### REMAPPING GRO-SEQ DATA #####

# 1) Liu et al., 2018 (PMID 29379150; GSE100010): 
# The most PolII-like track should be from nrpd1/e1 double mutant (because they lack PolIV and PolV);
# TruSeq Small RNA lib prep kit was used;
# FASTQ files are paired-end. Only R1 was remapped;
# 2) Zhu et al., 2018 (PMID 30374093; GSE109974 and GSE117014):
# A custom lib prep protocol from Wang et al., 2011 (PMID 21572438) was used;
# It gives the same flanking sequences as the Illumina Small RNA-Seq kit;
# FASTQ files are single-end.

# In both cases, the strand orientation of reads was not flipped;
# The 5'-terminal base of R1 was taken as proxy for the RNAPII position;
# Unlike plaNET-Seq and pNET-Seq, the GRO-Seq datasets do not contain UMIs. Thus, deduplication step was omitted


# Download paired-end FASTQ files from Liu et al., 2018:
echo -e "SRR5681055\tGROseq_Liu2018_rep1\nSRR5681056\tGROseq_Liu2018_rep2\n" > GROseq_Liu2018.txt

while IFS="\t" read -r line || [[ -n "$line" ]]; do acc="$(cut -f1 <<< $line)" && fname="$(cut -f2 <<< $line)" && echo $acc "->" $fname && fastq-dump --gzip --split-files $acc && mv ${acc}_1.fastq.gz ${fname}_R1.fq.gz && mv ${acc}_2.fastq.gz ${fname}_R2.fq.gz; done < GROseq_Liu2018.txt

# Remove R2 files and rename R1 files:
rm GROseq_Liu2018*_R2.fq.gz
for file in GROseq_Liu2018*_R1.fq.gz; do mv $file ${file/_R1/}; done

# Download single-end FASTQ files from Zhu et al., 2018:
echo -e "SRR6661079\tGROseq_Zhu2018_rep1\nSRR6661080\tGROseq_Zhu2018_rep2\n" > GROseq_Zhu2018.txt

while IFS="\t" read -r line || [[ -n "$line" ]]; do acc="$(cut -f1 <<< $line)" && fname="$(cut -f2 <<< $line)" && echo $acc "->" $fname && fastq-dump --gzip $acc && mv ${acc}.fastq.gz ${fname}.fq.gz; done < GROseq_Zhu2018.txt

# Align both Liu2018 and Zhu2018 data to TAIR10 in single-end mode:
for file in GROseq*fq.gz; do echo $file && STAR --genomeDir tair10 --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/.fq.gz/_} --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 --readFilesCommand zcat --clip3pAdapterSeq TGGAATTCTCGG --outSAMtype BAM Unsorted --sjdbGTFfile Plantsmart28.gff --sjdbGTFtagExonParentTranscript Parent --outSAMstrandField intronMotif; done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *out *tab *STARtmp *STARgenome

# Remove intersections with unwanted genes:
for file in GROseq*bam; do echo $file && bedtools intersect -v -abam $file -b Unwanted.bed > ${file/.bam/_clean.bam}; done

# Sort and filter by MAPQ values:
for file in GROseq*clean.bam; do echo $file && samtools view -hu -q 10 $file | samtools sort - -o ${file/.bam/_mapq.bam}; done

# Merge the replicates:
for f1 in GROseq*rep1*mapq.bam; do f2=${file/rep1/rep2} && echo $f1 $f2 && samtools merge ${f1/rep1/merged/} $f1 $f2; done

# Make strand-specific Bedgraph files (without strand switch; coverage of 5' bases only):
for str in "+" "-"; do [ "$str" = "+" ] && n="fw" || n="rev"; for file in GROseq*mapq.bam; do sample=${file/_clean_mapq.bam/} && echo $n $sample && bedtools genomecov -ibam $file -bg -5 -strand $str > ${sample}_${n}.bg; done; done

# Merge forward and reverse Bedgraph files for the same sample:
# (sense and antisense strand coverages are encoded by positive and negative values in column 4)
f_str="fw"; r_str="rev"; ext=".bg"; for file1 in GROseq*${f_str}${ext}; do file2=${file1/${f_str}/${r_str}} && outfile=${file1/${f_str}${ext}/fw_rev.bedgraph} && echo $file1 "+" $file2 "=" $outfile && awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | cat $file1 - | sort -k1,1 -k2,2n | sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' > $outfile; done

# Compress Bedgraph files:
for file in *bedgraph; do echo $file && gzip $file; done


##### REMAPPING STRAND-SPECIFIC RNA-SEQ #####

# Kohnen 2016 (PMID 27923878; GSE81202)
# The best proxy for Col-0 seedlings are samples GSM2144569 and GSM2144571 (TP0 white light cotyledon);

# Download single-end FASTQ files:
echo -e "SRR3480142\tRNAseq_Kohnen2016_rep1\nSRR3480144\tRNAseq_Kohnen2016_rep2\n" > RNAseq_Kohnen2016.txt

while IFS="\t" read -r line || [[ -n "$line" ]]; do acc="$(cut -f1 <<< $line)" && fname="$(cut -f2 <<< $line)" && echo $acc "->" $fname && fastq-dump --gzip $acc && mv ${acc}.fastq.gz ${fname}.fq.gz; done < RNAseq_Kohnen2016.txt

# Align to TAIR10 (Extend5pOfRead1 mode, TxDb-guided):
for file in RNAseq*fq.gz; do echo $file && STAR --genomeDir tair10 --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/.fq.gz/_} --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 --readFilesCommand zcat --clip3pAdapterSeq AGATCGGAAGAGC --outSAMtype BAM Unsorted --sjdbGTFfile ${basedir}/ann/Plantsmart28.gff --sjdbGTFtagExonParentTranscript Parent --outSAMstrandField intronMotif; done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARgenome *STARtmp *out *tab

# Sort BAM files, deduplicate on start coordinates, remove intersections with unwanted genes and low MAPQ reads:
for file in RNAseq*bam; do echo $file && samtools view -hu -q 10 | samtools sort - -o - | samtools rmdup -s - - | bedtools intersect -v -abam stdin -b Unwanted.bed > ${file/.bam/_clean.bam}; done

# Merge the replicates:
samtools merge RNAseq_Kohnen2016_merged_clean.bam RNAseq*rep?_clean.bam

# Make strand-specific Bedgraph files (with strand switch; coverage of the whole reads):
for str in "+" "-"; do [ "$str" = "+" ] && n="rev" || n="fw"; for file in RNAseq*clean.bam; do sample=${file/_clean.bam/} && echo $n $sample && bedtools genomecov -ibam $file -bg -split -strand $str > ${sample}_${n}.bg; done; done

# Merge forward and reverse Bedgraph files:
f_str="fw"; r_str="rev"; ext=".bg"; for file1 in RNAseq*${f_str}${ext}; do file2=${file1/${f_str}/${r_str}} && outfile=${file1/${f_str}${ext}/fw_rev.bedgraph.gz} && echo $file1 "+" $file2 "=" $outfile && awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | cat $file1 - | sort -k1,1 -k2,2n | sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | gzip > $outfile; done


##### REMAPPING DIRECT RNA-SEQ (DR-SEQ) ######
# 1) Schurch et al., 2014 (PMID 24722185): ERP003245 contains Helicos FASTQ files;
# 2) Sherstnev et al., 2012 (PMID 22820990): ERP001018 contains both Helicos FASTQ and BAM files;

# Download replicates of the wild type from Schurch et al., 2014:
base="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/ERA223/ERA223202"
wget ${base}/ERX267337/ERR294004.fastq.bz2 DRseq_Schurch2014_rep1.fq.bz2
wget ${base}/ERX267338/ERR294005.fastq.bz2 DRseq_Schurch2014_rep2.fq.bz2
wget ${base}/ERX267339/ERR294006.fastq.bz2 DRseq_Schurch2014_rep3.fq.bz2

# Remove trailing white spaces after read names:
for file in *Schurch2014*fq.bz2; do echo $file && bzcat $file | sed 's/ $//' | gzip > ${file/bz2/gz}; done

# Align Schurch 2014 data to TAIR10 by STAR (Local mode, not transcriptome-guided):
for file in *Schurch2014*fq.gz; do echo $file && STAR --genomeDir tair10 --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/.fq.gz/_} --outSAMmultNmax 1 --alignEndsType Local --readFilesCommand zcat --outSAMtype BAM Unsorted; done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rmdir *STARtmp; rm *out *tab

# Sort BAM files and remove low MAPQ reads:
for file in *Schurch2014*bam; do echo $file && samtools view -hu -q 10 $file | samtools sort - -o ${file/.bam/_mapq.bam}; done

# Merge all replicates:
samtools merge DRseq_Schurch2014_merged_mapq.bam *Schurch2014*rep?_mapq.bam

# Make strand-specific Bedgraph files (use only the first base; switch the strand orientation):
for str in "+" "-"; do [ "$str" = "+" ] && n="rev" || n="fw"; for file in *Schurch2014*mapq.bam; do sample=${file/_mapq.bam/} && echo $n $sample && bedtools genomecov -ibam $file -bg -5 -strand $str | sort -k1,1 -k2,2n > ${sample}_${n}.bg; done; done

# Download HeliScope BAM files from Sherstnev 2012:
base2="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR055"
base3="DRseq_Sherstnev2012_"
wget ${base2}/ERR055397/wt.t1.rep1.corr.bam ${base3}wt1_rep1.bam
wget ${base2}/ERR055398/wt.t1.rep2.corr.bam ${base3}wt1_rep2.bam
wget ${base2}/ERR055399/wt.t1.rep3.corr.bam ${base3}wt1_rep3.bam
wget ${base2}/ERR055400/wt.t2.rep1.corr.bam ${base3}wt2_rep1.bam
wget ${base2}/ERR055401/wt.t2.rep2.corr.bam ${base3}wt2_rep2.bam

# Adjust chromosome names in "wt2_rep1" file:
file=${base3}"wt2_rep1.bam"
samtools view -H $file | sed 's/chr/Chr/;s/ChrMt/mitochondria/;s/ChrPt/chloroplast/' > header.sam
samtools reheader header.sam $file > temp && mv temp $file

# Merge all samples from Sherstnev 2012:
file2=${base3}"merged.bam"
samtools merge $file2 ${base3}*rep?.bam

# Remove the individual replicates:
rm ${base3}*rep?.bam

# Fix the chromosome names:
samtools view -H $file2 | sed 's/Chr//;s/mitochondria/Mt/;s/chloroplast/Pt/' > header.sam
samtools reheader header.sam $file2 > temp && mv temp $file2

# Make strand-specific Bedgraph files (use only the first base; no strand switch!):
for str in "+" "-"; do [ "$str" = "-" ] && n="rev" || n="fw" && bedtools genomecov -ibam $file2 -bg -5 -strand $str | sort -k1,1 -k2,2n > ${file2/.bam/}_${n}.bg; done

# Merge forward and reverse Bedgraph files for both DR-Seq studies:
f_str="fw"; r_str="rev"; ext=".bg"; for file1 in DRseq*${f_str}${ext}; do file2=${file1/${f_str}/${r_str}} && outfile=${file1/${f_str}${ext}/fw_rev.bedgraph.gz} && echo $file1 "+" $file2 "=" $outfile && awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | cat $file1 - | sort -k1,1 -k2,2n | sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | gzip > $outfile; done

