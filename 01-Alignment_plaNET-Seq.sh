# Trim UMIs from both R1 and R2:
for f1 in *R1.fq.gz; do f2=${f1/R1/R2} && echo $f1 $f2 && umi_tools extract --stdin=${f1} --read2-in=${f2} --bc-pattern=NNNN --bc-pattern2=NNNN --stdout=${f1/.fq.gz/_UMI.fq.gz} --read2-out=${f2/.fq.gz/_UMI.fq.gz}; done

# (Optional) Export transcriptome annotation from Bioconductor:
# library(rtracklayer)
# library(TxDb.Athaliana.BioMart.plantsmart28)
# txdb <- TxDb.Athaliana.BioMart.plantsmart28
# export.gff3(txdb, "Plantsmart28.gff")

# Align to TAIR10 in SE mode (use only R2):
for file in *R2_UMI.fq.gz; do echo $file && STAR --genomeDir tair10 --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/_R2_UMI.fq.gz/_} --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 --readFilesCommand zcat --clip3pAdapterSeq GATCGTCGGACT --outSAMtype BAM Unsorted --sjdbGTFfile Plantsmart28.gff --sjdbGTFtagExonParentTranscript Parent --outSAMstrandField intronMotif; done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARgenome *STARtmp *out *tab

# Sort BAM files:
for file in *bam; do echo $file && samtools sort $file -o ${file/.bam/_sorted.bam} && rm $file; done

# Deduplicate (UMI-Tools):
for file in *sorted.bam; do echo $file && samtools index $file && umi_tools dedup --stdin=${file} --stdout=${file/.bam/_dedup.bam}; done

# (Optional) Extract coordinates of non-RNAPII-transcribed genes (rRNA, tRNA, snRNA, snoRNA) from Araport11:
# Download Araport11_GFF3_genes_transposons.201606.gff.gz from www.arabidopsis.org;
# zcat Araport11_GFF3_genes_transposons.201606.gff.gz | sed '/^#/d;s/^Chr//;s/^C/Pt/;s/^M/Mt/' | awk 'BEGIN{OFS="\t"}{if ($3=="rRNA" || $3=="tRNA" || $3=="snRNA" || $3=="snoRNA") print $1,$4-100,$5+100,$7}' | sort -k1,1 -k2,2n > Unwanted.bed

# Remove intersections with unwanted genes:
for file in *dedup.bam; do echo $file && bedtools intersect -v -abam $file -b Unwanted.bed > ${file/.bam/_clean.bam}; done

# Remove low MAPQ reads:
for file in *clean.bam; do echo $file && samtools view -hb -q 10 $file > ${file/.bam/_mapq.bam}; done

# Merge HiSeq and MiSeq BAM files:
for f1 in *HiSeq*mapq.bam; do f2=${file/HiSeq/MiSeq} && echo $f1 $f2 && samtools merge ${f1/HiSeq/full/} $f1 $f2; done
