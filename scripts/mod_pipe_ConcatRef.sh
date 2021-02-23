#! /bin/bash
#$ -S /bin/bash
#$ -cwd


AlignedPath=$1 #path/to/bam_file/

INPUT=$2 #bam_file.bam #aligned to ConcatRef

tag=$3 #tag="UCSC_mm10"

ID=${INPUT%%.bam*}


# Extract hg38 aligned reads 

#this command gets the input chromosomes from idxstats, cuts out the other columns, removes anything with _NCBI_GRCm38 or * 
chr_list=$(samtools idxstats  $AlignedPath$INPUT | cut -f 1 | grep -v $tag | grep -v \* )

echo $chr_list

#AL note: this next line just takes all reads that aligned with the human genome and puts them in a new file
#however it caused problems because there are pairs where one read is aligned to a normal chromosome and one is aligned to a 'junk chromosome' this lead to only a single read from the pair getting retained
#to fix this, I altered the chrom list
# samtools view -b $AlignedPath$INPUT chrM $chr_list chrX chrY > $AlignedPath$ID.human.bam #AL modified filenames to remove .hg19 and add .human_unfiltered
samtools view -b $AlignedPath$INPUT $chr_list > $AlignedPath$ID.human.bam #AL modified filenames to remove .hg19 and add .human_unfiltered


echo 'mouse reads removed'

#AL note: next step is to bai index the file
samtools index $AlignedPath$ID.human.bam $AlignedPath$ID.human.bam.bai

echo 'index created'


# Remove remained RNEXT mouse reads using tag information

#AL notes:
#index($2,"_NCBI_GRCm38")==0 is intended to remove the mouse contigs from the header #i guess awk is 1 indexed
#index($7,"_NCBI_GRCm38")==0 removes reads where the read pair is aligned to the mouse genome. This column (RNEXT) has an = for pairs on the same contig
samtools view -h $AlignedPath$ID.human.bam | awk -F '\t' -vTAG=${tag}  '{ if (index($2,TAG)==0 && index($7,TAG)==0) print $0 }' | samtools view -bS - > $AlignedPath$ID.cleaned.bam

echo 'reads with mouse-aligned mates removed'
samtools index $AlignedPath$ID.cleaned.bam $AlignedPath$ID.cleaned.bam.bai

echo 'index created'

