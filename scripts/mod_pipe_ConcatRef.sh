#! /bin/bash
#$ -S /bin/bash
#$ -cwd

# basic argument
# PROJECT=$1 #not sure what this is, not used in script
# AlignedPath=$2 #path/to/bam_file

AlignedPath=$1 #path/to/bam_file/

# Candidate=$3 #not sure what this is, not used in script

# INPUT=$4 #bam_file.bam

INPUT=$2 #bam_file.bam #aligned to ConcatRef

# DIR=$5 #also not sure what this is, not used in script
# SCRIPT=$DIR/script/tool_ConcatRef #not sure about this either, not used

ID=${INPUT%%.bam*}

# print start time
# date

# 1. Extract hg38 aligned reads #AL modified hg19 to hg38

#altered this to get all chroms
# for ((i=1;i<23;i++));
# do
# 	#chr=chr$i.hg19
# 	chr=chr$i #AL mod, chromosomes are just called chr1, chr2 etc
# 	chr_list="$chr_list $chr"
# done

#get all chroms from the bam
#this command gets the input chromosomes from idxstats, cuts out the other columns, removes anything with _NCBI_GRCm38 or * 
chr_list=$(samtools idxstats  $AlignedPath$INPUT | cut -f 1 | grep -v _NCBI_GRCm38 | grep -v \* )

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
# # 2. Replace BAM header

# AL notes: 
# I think what this is doing is changing the chrom names from chr#.hg19 to chr##
# see post on this forum using extremely similar code:
# http://seqanswers.com/forums/showthread.php?t=22504
# however, I did not change the names of my human chromosomes when making the ConcatRef genome so I will not do this step

# samtools view -H $AlignedPath$ID.human.bam | sed -e 's/SN:\(chr[0-9MXY]\).hg19/SN:\1/' -e 's/SN:\(chr[0-9][0-9]\).hg19/SN:\1/' | samtools reheader - $AlignedPath$ID.hg19.bam > $AlignedPath$ID.hg19.reheadered.bam
# samtools index $AlignedPath$ID.hg19.reheadered.bam $AlignedPath$ID.hg19.reheadered.bai


# # 3. Remove remained RNEXT mm10 reads

#AL notes:
#index($2,"_NCBI_GRCm38")==0 is intended to remove the mouse contigs from the header #i guess awk is 1 indexed
#index($7,"_NCBI_GRCm38")==0 removes reads where the read pair is aligned to the mouse genome. This column (RNEXT) has an = for pairs on the same contig

samtools view -h $AlignedPath$ID.human.bam | awk -F '\t' '{ if (index($2,"_NCBI_GRCm38")==0 && index($7,"_NCBI_GRCm38")==0) print $0 }' | samtools view -bS - > $AlignedPath$ID.cleaned.bam

echo 'reads with mouse-aligned mates removed'

samtools index $AlignedPath$ID.cleaned.bam $AlignedPath$ID.cleaned.bam.bai

echo 'index created'


# # print end time
# date
