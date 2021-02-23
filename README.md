## This pipeline consists of two parts :

    1) Unmap input BAM and subtract mouse reads ('unmap_and_subtract_mouse.snakefile')
    2) Realign mouse reads removed edited bam to human reference genome ('realign_bam_paired.snakefile')
    
    Please run 'unmap_and_subtract_mouse.snakefile' first and then run 'realign_bam_paired.snakefile' afterwards.

## How to run pipeline:

First, generate concatRef genome using this script
   
    ./scripts/concatenate_reference.py
    
    Additional info and details about the ConcatRef is in the read_me file:
    ./scripts/read_me.txt

Second, run 'unmap_and_subtract_mouse.snakefile'

    This snakemake unmaps Input BAM and remaps to concatRef, sorts the bam, indexes the bam, then subtracts mouse, and finally runs samtools idxstats (to get an idea of the percent mouse for each sample). 
    To run the snakemake, update the config files then follow the directions at the top of the snakefile. 

    This snakemake outputs two bam files:
    sample_name_ConcatRef_sorted.bam (an intermediate file containing reads aligned to the ConcatRef)
    sample_name _ConcatRef_sorted.cleaned.bam ('cleaned' bam with mouse reads removed and mouse chromosomes removed from the header - this is the one you will want for the next step)


Third, run 'realign_bam_paired.snakefile'
        
      After finish 'unmap_and_subtract_mouse.snakefile', the ‘cleaned.bam’ file needs to be realigned to the human reference. 
      Realign the ‘cleaned’ bam to the human reference genome using the ‘realign_bam_paired_snakemake’. 
      This is because the mouse subtraction snakemake produces an edited bam and I was unsure if this would work for tools that check the integrity of the bam file.


