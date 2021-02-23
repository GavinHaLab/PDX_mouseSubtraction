## Pipeline overview

    1) Prepare ConcatRef using ./scripts/concatenate_reference.py
    2) Unmap input BAM and subtract mouse reads using 'unmap_and_subtract_mouse.snakefile'
    3) Realign 'cleaned' bam with mouse reads removed to human reference genome using 'realign_bam_paired.snakefile'
    
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
    sample_name_ConcatRef_sorted.cleaned.bam ('cleaned' bam with mouse reads removed and mouse chromosomes removed from the header - this is the one you will want for the next step)


Third, run 'realign_bam_paired.snakefile'
        
      After finish 'unmap_and_subtract_mouse.snakefile', the ‘cleaned.bam’ file needs to be realigned to the human reference. 
      
      This is because the mouse subtraction snakemake produces an edited bam and I was unsure if this would work for tools that check the integrity of the bam file.


