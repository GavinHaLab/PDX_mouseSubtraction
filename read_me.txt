Notes about 'unmap_and_subtract_mouse_snakemake' and how to use it:

First, the concatenated reference genome is here:

    /fh/fast/ha_g/grp/reference/ConcatRef/fasta/GRCh38_plus_NCBI_GRCm38.fa


    This same folder also contains the bwa index files necessary for running the bwa aligner.


    I generated the concatRef genome using this script:

    /fh/fast/ha_g/grp/reference/ConcatRef/scripts/concatenate_GRCh38_plus_NCBI_GRCm38.py


    Currently this script is hard coded to use a specific mouse and human genome:

    human_ref_path = '/fh/fast/ha_g/grp/reference/GRCh38/GRCh38.fa'

    mouse_ref_path = '/fh/fast/ha_g/grp/reference/Mus_musculus_NCBI_GRCm38/NCBI/GRCm38/Sequence/WholeGenomeFasta/genome.fa'



    Additional info and details about the ConcatRef is in the read_me file:

    /fh/fast/ha_g/grp/reference/ConcatRef/read_me.txt



Second, to perform the mouse subtraction, I used the script here:

    /fh/fast/ha_g/user/adoebley/ha_lab_scripts/bin/PDX_mouse_subtraction/mod_pipe_ConcatRef.sh

    This is a modified version of the pipe_ConcatRef.sh that you showed me on github. In the script, I made notes on what everything does and all the changes I made so that it worked with my data.


Third, I made a snakemake:

    This snakemake removes unpaired reads, unmaps from the previous genome, remaps to concatRef, sorts the bam, indexes the bam, then subtracts mouse, and finally runs samtools idxstats (to get an idea of the percent mouse for each sample). The snakemake is located here:


    /fh/fast/ha_g/user/adoebley/ha_lab_scripts/src/snakemakes/unmap_and_subtract_mouse_snakemake

    This snakemake outputs two bam files:
    sample_name_ConcatRef_sorted.bam (an intermediate file containing reads aligned to the concat_ref)
    sample_name _ConcatRef_sorted.cleaned.bam ('cleaned' bam with mouse reads removed and mouse chromosomes removed from the header - this is the one you will want for the next step)


Fourth:

    After running this snakemake, I realigned the ‘cleaned’ bam to the human reference genome using the ‘realign_bam_paired_snakemake’ that I sent you on Monday. This is because the mouse subtraction snakemake produces an edited bam and I was unsure if this would work for tools that check the integrity of the bam file.



In summary:

    ConcatRef is here:
    /fh/fast/ha_g/grp/reference/ConcatRef/fasta/GRCh38_plus_NCBI_GRCm38.fa

    Mouse subtraction script is here:
    /fh/fast/ha_g/user/adoebley/ha_lab_scripts/bin/PDX_mouse_subtraction/mod_pipe_ConcatRef.sh

    Snakemake to run the mouse subtraction pipeline is here:
    /fh/fast/ha_g/user/adoebley/ha_lab_scripts/src/snakemakes/unmap_and_subtract_mouse_snakemake

    After running the pipeline, the ‘cleaned.bam’ file needs to be realigned to the human reference. 