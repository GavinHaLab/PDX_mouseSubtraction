## Pipeline overview

    1) Prepare ConcatRef using ./scripts/concatenate_reference.py. 
    2) Run subtract_mouse_and_realign.snakefile to unmap, subtract mouse, and realign the results (Separate realignment step no longer required)

## How to run pipeline

First, generate concatRef genome using this script (or use one of the existing concatRef files on fast)
   
    ./scripts/concatenate_reference.py
    
    Additional info and details about the ConcatRef is in the read_me file:
    ./scripts/read_me.txt

Second, run subtract_mouse_and_realign.snakefile to unmap, subtract mouse, and realign the results (Separate realignment step no longer required) 
    
    To run this step, there are 4 parameters that should be adjusted in the config:
    	1. input_reference_genome - if your input is a cram file, the reference genome for this file is required. If you have a bam file use 'null'
    	2. ConcatRef_genome - the concatenated mouse+human genome produced by ./scripts/concatenate_reference.py
    	3. human_reference_genome - the version of the human reference genome that will be used for the final realignment (recommended, but not required, to be the same genome version as the human part of the concatRef)
    	4. tag - the suffix that is used to denote mouse chromosomes in the concatRef. Any read mapping to a chromosome containing this tag in the name will be removed.
	
