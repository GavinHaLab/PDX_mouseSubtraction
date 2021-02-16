#unmap_and_subtract_mouse.snakefile
#Anna-Lisa Doebley
#Template made 2020-04-27
#Ha Lab
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml BWA/0.7.17-foss-2018b
ml SAMtools/1.10-GCCcore-8.3.0
# ml java/jdk1.8.0_31 #not available on bionic
ml picard/2.18.29-Java
export PATH="$PATH:/fh/fast/ha_g/user/adoebley/ha_lab_scripts/bin/PDX_mouse_subtraction/"

#command to run snakemake (remove -np at end when done validating):
snakemake -s unmap_and_subtract_mouse.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np
"""

configfile: "config/samples.yaml"
configfile: "config/config.yaml"

rule all:
	input:
		expand("results/{samples}/{samples}_ConcatRef_sorted.bam", samples=config["samples"]),
		expand("results/{samples}/{samples}_ConcatRef_sorted.bam.bai", samples=config["samples"]),
		expand("results/{samples}/{samples}_ConcatRef_sorted.cleaned.bam", samples=config["samples"]),
		expand("results/{samples}/{samples}_ConcatRef_sorted.cleaned.bam.bai", samples=config["samples"]),
		expand("results/{samples}/{samples}_ConcatRef_wgs_metrics.txt", samples=config["samples"]),
		expand("results/{samples}/{samples}_ConcatRef_ASM_metrics.txt", samples=config["samples"])


rule unmap:
	input:
		bam_file = lambda wildcards: config["samples"][wildcards.samples]
	output:
		fastq1 = temp("results/fastqs/{samples}_fastq1.fq.gz"),
		fastq2 = temp("results/fastqs/{samples}_fastq2.fq.gz"),
		fastqUnpaired = temp("results/fastqs/{samples}_fastq_unpaired.fq.gz")
	params:
		java=config["java"],
		picard_jar = config["picard_jar"]
	log:
		"logs/unmap/{samples}_unmap.txt"
	shell:
		"{params.java} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx70G -jar \
		{params.picard_jar} SamToFastq \
		I={input.bam_file} \
		TMP_DIR=results/tmps/{wildcards.samples}_tmp \
		FASTQ={output.fastq1} \
		SECOND_END_FASTQ={output.fastq2} \
		UNPAIRED_FASTQ={output.fastqUnpaired}"


rule map_to_ConcatRef:
	input:
		fastq1="results/fastqs/{samples}_fastq1.fq.gz",
		fastq2="results/fastqs/{samples}_fastq2.fq.gz"
	output:
		temp("results/{samples}/{samples}_ConcatRef_unsorted.bam")
	params:
		reference_genome=config["ConcatRef_genome"],
		bwa=config["bwa"],
		samtools=config["samtools"],
		bwa_threads=config["bwa_threads"]
	log:
		"logs/map_to_reference/{samples}_map_to_ConcatRef_reference.txt"
	shell:
		"({params.bwa} mem -t {params.bwa_threads} -M \
		-R '@RG\\tID:no_id\\tLB:no_library\\tPL:no_platform\\tPU:no_unit\\tSM:{wildcards.samples}' \
		{params.reference_genome} \
		{input.fastq1} {input.fastq2} | {params.samtools} view -b - > {output}) 2> {log}"


rule sort_by_coord:
	input:
		"results/{samples}/{samples}_ConcatRef_unsorted.bam"
	output:
		"results/{samples}/{samples}_ConcatRef_sorted.bam"
	params:
		samtools=config["samtools"]
	log:
		"logs/sort_by_coord/{samples}_ConcatRef_sort_by_coord.txt"
	shell:
		"({params.samtools} sort -o {output} {input}) 2> {log}"

rule index_bam:
	input:
		sorted_bam = "results/{samples}/{samples}_ConcatRef_sorted.bam"
	output:
		sorted_bam_index = "results/{samples}/{samples}_ConcatRef_sorted.bam.bai"

	params:
		samtools=config["samtools"]
	shell:
		"{params.samtools} index {input.sorted_bam} {output.sorted_bam_index}"

rule remove_mouse:
	input:
		sorted_bam = "results/{samples}/{samples}_ConcatRef_sorted.bam",
		bai_index = "results/{samples}/{samples}_ConcatRef_sorted.bam.bai"

	output:
		sorted_bam_step_1 = temp("results/{samples}/{samples}_ConcatRef_sorted.human.bam"),
		sorted_bam_step_1_index = temp("results/{samples}/{samples}_ConcatRef_sorted.human.bam.bai"),
		sorted_clean_bam = "results/{samples}/{samples}_ConcatRef_sorted.cleaned.bam",
		sorted_clean_bam_index = "results/{samples}/{samples}_ConcatRef_sorted.cleaned.bam.bai"
	params:
		path = "results/{samples}/",
		filename = "{samples}_ConcatRef_sorted.bam"
	shell:
		"mod_pipe_ConcatRef.sh {params.path} {params.filename}"

rule get_wgs_metrics:
	input:
		sorted_bam = "results/{samples}/{samples}_ConcatRef_sorted.bam",
		bai_index = "results/{samples}/{samples}_ConcatRef_sorted.bam.bai"
	output:
		"results/{samples}/{samples}_ConcatRef_wgs_metrics.txt"
	params:
		reference_genome=config["ConcatRef_genome"],
		java=config["java"],
		picard_jar=config["picard_jar"]
	run:
		shell("{params.java} -jar {params.picard_jar} CollectWgsMetrics \
		I={input.sorted_bam} \
		O={output} \
		R={params.reference_genome} \
		TMP_DIR=results/tmps/{wildcards.samples}_tmp \
		COUNT_UNPAIRED=true \
		USE_FAST_ALGORITHM=true \
		INCLUDE_BQ_HISTOGRAM=true")
		
rule get_alignment_metrics:
	input:
		sorted_bam = "results/{samples}/{samples}_ConcatRef_sorted.bam",
		bai_index = "results/{samples}/{samples}_ConcatRef_sorted.bam.bai"
	output:
		"results/{samples}/{samples}_ConcatRef_ASM_metrics.txt"
	params:
		reference_genome=config["ConcatRef_genome"],
		java=config["java"],
		picard_jar=config["picard_jar"]

	shell:
		"{params.java} -jar {params.picard_jar} CollectAlignmentSummaryMetrics \
		R={params.reference_genome} \
		I={input.sorted_bam} \
		O={output} \
		TMP_DIR=results/tmps/{wildcards.samples}_tmp"