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

#command to run snakemake (remove -np at end when done validating):
snakemake -s unmap_and_subtract_mouse.snakefile --latency-wait 60 --restart-times 3 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np
"""

configfile: "config/samples.yaml"
configfile: "config/config.yaml"

rule all:
	input:
		expand("results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.bam", samples=config["samples"]),
		expand("results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.bam.bai", samples=config["samples"]),
		expand("results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.cleaned.bam", samples=config["samples"]),
		expand("results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.cleaned.bam.bai", samples=config["samples"]),
		expand("results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted_idxstats.txt", samples=config["samples"])

rule unmap:
	input:
		bam_file = lambda wildcards: config["samples"][wildcards.samples]
	output:
		fastq1 = "results/subtract_mouse/fastqs/{samples}_fastq1.fq.gz",
		fastq2 = "results/subtract_mouse/fastqs/{samples}_fastq2.fq.gz",
		unpaired_fastq = "results/subtract_mouse/fastqs/{samples}_unpaired.fq.gz"
	params:
		java=config["java"],
		picard_jar = config["picard_jar"],
		reference_genome = config["reference_genome"]
	log:
		"logs/subtract_mouse/unmap/{samples}_unmap.txt"
	shell:
		"{params.java} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx85G \
		-jar {params.picard_jar} SamToFastq I={input.bam_file} \
		TMP_DIR=results/tmps/{wildcards.samples}_tmp \
		REFERENCE_SEQUENCE= {params.reference_genome} \
		FASTQ={output.fastq1} \
		SECOND_END_FASTQ={output.fastq2} \
		UNPAIRED_FASTQ={output.unpaired_fastq}"


rule map_to_ConcatRef:
	input:
		fastq1="results/subtract_mouse/fastqs/{samples}_fastq1.fq.gz",
		fastq2="results/subtract_mouse/fastqs/{samples}_fastq2.fq.gz"
	output:
		temp("results/subtract_mouse/{samples}/{samples}_ConcatRef_unsorted.bam")
	params:
		reference_genome=config["ConcatRef_genome"],
		bwa=config["bwa"],
		samtools=config["samtools"],
		bwa_threads=config["bwa_threads"]
	log:
		"logs/subtract_mouse/map_to_reference/{samples}_map_to_ConcatRef_reference.txt"
	shell:
		"({params.bwa} mem -t {params.bwa_threads} -M \
		-R '@RG\\tID:no_id\\tLB:no_library\\tPL:no_platform\\tPU:no_unit\\tSM:{wildcards.samples}' \
		{params.reference_genome} \
		{input.fastq1} {input.fastq2} | {params.samtools} view -b - > {output}) 2> {log}"


rule sort_by_coord:
	input:
		"results/subtract_mouse/{samples}/{samples}_ConcatRef_unsorted.bam"
	output:
		"results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.bam"
	params:
		samtools=config["samtools"]
	log:
		"logs/subtract_mouse/sort_by_coord/{samples}_ConcatRef_sort_by_coord.txt"
	shell:
		"({params.samtools} sort -o {output} {input}) 2> {log}"

rule index_bam:
	input:
		sorted_bam = "results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.bam"
	output:
		sorted_bam_index = "results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.bam.bai"

	params:
		samtools=config["samtools"]
	shell:
		"{params.samtools} index {input.sorted_bam} {output.sorted_bam_index}"

rule remove_mouse:
	input:
		sorted_bam = "results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.bam",
		bai_index = "results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.bam.bai"

	output:
		sorted_bam_step_1 = temp("results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.human.bam"),
		sorted_bam_step_1_index = temp("results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.human.bam.bai"),
		sorted_clean_bam = "results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.cleaned.bam",
		sorted_clean_bam_index = "results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.cleaned.bam.bai"
	params:
		path = "results/subtract_mouse/{samples}/",
		filename = "{samples}_ConcatRef_sorted.bam",
		tag = config["tag"]
	shell:
		"sh ./scripts/mod_pipe_ConcatRef.sh {params.path} {params.filename} {params.tag}"

rule get_idxstats:
	input:
		sorted_bam = "results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.bam",
		bai_index = "results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.bam.bai"
	output:
		idx_stats = "results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted_idxstats.txt"
	shell:
		"samtools idxstats {input.sorted_bam} > {output.idx_stats}"

rule get_wgs_metrics:
	input:
		"results/subtract_mouse/{samples}/{samples}_ConcatRef_sorted.cleaned.bam"
	output:
		protected("results/subtract_mouse/{samples}/{samples}_wgs_metrics.txt")
	params:
		reference_genome=config["reference_genome"],
		is_wgs=config["is_wgs"],
		java=config["java"],
		picard_jar=config["picard_jar"]
	log:
		"logs/wgs_metrics/subtract_mouse/{samples}_get_wgs_metrics.txt"
	run:
		if params.is_wgs:
			shell("({params.java} -jar {params.picard_jar} CollectWgsMetrics \
			I={input} \
			O={output} \
			R={params.reference_genome} \
			COUNT_UNPAIRED=true \
			USE_FAST_ALGORITHM=true \
			INCLUDE_BQ_HISTOGRAM=true) 2> {log}")
			