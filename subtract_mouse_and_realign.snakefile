#subtract_mouse_and_realign.snakefile
#Anna-Lisa Doebley
#Template made 2021-09-24
#Ha Lab
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml BWA/0.7.17-foss-2018b
ml SAMtools/1.10-GCCcore-8.3.0
ml picard/2.18.29-Java
ml Java/11.0.2
ml GATK/4.1.4.1-GCCcore-8.3.0-Java-11
ml R/3.6.2-foss-2019b-fh1


#command to run snakemake (remove -np at end when done validating):
snakemake -s subtract_mouse_and_realign.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np

#output file marked as temp is deleted after all rules that use it as an input are completed
"""

configfile: "config/samples.yaml"
configfile: "config/config.yaml"

rule all:
	input:
		#concatRef metrics
		expand(config['results_path']+"/{samples}/{samples}_ConcatRef_ASM.txt", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_ConcatRef_wgs_metrics.txt", samples=config["samples"]),
		#results from realignment
		expand(config['results_path']+"/{samples}/{samples}_marked_dup_metrics.txt", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_recalibration_data.table", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_recalibrated.bam", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_recalibrated.bam.bai", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_ASM.txt", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_wgs_metrics.txt", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_isize.txt", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_isize.pdf", samples=config["samples"])


rule unmap:
	input:
		bam_file = lambda wildcards: config["samples"][wildcards.samples]
	output:
		fastq1 = temp(config['results_path']+"/concatRef_fastqs/{samples}_fastq1.fq.gz"),
		fastq2 = temp(config['results_path']+"/concatRef_fastqs/{samples}_fastq2.fq.gz"),
		unpaired_fastq = temp(config['results_path']+"/concatRef_fastqs/{samples}_unpaired.fq.gz")
	params:
		java=config["java"],
		picard_jar = config["picard_jar"],
		input_reference_genome = config["input_reference_genome"] #required for cram input
	shell:
		"{params.java} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx85G \
		-jar {params.picard_jar} SamToFastq I={input.bam_file} \
		TMP_DIR="+config['results_path']+"/concatRef_tmps \
		REFERENCE_SEQUENCE={params.input_reference_genome} \
		FASTQ={output.fastq1} \
		SECOND_END_FASTQ={output.fastq2} \
		UNPAIRED_FASTQ={output.unpaired_fastq}"


rule map_to_ConcatRef:
	input:
		fastq1=config['results_path']+"/concatRef_fastqs/{samples}_fastq1.fq.gz",
		fastq2=config['results_path']+"/concatRef_fastqs/{samples}_fastq2.fq.gz"
	output:
		temp(config['results_path']+"/{samples}/{samples}_ConcatRef_unsorted.bam")
	params:
		reference_genome=config["ConcatRef_genome"],
		bwa=config["bwa"],
		samtools=config["samtools"],
		bwa_threads=config["bwa_threads"]
	shell:
		"({params.bwa} mem -t {params.bwa_threads} -M \
		-R '@RG\\tID:no_id\\tLB:no_library\\tPL:no_platform\\tPU:no_unit\\tSM:{wildcards.samples}' \
		{params.reference_genome} \
		{input.fastq1} {input.fastq2} | {params.samtools} view -b - > {output})"


rule sort_concatRef_by_coord:
	input:
		config['results_path']+"/{samples}/{samples}_ConcatRef_unsorted.bam"
	output:
		temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam")
	params:
		samtools=config["samtools"]
	shell:
		"({params.samtools} sort -o {output} {input})"


rule index_ConcatRef_bam:
	input:
		sorted_bam = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam"
	output:
		sorted_bam_index = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam.bai")
	params:
		samtools=config["samtools"]
	shell:
		"{params.samtools} index {input.sorted_bam} {output.sorted_bam_index}"


rule get_ConcatRef_wgs_metrics:
	input:
		sorted_bam = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam",
		bai_index = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam.bai"
	output:
		config['results_path']+"/{samples}/{samples}_ConcatRef_wgs_metrics.txt"
	params:
		reference_genome=config["ConcatRef_genome"],
		java=config["java"],
		picard_jar=config["picard_jar"]
	run: 
		shell("{params.java} -jar {params.picard_jar} CollectWgsMetrics \
		I={input.sorted_bam} \
		O={output} \
		R={params.reference_genome} \
		TMP_DIR="+config['results_path']+"/concatRef_tmps \
		COUNT_UNPAIRED=true \
		USE_FAST_ALGORITHM=true \
		INCLUDE_BQ_HISTOGRAM=true")
		

rule get_ConcatRef_alignment_metrics:
	input:
		sorted_bam = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam",
		bai_index = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam.bai"
	output:
		config['results_path']+"/{samples}/{samples}_ConcatRef_ASM.txt"
	params:
		reference_genome=config["ConcatRef_genome"],
		java=config["java"],
		picard_jar=config["picard_jar"]
	shell:
		"{params.java} -jar {params.picard_jar} CollectAlignmentSummaryMetrics \
		R={params.reference_genome} \
		I={input.sorted_bam} \
		O={output} \
		TMP_DIR="+config['results_path']+"/concatRef_tmps"

rule remove_mouse:
	input:
		sorted_bam = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam",
		bai_index = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam.bai"
	output:
		sorted_bam_step_1 = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.human.bam"),
		sorted_bam_step_1_index = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.human.bam.bai"),
		sorted_clean_bam = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.cleaned.bam"),
		sorted_clean_bam_index = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.cleaned.bam.bai")
	params:
		path = config['results_path']+"/{samples}/",
		filename = "{samples}_ConcatRef_sorted.bam",
		tag = config["tag"]
	shell:
		"sh ./scripts/mod_pipe_ConcatRef.sh {params.path} {params.filename} {params.tag}"


#the following rules realign post mouse subtraction
rule unmap_cleaned:
	input:
		bam_file = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.cleaned.bam"
	output:
		fastq1 = temp(config['results_path']+"/human_fastqs/{samples}_fastq1.fq.gz"),
		fastq2 = temp(config['results_path']+"/human_fastqs/{samples}_fastq2.fq.gz"),
		fastqUnpaired = temp(config['results_path']+"/human_fastqs/{samples}_fastq_unpaired.fq.gz")
	params:
		java=config["java"],
		picard_jar = config["picard_jar"]
	shell:
		"{params.java} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx70G -jar \
		{params.picard_jar} SamToFastq \
		I={input.bam_file} \
		TMP_DIR="+config['results_path']+"/human_tmps \
		FASTQ={output.fastq1} \
		SECOND_END_FASTQ={output.fastq2} \
		UNPAIRED_FASTQ={output.fastqUnpaired}"

rule map_to_human_ref:
	input:
		fastq1=config['results_path']+"/human_fastqs/{samples}_fastq1.fq.gz",
		fastq2=config['results_path']+"/human_fastqs/{samples}_fastq2.fq.gz"
	output:
		temp(config['results_path']+"/{samples}/{samples}_unsorted.bam")
	params:
		reference_genome=config["human_reference_genome"],
		bwa=config["bwa"],
		samtools=config["samtools"],
		bwa_threads=config["bwa_threads"]
	shell:
		"({params.bwa} mem -t {params.bwa_threads} -M \
		-R '@RG\\tID:no_id\\tLB:no_library\\tPL:no_platform\\tPU:no_unit\\tSM:{wildcards.samples}' \
		{params.reference_genome} \
		{input.fastq1} {input.fastq2} | {params.samtools} view -b - > {output})"

rule sort_by_coord_cleaned:
	input:
		config['results_path']+"/{samples}/{samples}_unsorted.bam"
	output:
		temp(config['results_path']+"/{samples}/{samples}_sorted.bam")
	params:
		samtools=config["samtools"]
	shell:
		"({params.samtools} sort -o {output} {input})"


rule mark_dups_cleaned:
	input:
		config['results_path']+"/{samples}/{samples}_sorted.bam"
	output:
		bam=temp(config['results_path']+"/{samples}/{samples}_dups_marked.bam"),
		metrics=protected(config['results_path']+"/{samples}/{samples}_marked_dup_metrics.txt")
	params:
		java=config["java"],
		picard_jar = config["picard_jar"]
	shell:
		"({params.java} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx72G \
		-jar {params.picard_jar} MarkDuplicates \
		I={input} \
		O={output.bam} \
		M={output.metrics} \
		TMP_DIR="+config['results_path']+"/human_tmps)"


rule build_recalibrator_model_cleaned:
	input:
		config['results_path']+"/{samples}/{samples}_dups_marked.bam"
	output:
		protected(config['results_path']+"/{samples}/{samples}_recalibration_data.table")
	params:
		gatk=config["gatk"],
		reference_genome=config["human_reference_genome"],
		known_polymorphic_sites1=config["known_polymorphic_sites1"],
		known_polymorphic_sites2=config["known_polymorphic_sites2"],
		base_recalibrator_gap_open_penalty=config["base_recalibrator_gap_open_penalty"]
	shell:
		"{params.gatk} BaseRecalibrator \
		-I {input} \
		-R {params.reference_genome} \
		--known-sites {params.known_polymorphic_sites1} \
		--known-sites {params.known_polymorphic_sites2} \
		--bqsr-baq-gap-open-penalty {params.base_recalibrator_gap_open_penalty} \
		-O {output} \
		--tmp-dir "+config['results_path']+"/human_tmps"


rule apply_recalibration_cleaned:
	input:
		bam=config['results_path']+"/{samples}/{samples}_dups_marked.bam",
		model=config['results_path']+"/{samples}/{samples}_recalibration_data.table"
	output:
		bam=protected(config['results_path']+"/{samples}/{samples}_recalibrated.bam"),
		index=config['results_path']+"/{samples}/{samples}_recalibrated.bai"
	params:
		gatk=config["gatk"],
		reference_genome=config["human_reference_genome"]
	shell:
		"({params.gatk} ApplyBQSR \
		-R {params.reference_genome} \
		-I {input.bam} \
		--bqsr-recal-file {input.model} \
		-O {output.bam})"


rule rename_index_files_cleaned:
	input:
		config['results_path']+"/{samples}/{samples}_recalibrated.bai"
	output:
		protected(config['results_path']+"/{samples}/{samples}_recalibrated.bam.bai")
	shell:
		"(mv {input} {output})"


rule get_alignment_metrics_cleaned:
	input:
		config['results_path']+"/{samples}/{samples}_recalibrated.bam"
	output:
		protected(config['results_path']+"/{samples}/{samples}_ASM.txt")
	params:
		reference_genome=config["human_reference_genome"],
		java=config["java"],
		picard_jar=config["picard_jar"]
	shell:
		"({params.java} -jar {params.picard_jar} CollectAlignmentSummaryMetrics \
		R={params.reference_genome} \
		I={input} \
		O={output})"


rule get_wgs_metrics_cleaned:
	input:
		config['results_path']+"/{samples}/{samples}_recalibrated.bam"
	output:
		protected(config['results_path']+"/{samples}/{samples}_wgs_metrics.txt")
	params:
		reference_genome=config["human_reference_genome"],
		is_wgs=config["is_wgs"],
		java=config["java"],
		picard_jar=config["picard_jar"]
	run:
		if params.is_wgs:
			shell("({params.java} -jar {params.picard_jar} CollectWgsMetrics \
			I={input} \
			O={output} \
			R={params.reference_genome} \
			COUNT_UNPAIRED=true \
			USE_FAST_ALGORITHM=true \
			INCLUDE_BQ_HISTOGRAM=true)")

rule get_isize_cleaned:
	input:
		config['results_path']+"/{samples}/{samples}_recalibrated.bam"
	output:
		isize_txt = config['results_path']+"/{samples}/{samples}_isize.txt",
		isize_pdf = config['results_path']+"/{samples}/{samples}_isize.pdf"
	params:
		java=config["java"],
		picard_jar=config["picard_jar"]
	shell:
 		"time {params.java} -jar {params.picard_jar} CollectInsertSizeMetrics \
 		I={input} \
 		O={output.isize_txt} \
 		H={output.isize_pdf} \
 		HISTOGRAM_WIDTH=500 \
 		TMP_DIR="+config['results_path']+"/human_tmps"

