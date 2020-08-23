import glob 
import re

configfile: "config.json"

# from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider(access_key_id="", secret_access_key="")
#TODO add resources 

def split_f(string):
	return '/'.join(map(str,string.split('/')[-2:]))[:-16];

files=list(map(split_f,glob.glob("{}/*/*.fastq.gz".format(config["folder"]))))

# very ugly but works
rule all:
	input:
		expand("VCF/{files}.g.vcf",files=files[:5]),

#rule trim:
#        input:
#                R1=expand("{folder}/{{samples}}/{{file}}_R1_001.fastq.gz",folder = config["folder"]),
#                R2=expand("{folder}/{{samples}}/{{file}}_R2_001.fastq.gz",folder = config["folder"]),
#        output:
#                R1="trimmed/{samples}/{file}_R1_001.fastq.gz",
#                R2="trimmed/{samples}/{file}_R2_001.fastq.gz",
#        threads: 
#                workflow.cores * config["cores_per"] 
#        message: 
#                "trimming {wildcards.samples}/{wildcards.file}"
#        shell:
#                """
#                        fastp --detect_adapter_for_pe --overrepresentation_analysis --correction -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2}
#                """

rule trim_test:
	input:
		R1=expand("{folder}/{{samples}}/{{file}}_R1_001.fastq.gz",folder = config["folder"]),
		R2=expand("{folder}/{{samples}}/{{file}}_R2_001.fastq.gz",folder = config["folder"]),
	output:
		R1="trimmed/{samples}/{file}_R1_001.fastq.gz",
		R2="trimmed/{samples}/{file}_R2_001.fastq.gz",
	params:
		trimmomatic="/home/hubner/Trimmomatic-0.36/trimmomatic-0.36.jar",
		adapters="/home/hubner/Trimmomatic-0.36/adapters/AllAdapt.fa:2:30:10",
		leadin=3,
		trailing=3,
		slidingwindow="4:10",
		minlen=36,
	threads: 
		workflow.cores * config["cores_per"] 
	message: 
		"trimming {wildcards.samples}/{wildcards.file}"
	benchmark:
		"benchmarks/trim/{samples}/{file}.tsv"
	log:
		"logs/trim/{samples}/{file}.log"
	shell:
		"""
			java -jar {params.trimmomatic} PE\
					-threads {threads}\
					{input.R1} {input.R2}\
					{output.R1} /dev/null\
					{output.R2} /dev/null\
					ILLUMINACLIP:{params.adapters}\
					LEADIN:{params.leadin}\
					TRAILING:{params.trailing}\
					SLIDINGWINDOW:{params.slidingwindow}\
					MINLEN:{params.minlen}\
					> {log} 2>&1	
		"""

# check if file size is bigger

rule map:
	input:
		R1="trimmed/{samples}/{file}_R1_001.fastq.gz",
		R2="trimmed/{samples}/{file}_R2_001.fastq.gz",
		fasta=config["fasta"],
	output:
		"mapping/{samples}/{file}.bam"
	threads: 
		workflow.cores * config["cores_per"] 
	message: 
		"mapping {wildcards.file} to bam"
	benchmark:
		"benchmarks/map/{samples}/{file}.tsv"
	log:
		"logs/map/{samples}/{file}.log"
	shell:
		"""
			bwa mem -t {threads} -R "$(bash RG.sh {wildcards.file})" {input.fasta} {input.R1} {input.R2} |\
			samtools sort -@ {threads} --reference {input.fasta} -o {output}\
			> {log} 2>&1	
		"""

rule index:
	input:
		"mapping/{samples}/{file}.bam"
	output:
		"mapping/{samples}/{file}.bam.bai"
	threads: 
		workflow.cores * config["cores_per"] 
	message: 
		"indexing {wildcards.file}.bam"	
	benchmark:
		"benchmarks/index/{samples}/{file}.tsv"
	log:
		"logs/index/{samples}/{file}.log"
	shell:
		"""
			samtools.0.1.19 index {input}\
			> {log} 2>&1
		"""

rule mark_duplicates:
	input:
		bam="mapping/{samples}/{file}.bam",
		index_bam="mapping/{samples}/{file}.bam.bai"
	params:
		gatk="/mnt/data/Tools/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar"
	output:
		bam="mapping/{samples}/{file}.MarkedDup.bam",
		txt="mapping/{samples}/{file}-marked_dup_metrics.txt"
	threads: 
		workflow.cores * config["cores_per"] 
	message: 
		"MarkDuplicates {wildcards.file}"
	benchmark:
		"benchmarks/mark_dup/{samples}/{file}.tsv"
	log:
		"logs/mark_dup/{samples}/{file}.log"
	shell:
		# java -Xmx50G -jar {params.gatk} MarkDuplicatesSpark --spark-master local[{threads}]  -I {params.bam} -O {output.bam} -M {output.txt} 
		"""
			java -Xmx50G -jar {params.gatk} MarkDuplicates -I {input.bam} -O {output.bam} -M {output.txt}\
			> {log} 2>&1
		"""

rule haplotypecaller:
	input:
		MarkDup="mapping/{samples}/{file}.MarkedDup.bam",
		MarkDupIndex="mapping/{samples}/{file}.MarkedDup.bam.bai",
		fasta=config["fasta"]
	params:
		gatk="/mnt/data/tools/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar"
	output:
		"VCF/{samples}/{file}.g.vcf"
	threads: 
		workflow.cores * config["cores_per"] 
	message:
		"running HaplotypeCaller on {input.fasta}"
	benchmark:
		"benchmarks/haplotypecaller/{samples}/{file}.tsv"
	log:
		"logs/haplotypecaller/{samples}/{file}.log"
	shell:
		# java -Xmx50G -jar {params.gatk} HaplotypeCallerSpark --spark-master local[{threads}] --tmp-dir tmp  -R {params.fasta} -I {input.MarkDup} -L chr6H -O {output} -ERC GVCF
		""" 
			java -Xmx50G -jar {params.gatk} HaplotypeCaller -R {input.fasta} -I {input.MarkDup} -L chr6H1-100 -O {output} -ERC GVCF\
			> {log} 2>&1
		"""
#TODO:
#1 	-L chr6H 
#2 	 adjust the heterozygosity parameter  -- to ????? 

rule sample_map:
	input:
		vcfs=expand("VCF/{{samples}}/{vcf}.g.vcf",vcf=files)
	output:
		map_vcf="VCF/{samples}/cohort.sample_map",
	threads: 
		workflow.cores * config["cores_per"] ,
	message: 
		"creating sample map"
	shell:
		""" 
			echo {input} > {output}
		"""


rule GenomicsDBImport:
	input:
		map_vcf="VCF/{samples}/cohort.sample_map",
	output:
		path=directory("genomicsdb/{samples}/"),
	params:
		gatk="/mnt/data/Tools/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar",
		batch_size=2000,
	threads: 
		workflow.cores * config["cores_per"] ,
	message: 
		""
	shell:
		""" 
			java -jar {params.gatk} --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
					--genomicsdb-workspace-path {output.path}\
					--batch-size {params.batch_size} \
					-L {} \
					--sample-name-map {input.map_vcf}\
					--tmp-dir={} \
					--reader-threads {threads}
		"""

#java -Xmx50G -jar {params.gatk} -V 3GenomicsDBImport --genomicsdb-workspace-path {output.path} --samples-name-map {output.name_map}

#rule GenotypeGVCFs:
	#input:
		#output:
			#params:
				#gatk="/mnt/data/Tools/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar"
	#threads: 
		#workflow.cores * config["cores_per"] 
	#shell:
		#""" 
			#java -Xmx50G -jar {params.gatk} GenotypeGVCFs -R {input.ref} -V gendb://{input.??} -G StandardAnnotation -O {output} 
		#"""


#rule splitGVCF:
	#input:
		#bed=config["bed"],
		#id_list=config["id"],
		#vcf_folder=config["vcf_path"],
		#ref=config["fasta_ref"],
	#output:
		#config["output_path"],	
	#params:
		#splits=30000,

	#threads: workflow.cores * config["cores_per"] 

	#shell:
		#"""
			#bash SplitGVCF.v3.sh {params.splits} {input.bed} {input.id_list} {input.vcf_folder} {output} {threads} {input.ref}	
		#"""
