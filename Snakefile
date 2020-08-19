import glob 
import re

configfile: "config.json"

# from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider(access_key_id="", secret_access_key="")

#TODO add resources 
def split_f(string):
	return '/'.join(map(str,string.split('/')[-2:]))[:-9];

files=glob.glob("{}/*/*.fastq.gz".format(config["folder"]))
files=list(map(split_f,files))


rule all:
	input:
		expand("VCF/{files}.g.vcf",files=files)

#TODO fix rule all input is acausally output

rule trim:
	input:
		expand("{folder}/{{samples}}/{{file}}.fastq.gz",folder = config["folder"]),
	output:
		"trimmed/{samples}/{file}.fastq.gz"
	threads: 
		workflow.cores * config["cores_per"] 
	message: 
		"trimming {wildcards.samples}/{wildcards.file}"
	shell:
		"""
			fastp --detect_adapter_for_pe --overrepresentation_analysis -i {input} -o {output} 
		"""

#fastp --detect_adapter_for_pe --overrepresentation_analysis --correction -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]}

rule map:
	input:
		files="trimmed/{samples}/{file}.fastq.gz",
		fasta=config["fasta"],
	output:
		"mapping/{samples}/{file}.bam"
	threads: 
		workflow.cores * config["cores_per"] 
	message: 
		"mapping {wildcards.file}.fastq.gz to bam"
	shell:
		"""
			bwa mem -t {threads} -R "$(bash RG.sh {input.files})" {input.fasta} {input.files} | samtools sort -@ {threads} --reference {input.fasta} -o {output}
		"""

#TODO change to correct sample name in @RG -- done

rule index:
	input:
		"mapping/{samples}/{file}.bam"
	output:
		"mapping/{samples}/{file}.bam.bai"
	threads: 
		workflow.cores * config["cores_per"] 
	message: 
		"indexing {wildcards.file}.bam"	
	shell:
		"""
			samtools.0.1.19 index {input}
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
	shell:
		# java -Xmx50G -jar {params.gatk} MarkDuplicatesSpark --spark-master local[{threads}]  -I {params.bam} -O {output.bam} -M {output.txt} 
		"""
			java -Xmx50G -jar {params.gatk} MarkDuplicates -I {input.bam} -O {output.bam} -M {output.txt} -VALIDATION_STRINGENCY SILENT
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
	shell:
		# java -Xmx50G -jar {params.gatk} HaplotypeCallerSpark --spark-master local[{threads}] --tmp-dir tmp  -R {params.fasta} -I {input.MarkDup} -L chr6H -O {output} -ERC GVCF
		""" 
			java -Xmx50G -jar {params.gatk} HaplotypeCaller -R {input.fasta} -I {input.MarkDup} -O {output} -ERC GVCF
		"""
#-L chr6H 
# adjust the heterozygosity parameter  -- to ????? 

#rule GenomicsDBImport:
	#input:
	#output:
		#path="",
		#name_map="",
	#params:
		#gatk="/mnt/data/Tools/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar"
	#threads: 
		#workflow.cores * config["cores_per"] 
	##message: ""
	#shell:
		#""" 
			#java -Xmx50G -jar {params.gatk} GenomicsDBImport --genomicsdb-workspace-path {output.path} --samples-name-map {output.name_map}
		#"""

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
