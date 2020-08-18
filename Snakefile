configfile: "config.json"

# from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider(access_key_id="", secret_access_key="")

rule all:
	input:
		expand("{folder}/data/{file}.fastq.gz",folder = config["folder"],file = config["files"]),
		config["fasta"],
rule trim:
	input:
		expand("{folder}/data/{file}.fastq.gz",folder = config["folder"],file = config["files"]),
	output:
		expand("trimmed/{file}.fastq.gz",file = config["files"]),

	threads: workflow.cores * config["cores_per"] 
	#message: ""
	shell:
		"""
			fastp --detect_adapter_for_pe --overrepresentation_analysis --correction -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]}
		"""
rule map:
	input:
		files=expand("trimmed/{file}.fastq.gz",file = config["files"]),
		fasta=config["fasta"]
	output:
		"mapping/sorted.bam"
	threads: workflow.cores * config["cores_per"] 
	#message: ""
	shell:
		"""
			bwa mem -t {threads} -R '@RG\\tID:ERR753110\\tSM:ERR753110\\tPL:ILLUMINA' {input.fasta} {input.files} | samtools sort -@ {threads} --reference {input.fasta} -o {output}
		"""
rule index:
	input:
		"mapping/{samples}.bam"
	output:
		"mapping/{samples}.bam.bai"
	threads: workflow.cores * config["cores_per"] 
	#message: ""
	shell:
		"""
			samtools.0.1.19 index {input}
		"""

rule mark_duplicates:
	input:
		bam="mapping/sorted.bam",
		index_bam="mapping/sorted.bam.bai"
	params:
		gatk="/mnt/data/Tools/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar"
	output:
		bam="mapping/sorted.MarkedDup.bam",
		txt="mapping/marked_dup_metrics.txt"
	threads: workflow.cores * config["cores_per"] 
	#message: ""
	shell:
		# java -Xmx50G -jar {params.gatk} MarkDuplicatesSpark --spark-master local[{threads}]  -I {params.bam} -O {output.bam} -M {output.txt} 
		"""
			java -Xmx50G -jar {params.gatk} MarkDuplicates -I {input.bam} -O {output.bam} -M {output.txt} -VALIDATION_STRINGENCY SILENT
		"""
rule haplotypecaller:
	input:
		MarkDup="mapping/sorted.MarkedDup.bam",
		MarkDupIndex="mapping/sorted.MarkedDup.bam.bai",
		fasta=config["fasta"]
	params:
		gatk="/mnt/data/tools/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar"
	output:
		"VCF/ERR753110.vcf"
	threads: workflow.cores * config["cores_per"] 
	#message: ""
	shell:
		# java -Xmx50G -jar {params.gatk} HaplotypeCallerSpark --spark-master local[{threads}] --tmp-dir tmp  -R {params.fasta} -I {input.MarkDup} -L chr6H -O {output} -ERC GVCF
		""" 
			java -Xmx50G -jar {params.gatk} HaplotypeCaller -R {input.fasta} -I {input.MarkDup} -L chr6H -O {output} -ERC GVCF
		"""
rule splitGVCF:
	input:
		bed=config["bed"],
		id_list=config["id"],
		vcf_folder=config["vcf_path"],
		ref=config["fasta_ref"],
	output:
		config["output_path"],	
	params:
		splits=30000,

	threads: workflow.cores * config["cores_per"] 

	shell:
		"""
			bash SplitGVCF.v3.sh {params.splits} {input.bed} {input.id_list} {input.vcf_folder} {output} {threads} {input.ref}	
		"""

# N=$1		# number of splits 
#BED=$2		# path to BED file with start and end coordinates for each contig
#ID=$3		# path to sample id list
#In=$4		# path to input g.vcf.gz folder, including idx files. (index)
#O=$5		# path to output folder
#Tr=$6		# number of threads to use
#Ref=$7		# path to reference fasta

