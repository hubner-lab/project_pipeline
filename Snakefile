import glob 
import re

configfile: "config.json"

#onerror:
#		shell(
#onsuccess:
#		shell(

# from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider(access_key_id="", secret_access_key="")

# TODO add resources 

def split_f(string):
	return '/'.join(map(str,string.split('/')[-2:]))[:-16]; # get the last 2 paths the remove the "_R{1|2}_001.fastq.gz" ending  

files=list(map(split_f,glob.glob("{}/*/*.fastq.gz".format(config["folder"])))) # very ugly but works

rule all:
	input:
		expand("VCF/{files}.list",files=files[:5]),

rule trim:
        input:
                R1=expand("{folder}/{{samples}}/{{file}}_R1_001.fastq.gz",folder = config["folder"]),
                R2=expand("{folder}/{{samples}}/{{file}}_R2_001.fastq.gz",folder = config["folder"]),
        output:
                R1="trimmed/{samples}/{file}_R1_001.fastq.gz",
                R2="trimmed/{samples}/{file}_R2_001.fastq.gz",
        threads: 
                workflow.cores * config["cores_per"] 
        message: 
                "trimming {wildcards.samples}/{wildcards.file}",
	benchmark:
		"benchmarks/trim/{samples}/{file}.tsv",
	log:
		"logs/trim/{samples}/{file}.log",
	shell:
		"""
		        fastp --detect_adapter_for_pe --overrepresentation_analysis --correction -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --json /dev/null --html /dev/null \
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
		gatk="/mnt/data/Tools/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar",
	output:
		bam=protected("mapping/{samples}/{file}.MarkedDup.bam"),
		txt=protected("mapping/{samples}/{file}-marked_dup_metrics.txts"),
	threads: 
		workflow.cores * config["cores_per"] 
	message: 
		"MarkDuplicates {wildcards.file}"
	benchmark:
		"benchmarks/mark_dup/{samples}/{file}.tsv"
	log:
		"logs/mark_dup/{samples}/{file}.log"
	shell:
		# java -jar {params.gatk} MarkDuplicates -I {input.bam} -O {output.bam} -M {output.txt} --VALIDATION_STRINGENCY LENIENT\
		# > {log} 2>&1
				#--spark-master local[{threads}]\
		"""
			java -jar {params.gatk} MarkDuplicatesSpark\
					-I {input.bam}\
					-O {output.bam}\
					-M {output.txt}\
		 			--VALIDATION_STRINGENCY LENIENT\
			> {log} 2>&1
		"""

# TODO: try the checkpoint feature
		# Unstable --  java -Xmx50G -jar {params.gatk} HaplotypeCallerSpark --spark-master local[{threads}] --tmp-dir tmp  -R {params.fasta} -I {input.MarkDup} -L chr6H -O {output} -ERC GVCF
checkpoint haplotypecaller: 
	input:
		MarkDup="mapping/{samples}/{file}.MarkedDup.bam",
		MarkDupIndex="mapping/{samples}/{file}.MarkedDup.bam.bai",
		fasta=config["fasta"],	
	output:
		dynamic("VCF/{samples}/{file}.g.vcf.{split}",
	params:
		gatk="/mnt/data/tools/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar",
		samples=config["splits"],
		batch=config["batch_size"],
	threads: 
		workflow.cores * config["cores_per"] 
	message:
		"running HaplotypeCaller on {input.fasta}"
	benchmark:
		"benchmarks/haplotypecaller/{samples}/{file}_{split}.tsv"
	shell:
		""" 
			source n_batch.sh

			N=$(( {threads} ))
			echo $N
			open_sem $N

			chrs=($(samtools view -H {input.MarkDup} | grep -oP "(?<=SN:)([^\s]*)"))

			for chr in "${{chrs[@]}}" 
			do
				len="$(samtools view -H {input.MarkDup} | grep "$chr" | grep -oP "(?<=LN:).*" )"

				l_sample=$(( len / {params.samples} ))

				f_break=0
				current=1
 				while true 
				do 
					next=$(( current + l_sample ))

					if [[ "$next" -ge "$len" ]] 
					then
						next=$len	
						f_break=1
					fi
					
					if [[ ! -f "VCF/{wildcards.samples}/{wildcards.file}.g.vcf.$chr.$current.$next" ]]
					then
						run_with_lock java -jar {params.gatk} HaplotypeCaller\*/
						-R {input.fasta}\
						-I {input.MarkDup}\
						-L "$chr:$current-$next"\
						-O "VCF/{wildcards.samples}/{wildcards.file}.g.vcf.$chr.$current.$next" \
						-ERC GVCF\
						> /dev/null 2>&1
					fi

					[[ "$f_break" -eq "1" ]] && break

					current=$next
				done
			done
		"""

#/*def list_files_input(wildcards):*/
#    /*checkpoint_output = checkpoints.haplotypecaller.get(**wildcards).output[0]*/

rule list_files:
	input:
		list_files_input
	output:
		"VCF/{samples}/{file}.list",
	threads: 
		workflow.cores * config["cores_per"] ,
	shell:
		"""
			echo {input} > {output}
		"""


#bcftools merge {output}* -g {input.fasta} -O v -o {output}  # bgzip all files

#c=$((c+1))
#c=$((c%N))
#[[ "$c" -eq "0" ]] && echo "waiting" && wait

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
#rule trim_test:
	#input:
		#R1=expand("{folder}/{{samples}}/{{file}}_R1_001.fastq.gz",folder = config["folder"]),
		#R2=expand("{folder}/{{samples}}/{{file}}_R2_001.fastq.gz",folder = config["folder"]),
	#output:
		#R1="trimmed/{samples}/{file}_R1_001.fastq.gz",
		#R2="trimmed/{samples}/{file}_R2_001.fastq.gz",
	#params:
		#trimmomatic="/home/hubner/Trimmomatic-0.36/trimmomatic-0.36.jar",
		#adapters="/home/hubner/Trimmomatic-0.36/adapters/AllAdapt.fa:2:30:10",
		#leadin=3,
		#trailing=3,
		#slidingwindow="4:10",
		#minlen=36,
	#threads: 
		#workflow.cores * config["cores_per"] 
	#message: 
		#"trimming {wildcards.samples}/{wildcards.file}"
	#benchmark:
		#"benchmarks/trim/{samples}/{file}.tsv"
	#log:
		#"logs/trim/{samples}/{file}.log"
	#shell:
		#"""
			#java -jar {params.trimmomatic} PE\
					#-threads {threads}\
					#{input.R1} {input.R2}\
					#{output.R1} /dev/null\
					#{output.R2} /dev/null\
					#ILLUMINACLIP:{params.adapters}\
					#LEADIN:{params.leadin}\
					#TRAILING:{params.trailing}\
					#SLIDINGWINDOW:{params.slidingwindow}\
					#MINLEN:{params.minlen}\
					#> {log} 2>&1	
	       #
	       #def aggregate_haplotypecaller(wildcards):
#    checkpoint_output = checkpoints.haplotypecaller.get(**wildcards).output[0]
#    print(checkpoint_output)
#    return file_names
#"""
