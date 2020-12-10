/*nextflow.preview.dsl=2*/
/*params.genome     = "$baseDir/data/genome.fa"*/
/*params.variants   = "$baseDir/data/known_variants.vcf.gz"*/
/*params.denylist   = "$baseDir/data/denylist.bed" */
/*params.reads      = "$baseDir/data/reads/rep1_{1,2}.fq.gz"*/
/*params.reads = "input/*_R{1,2}_001.fastq.gz"*/

params.threads = 100

println file('s3://lab.migal.bucket.us/ref/*.dict').text



/*Channel*/
    /*.fromFilePairs( params.reads )*/
    /*.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }*/
    /*.set { read_pairs_ch }*/

/*process trim {*/

    /*input:*/
    /*val threads from params.threads*/
    /*tuple val(pair_id), path(reads) from read_pairs_ch*/

    /*output:*/
    /*val "${reads.baseName}.out" into exp_channel*/

    /*script:*/
    /*trim_1="${reads[0]}.out"*/
    /*trim_2="${reads[1]}.out"*/
    /*"""*/
    /*fastp \*/
        /*--thread ${threads} \*/
        /*--detect_adapter_for_pe \*/
        /*--overrepresentation_analysis \*/
        /*--correction \*/
        /*-i ${reads[0]} \*/
        /*-I ${reads[1]} \*/
        /*-o ${trim_1} \*/
        /*-O ${trim_2} \*/
        /*--json /dev/null \*/
        /*--html /dev/null \*/
	/*> ${log} 2>&1	*/

    /*"""*/
/*}*/


/*[>[>read_pairs_ch.subscribe { println it }<]<]*/


/*process splitLetters {*/

    /*output:*/
    /*file 'chunk_*' into letters*/

    /*'''*/
    /*printf 'Hola' | split -b 1 - chunk_*/
    /*'''*/
/*}*/

/*letters*/
    /*.flatMap()*/
    /*.subscribe { println "File: ${it.name} => ${it.text}" }*/

