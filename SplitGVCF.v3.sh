N=$1		# number of splits 
BED=$2		# path to BED file with start and end coordinates for each contig
ID=$3		# path to sample id list
In=$4		# path to input g.vcf.gz folder, including idx files. (index)
O=$5		# path to output folder
Tr=$6		# number of threads to use
Ref=$7		# path to reference fasta

mkdir -p $O"Maps/" $O"LOG/" $O"DB/" $O"VCF/" $O"Run/"

while read p ; do 
	chr=$(echo "$p" | cut -f1) ; 
	a=$(echo "$p" | cut -f3) ; 
	b=$(echo "$a / $N" | bc) ; 
	x=1 ; 
	n=1 ; 
	while ((x<=$a)) ; do 
		let c=$x+$b ; 
		while read q ; do 
			echo -e $q"\t""$In"$q.$chr.g.vcf.gz   
		done < $ID > $O"Maps/"$chr.map ; 
			
			echo "java -jar /mnt/data/Tools/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar GenomicsDBImport --genomicsdb-workspace-path $O"DB/"$chr.$x-$c -L $chr:$x-$c --sample-name-map $O"Maps"/$chr.map" >> $O"Run/"GenomicDB.par.sh
		
			echo "java -jar /mnt/data/Tools/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar GenotypeGVCFs -R $Ref -V gendb://$O"DB/"$chr.$x-$c/ -G StandardAnnotation -O $O"VCF/"$n.$chr.$x-$c.vcf.gz" >> $O"Run/"GenotypeGVCF.par.sh 

		let x=x+$b ; 
		let n=n+1 ; 
	done ; 
done < $BED  ; 

parallel -j $Tr :::: $O"Run/"GenomicDB.par.sh
parallel -j $Tr :::: $O"Run/"GenotypeGVCF.par.sh
