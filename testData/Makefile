
all: HG02891.summary.txt

HG02891.21.bam:
	samtools view -b ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02891/alignment/HG02891.mapped.ILLUMINA.bwa.GWD.low_coverage.20121211.bam 21 > HG02891.21.bam 

HG02891.21.bam.bai: HG02891.21.bam 
	samtools index HG02891.21.bam 

human_g1k_v37.masked.fa.gz:
	wget https://www.dropbox.com/s/4ayoefse8a5ri1n/human_g1k_v37.masked.fa.gz?dl=0 -O human_g1k_v37.masked.fa.gz

human_g1k_v37.masked.fa: human_g1k_v37.masked.fa.gz
	gunzip -f human_g1k_v37.masked.fa.gz

human_g1k_v37.masked.fa.fai: human_g1k_v37.masked.fa
	samtools faidx human_g1k_v37.masked.fa


HG02891.summary.txt: human_g1k_v37.masked.fa.fai HG02891.21.bam.bai
	../src/rohan  --nogl --auto autosomes  --err  ../illuminaProf/error.prof  -t  4  -o  HG02891  human_g1k_v37.masked.fa HG02891.21.bam	
