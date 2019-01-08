
all: HG02891.summary.txt Loschbour.hg19_1000g.summary.txt

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
	../src/rohan   --auto autosomes  --err  ../illuminaProf/error.prof  -t  4  -o  HG02891  human_g1k_v37.masked.fa HG02891.21.bam	


Loschbour.hg19_1000g.5p.prof: Loschbour.hg19_1000g.chr21.bam.bai
	samtools view -F1 -h -q 30 Loschbour.hg19_1000g.chr21.bam| awk '{if( substr($$1,1,1)=="@"  || length($$10)>45){ print $$0; }}'  | head -n 1000000 | samtools view -bS - | ../bam2prof/bam2prof -5p Loschbour.hg19_1000g.5p.prof -3p Loschbour.hg19_1000g.3p.prof -minq 30  -length 40  - 


Loschbour.hg19_1000g.chr21.bam.bai:
	wget -O Loschbour.hg19_1000g.chr21.bam     https://www.dropbox.com/s/1n24l8jijh7a08m/Loschbour.hg19_1000g.chr21.bam?dl=0
	wget -O Loschbour.hg19_1000g.chr21.bam.bai https://www.dropbox.com/s/nhdvjcr1no71re3/Loschbour.hg19_1000g.chr21.bam.bai?dl=0


Loschbour.hg19_1000g.summary.txt: Loschbour.hg19_1000g.chr21.bam.bai Loschbour.hg19_1000g.5p.prof Loschbour.hg19_1000g.chr21.bam.bai  human_g1k_v37.masked.fa.fai
	../src/rohan --rohmu 5e-5 --deam5p Loschbour.hg19_1000g.5p.prof --deam3p Loschbour.hg19_1000g.3p.prof --out Loschbour.hg19_1000g --auto autosomes  --err ../illuminaProf/error.prof human_g1k_v37.masked.fa Loschbour.hg19_1000g.chr21.bam
