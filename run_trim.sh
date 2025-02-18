#!/bin/bash
f1=( $(ls *R1.fq.gz) )
f2=( $(ls *R2.fq.gz) )
g1=( $(ls *R1.fq.gz) )
g2=( $(ls *R2.fq.gz) )

for i in {0..5}
do
	echo ${f1[$i]}
	echo ${f2[$i]}
	java -Djava.lirary.path=/home/data/tapm/transcriptomics/bbmap/jni/ -ea -Xmx16g -Xms16g -cp /home/data/tapm/transcriptomics/bbmap/current/ jgi.BBDukF -Xmx16g in1=${f1[$i]} in2=${f2[$i]} out1=${g1[$i]/.fq.gz/.trim.fq.gz} out2=${g2[$i]/.fq.gz/.trim.fq.gz}  k=31 mink=5 hdist=1 ktrim=r ref=/home/data/tapm/transcriptomics/bbmap/resources/adapters.fa trimq=20 qtrim=rl stats=./${g1[$i]}.stats
done

