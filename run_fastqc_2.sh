#!/bin/tcsh

foreach f(*trim.fq.gz)

        fastqc -q $f --outdir=./FASTQC

end

