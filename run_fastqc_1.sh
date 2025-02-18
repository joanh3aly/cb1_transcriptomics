#!/bin/tcsh

foreach f(*.fq.gz)

        fastqc -q $f --outdir=./FASTQC_1

end

