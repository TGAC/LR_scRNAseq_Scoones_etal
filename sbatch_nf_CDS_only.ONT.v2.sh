#!/bin/sh
# run the pipeline on a very cut-down version of data
nextflow run ./single_tx_geneBodyCoverage.nf.ONT.v1.nf \
    -w $./work \
    -with-report -resume -with-trace -with-dag flowchart.html \
    -c ./GENANNO-579.nf.ONT.config
        

