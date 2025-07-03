// standard run case:'nextflow run scqc_nf.sh -c PIP-3144.config -with-dag flowchart.png -with-report -resume'

nextflow.enable.dsl=2

process ListMolecule {
    // search pigeon classification output, scisoseq.annotated.info.csv, for a list of molecules mapped to the transcript, tx. 
    tag "$tx"
    errorStrategy 'retry'
    maxRetries 3 

    input: 
        val tx
    output:
        path "molecules.${tx}.list", emit: listfile
        path "ListMolecule.done"
        val tx, emit: tx

    script:

    """
    awk -v pattern="${tx}" -e '\$2==pattern {print \$1}' ${params.infocsv} | sort|uniq > molecules.${tx}.list && sleep 20 && \
    echo done > ListMolecule.done
    """
}

process PrintEmptyMoleculeList{
    tag "$tx"
    // errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' } 
    errorStrategy 'retry'
    maxRetries 3 
    publishDir "${params.emptylistout}", mode: 'copy', pattern: "molecules*list"

    input:
        path listfile
        val tx
    output:
        path "molecules.${tx}.list"
        path "PrintEmptyMoleculeList.done"
    when: 
        listfile.size() == 0

    script:
    """
    echo done > PrintEmptyMoleculeList.done && sleep 20
    """
}

process ExtractMoleculeBam {
    // make bam file containing only the listed molicues
    beforeScript 'source package 638df626-d658-40aa-80e5-14a275b7464b'
    tag "$tx"
    // errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' } 
    errorStrategy 'retry'
    maxRetries 3 
    publishDir "${params.moleculeout}", mode: 'copy', pattern: "molecules*list"

    input:
        path listfile
        val tx
    output:
        val tx, emit: tx
        path "molecules.${tx}.bam", emit: bam
        path "molecules.${tx}.bam.bai", emit: bambai
        path "molecules.${tx}.list"
        path "ExtractMoleculeBam.done"
    when: 
        listfile.size() > 0

    script:
    """
    samtools view --threads 4 --qname-file ${listfile} -b ${params.bamfile} | samtools sort --threads 4 -o molecules.${tx}.bam && \
        samtools index molecules.${tx}.bam && sleep 30 && \
        echo done > ExtractMoleculeBam.done

    """
}

process GeneBodyCoverage {
    publishDir "${params.rseqcout}", mode: 'copy', pattern: "*geneBodyCoverage.txt"
    tag "$tx"
    errorStrategy "retry"
    maxRetries 3 
    label "singularity_eimap"

    input: 
        path bam
        path bambai
        val tx
    output:
        path "${tx}.geneBodyCoverage.txt", emit: txt 
        path "GeneBodyCoverage.done"


    script:
    """
    geneBody_coverage.py -i ${bam} -r ${params.txdir}/${tx}.bed -o ${tx} && echo done > GeneBodyCoverage.done && sleep 30
    """
}

process CatGeneBodyCoverage {
    publishDir "${params.genebodyout}", mode: 'copy', pattern: "geneBodyCoverage*txt"
    input:
        path all_txt_files

    output:
        path "geneBodyCoverage.all.txt"
        path "geneBodyCoverage.sum.txt"
        path "geneBodyCoverage.norm.txt"
        path "CatGeneBodyCoverage.done"


    // normalisation: (i -min(dat))/(max(dat) - min(dat))
    script:
    """
        # cat ${all_txt_files} > GeneBodyCoverage.all.txt 
        if [ -f geneBodyCoverage.all.txt ]; then rm geneBodyCoverage.all.txt; fi && \
        for file in ${all_txt_files}; do cat \$file | grep -v Percentile >> geneBodyCoverage.all.txt; done && \
        awk '
        {
            for (i=2;i<=NF;i++) bin[i-1]=bin[i-1]+\$i
        }
        END{
            for(i=1; i<=100; i++){printf "%d\\t", bin[i]} 
            print "\\n"
        }' geneBodyCoverage.all.txt > geneBodyCoverage.sum.txt && \
        awk '
        {   for (i=1;i<=NF;i++) {
                if (min>\$i){min=\$i} 
                if (max<\$i){max=\$i} 
            }
            for (i=1;i<=NF;i++) {
                norm[i]=(\$i-min)/(max-min)
            }
        }
        END{
            for (i=1;i<=100;i++){
                printf "%f\\t", norm[i];
            }
            print "\\n";
        }' geneBodyCoverage.sum.txt > geneBodyCoverage.norm.txt && \
        echo done > CatGeneBodyCoverage.done

    """
    
}

workflow {
    
    // load list of transcripts
    tx_list = Channel.from(file(params.txlist).readLines())
    // tx_list.view()

    // search pigeon classification output, scisoseq.annotated.info.csv, for a list of molecules mapped to the transcript, tx. 
    ListMolecule(tx_list)

    // from pbmm2 mapped bamfile, scisoseq.mapped.bam, extract bams listed from ListMolecule()
    ExtractMoleculeBam(ListMolecule.out.listfile, ListMolecule.out.tx)

    PrintEmptyMoleculeList(ListMolecule.out.listfile, ListMolecule.out.tx)

    // compute coverage for each pair of bam and tx
    GeneBodyCoverage(ExtractMoleculeBam.out.bam, ExtractMoleculeBam.out.bambai, ExtractMoleculeBam.out.tx)

    // make coverage summary files
    CatGeneBodyCoverage(GeneBodyCoverage.out.txt.collect())

    


}
    







