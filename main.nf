#!/usr/bin/env nextflow

// Simple Nextflow Script for Testing
// File: main.nf

nextflow.enable.dsl=2

// Workflow definition
workflow {
    def trait_channel = Channel.of(*params.traits)
    trait_channel.view { trait -> println("trait: ${trait}") }
    def pair_channel = trait_channel.collect().map { traits -> [traits[0], traits[1]] }
    pair_channel.view { pair -> println("pair: ${pair[0]} and ${pair[1]}") }
    MungeGWAS(trait_channel)
    Heritability(trait_channel)
    GeneticCorrelation(pair_channel)
    Clump(pair_channel) 
    LD_Matrix(pair_channel)
    MR(pair_channel)
}

// Munge GWAS to hm3
process MungeGWAS {
    input:
        val trait

    script:
    """
    if [ ! -f ${params.path_munge}${trait}.sumstats.gz ]; then
        conda run -n ldsc bash -c " \
        ${baseDir}/ldsc/munge_sumstats.py \
        --sumstats ${params.path_gwas}${trait}.txt.gz \
        --out ${params.path_munge}${trait} \
        --merge-alleles ${params.path_hm3}"
    fi
    """
}

// Heritability estimation via single-trait LDSC
process Heritability {
    input:
        val trait

    script:
    """
    if [ ! -f ${params.path_result}${trait}.log ]; then
        conda run -n ldsc bash -c "\
        ${baseDir}/ldsc/ldsc.py \
        --h2 ${params.path_munge}${trait}.sumstats.gz \
        --ref-ld-chr ${params.path_LD} \
        --w-ld-chr ${params.path_LD} \
        --out ${params.path_result}${trait}"
    fi
    """
}

// Genetic correlation estimation vis Cross-trait LDSC
process GeneticCorrelation {
    input:
        val pair

    script:
    """
    if [ ! -f ${params.path_result}${pair[0]}_AND_${pair[1]}.log ]; then
        conda run -n ldsc bash -c "\
        ${baseDir}/ldsc/ldsc.py \
        --rg ${params.path_munge}${pair[0]}.sumstats.gz,${params.path_munge}${pair[1]}.sumstats.gz \
        --ref-ld-chr ${params.path_LD} \
        --w-ld-chr ${params.path_LD} \
        --out ${params.path_result}${pair[0]}_AND_${pair[1]}"
    fi
    """
}

// Clump significant and independent SNP for MR
process Clump {
    input:
        val pair

    script:
    """
    zcat ${params.path_gwas}${pair[0]}.txt.gz | awk 'NR==1 {print \$0} \$9<5e-8 {print \$0}' > ${params.path_clump}${pair[0]}.forClump

    if [ ! -f ${params.path_clump}${pair[0]}.iv ]; then
        plink --allow-no-sex \
            --clump ${params.path_clump}${pair[0]}.forClump \
            --bfile ${params.path_bfile} \
            --clump-kb 1000 --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.05 \
            --out ${params.path_clump}${pair[0]}

        awk 'NR>1 {print \$3}' ${params.path_clump}${pair[0]}.clumped > ${params.path_clump}${pair[0]}.iv
    fi
    """
}

// Make LD matrix for GSMR
process LD_Matrix {
    input:
        val pair

    script:
    """
    if [ ! -f ${params.path_clump}${pair[0]}.ld ]; then
        plink --allow-no-sex \
            --extract ${params.path_clump}${pair[0]}.iv \
            --bfile ${params.path_bfile} \
            --r2 square \
            --write-snplist \
            --out ${params.path_clump}${pair[0]}
    fi
    """
}

 process MR {
    input:
        val pair

    script:
    """
    conda run -n r421_mr bash -c "Rscript ${baseDir}/MR.r \
        ${params.path_gwas}${pair[0]}.txt.gz ${params.path_gwas}${pair[1]}.txt.gz \
        ${params.path_clump}${pair[0]}.iv \
        ${params.path_clump}${pair[0]}.ld \
        ${params.path_result}${pair[0]}_TO_${pair[1]}"
    """
}