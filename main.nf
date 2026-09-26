#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process PREPARE_DONORS {
    tag "${donors.baseName}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd4d58047a5a9c36d1ab9f6b9ac3242d1ee123e773f96548a3ddf5a26fb15453/data'

    input:
    path(donors)

    output:
    path("${donors.baseName}_vcf.txt"), emit: donors_vcf

    script:
    """
    prepare_donors.py ${donors} ${donors.baseName}_vcf.txt
    """
}

process BCFTOOLS_EXTRACT_DONORS {
    tag "${vcf.getBaseName(2)}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a08933da914fc6b3650dbe842d7c2e4d95155d9a57a8c49c3756ad7248b7fc2c/data'

    input:
    path(vcf)
    path(donors_vcf)
    
    output:
    path("${vcf.getBaseName(2)}_donors.vcf.gz")

    script:
    """
    bcftools view ${vcf} \\
        --samples-file <(cut -d' ' -f1 ${donors_vcf}) \\
        --output-type z \\
        --output ${vcf.getBaseName(2)}_extracted.vcf.gz
    
    bcftools reheader ${vcf.getBaseName(2)}_extracted.vcf.gz \\
        --samples ${donors_vcf} \\
        --output ${vcf.getBaseName(2)}_donors.vcf.gz
    """
}

process BCFTOOLS_LIFTOVER {
    tag "${vcf.getBaseName(2)}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a08933da914fc6b3650dbe842d7c2e4d95155d9a57a8c49c3756ad7248b7fc2c/data'

    input:
    path(vcf)
    path(source_fasta)
    path(target_fasta)
    path(chain_file)

    output:
    path("${vcf.getBaseName(2)}_lift_sorted.vcf.gz"), emit: lifted
    path("${vcf.getBaseName(2)}_lift_rejected_sorted.vcf.gz"), emit: rejected

    script:
    """
    bcftools norm ${vcf} \\
        --check-ref s \\
        --do-not-normalize \\
        --fasta-ref ${source_fasta} \\
        --output-type u \\
        | bcftools +liftover \\
            --output-type u \\
            -- \\
            --src-fasta-ref ${source_fasta} \\
            --fasta-ref ${target_fasta} \\
            --chain ${chain_file} \\
            --write-reject \\
            --reject ${vcf.getBaseName(2)}_lift_rejected_unsorted.vcf.gz \\
            --reject-type z \\
        | bcftools sort \\
            --output-type z \\
            --output ${vcf.getBaseName(2)}_lift_sorted.vcf.gz

    bcftools sort ${vcf.getBaseName(2)}_lift_rejected_unsorted.vcf.gz \\
        --output-type z \\
        --output ${vcf.getBaseName(2)}_lift_rejected_sorted.vcf.gz
    """
}

process BCFTOOLS_VARIANT_QC {
    tag "${vcf.getBaseName(2)}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a08933da914fc6b3650dbe842d7c2e4d95155d9a57a8c49c3756ad7248b7fc2c/data'

    input:
    path(vcf)

    output:
    path("${vcf.getBaseName(2)}_filt.vcf.gz"), emit: vcf

    script:
    """
    bcftools view ${vcf} \\
        --min-ac ${params.min_mac}:minor \\
        --include 'F_MISSING <= ${params.max_missing}' \\
        --min-alleles 2 \\
        --max-alleles 2 \\
        --types snps \\
        --output-type z \\
        --output ${vcf.getBaseName(2)}_filt.vcf.gz
    """
}

process BCFTOOLS_INDEX {
    tag "${vcf.getBaseName(2)}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a08933da914fc6b3650dbe842d7c2e4d95155d9a57a8c49c3756ad7248b7fc2c/data'

    input:
    path(vcf)

    output:
    tuple path("${vcf}"), path("${vcf}.csi"), emit: vcf

    script:
    """
    bcftools index --csi ${vcf}
    """
}

process PLINK_MAKE_PFILE {
    tag "${vcf.getBaseName(2)}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fbbbd635adc17f214e69145009a0d1d0411c350b5e70eb14b5aa68d79a3fa1b/data'

    input:
    path(vcf)

    output:
    tuple path("${vcf.getBaseName(2)}.pgen"), path("${vcf.getBaseName(2)}.psam"), path("${vcf.getBaseName(2)}.pvar"), emit: pfile

    script:
    """
    plink2 \\
        --vcf ${vcf} \\
        --set-all-var-ids '@:#:\$r:\$a' \\
        --make-pgen \\
        --out ${vcf.getBaseName(2)}
    """
}

process PLINK_INDEP_PAIRWISE {
    tag "${pgen.baseName}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fbbbd635adc17f214e69145009a0d1d0411c350b5e70eb14b5aa68d79a3fa1b/data'

    input:
    tuple path(pgen), path(psam), path(pvar)
    path(ld_exclude_bed)

    output:
    tuple path(pgen), path(psam), path(pvar), path("${pgen.baseName}.prune.in"), emit: pfile

    script:
    """
    plink2 \\
        --pfile ${pgen.baseName} \\
        --exclude bed1 ${ld_exclude_bed} \\
        --bad-ld \\
        --indep-pairwise 200kb 0.2 \\
        --out ${pgen.baseName}
    """
}

process PLINK_PCA {
    tag "${pgen.baseName}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fbbbd635adc17f214e69145009a0d1d0411c350b5e70eb14b5aa68d79a3fa1b/data'

    publishDir "${params.outdir}/plink2", mode: 'copy', pattern: "*.eigenvec"
    input:
    tuple path(pgen), path(psam), path(pvar), path(pruned)

    output:
    path("${pgen.baseName}.eigenvec")

    script:
    """
    plink2 \
        --pfile ${pgen.baseName} \\
        --extract ${pruned} \\
        --bad-freqs \\
        --pca ${params.n_genotype_pcs} \\
        --out ${pgen.baseName}
    """
}

process PREPARE_GENE_REGIONS {
    tag "${gtf.baseName}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd4d58047a5a9c36d1ab9f6b9ac3242d1ee123e773f96548a3ddf5a26fb15453/data'
    publishDir "${params.outdir}/gene_regions", mode: 'copy'
    
    input:
    path(gtf)
    path(fai)

    output:
    path("gene_regions.tsv"), emit: regions

    script:
    """
    prepare_gene_regions.py ${gtf} ${fai} gene_regions.tsv ${params.cis_window}
    """
}

process PREPARE_COVARIATES {
    tag "${donors}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5c688ea7f743de8aa32394006d13bd87e1d5d03df3e1e3b0443907de1ac786c9/data'
    publishDir "${params.outdir}/covariates"

    input:
    path(donors)
    path(pca)

    output:
    path("covariates.tsv")

    script:
    """
    prepare_covariates.py ${donors} ${pca}
    """
}

process PREPARE_QTL_PHENOTYPES {
    tag "${adata.baseName}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/838f95344b08f1f01048ac394d5e0bc1186657d47f6e9a45398e3661f91c41fb/data'
    publishDir "${params.outdir}/phenotypes"

    input:
    path(adata)
    path(regions)

    output:
    path("*.phenotype.bed.gz"), emit: phenotypes
    path("*.offset.tsv"), emit: offsets

    script:
    """
    export NUMBA_CACHE_DIR="\$PWD/.numba_cache"
    mkdir -p "\$NUMBA_CACHE_DIR"

    prepare_qtl_phenotypes.py \\
        ${adata} \\
        ${regions} \\
        ${params.group_col}
    """
}

process JAXQTL_COMPUTE_PCS {
    tag "${phenotype_name}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/84da39d64268bd3e4e759a8de83d38408f3edf0da682eafb3e8197d3e0457dc3/data'
    publishDir "${params.outdir}/covariates", mode: 'copy', pattern: '*pcs.tsv'
    input:
    tuple val(phenotype_name), path(phenotype),  path(offset)
    path(covariates)

    output:
    tuple val(phenotype_name), path(phenotype), path(offset), path("${phenotype_name}.covariates_pcs.tsv")

    script:
    """
    jaxqtl compute-pcs \\
        --pheno ${phenotype} \\
        --covar ${covariates} \\
        --num-pcs ${params.n_expr_pcs} \\
        --out ${phenotype_name}.covariates_pcs.tsv
    """
}

process JAXQTL_CIS_GENES {
    tag "${phenotype_name}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/84da39d64268bd3e4e759a8de83d38408f3edf0da682eafb3e8197d3e0457dc3/data'
    publishDir "${params.outdir}/jaxqtl", mode: 'copy'

    input:
    tuple val(phenotype_name), path(phenotype),  path(offset), path(covariates)
    tuple path(vcf), path(vcf_index)

    output:
    tuple val(phenotype_name), path("${phenotype_name}.cis.score.spa.acat.parquet.gz"), emit: cis_genes

    script:
    """
    jaxqtl cis \\
        --vcf ${vcf} \\
        --pheno ${phenotype} \\
        --covar ${covariates} \\
        --offset ${offset} \\
        --model nb \\
        --test score \\
        --window ${params.cis_window} \\
        --tss-centered \\
        --min-gene-expr-pct ${params.min_gene_expr_pct} \\
        --one-hot \\
        --normalize-covar \\
        --spa \\
        --acat \\
        --out ${phenotype_name}
    """
}

process JAXQTL_NOMINAL {
    tag "${phenotype_name}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/84da39d64268bd3e4e759a8de83d38408f3edf0da682eafb3e8197d3e0457dc3/data'
    publishDir "${params.outdir}/jaxqtl", mode: 'copy'
   
    input:
    tuple val(phenotype_name), path(phenotype),  path(offset), path(covariates)
    tuple path(vcf), path(vcf_index)

    output:
    tuple val(phenotype_name), path("${phenotype_name}.nominal.score.parquet.gz"), emit: nominal

    script:
    """
    jaxqtl nominal \\
        --vcf ${vcf} \\
        --pheno ${phenotype} \\
        --covar ${covariates} \\
        --offset ${offset} \\
        --model nb \\
        --test wald \\
        --window ${params.cis_window} \\
        --tss-centered \\
        --min-gene-expr-pct ${params.min_gene_expr_pct} \\
        --one-hot \\
        --normalize-covar \\
        --out ${phenotype_name}
    """
}

workflow {
    if (!params.group_col) {
        error 'Missing required parameter: --group_col'
    }

    // Plain file inputs keep shared preparation outputs as reusable value channels.
    adata = file(params.adata, checkIfExists: true)
    donors = file(params.donors, checkIfExists: true)
    vcf = file(params.vcf, checkIfExists: true)
    source_fasta = file(params.source_fasta, checkIfExists: true)
    target_fasta = file(params.target_fasta, checkIfExists: true)
    target_fai = file("${params.target_fasta}.fai", checkIfExists: true)
    chain_file = file(params.chain_file, checkIfExists: true)
    gene_gtf = file(params.gene_gtf, checkIfExists: true)
    ld_exclude_bed = file(params.ld_exclude_bed, checkIfExists: true)

    PREPARE_DONORS(donors)
    BCFTOOLS_EXTRACT_DONORS(vcf, PREPARE_DONORS.out.donors_vcf)
    BCFTOOLS_LIFTOVER(
        BCFTOOLS_EXTRACT_DONORS.out,
        source_fasta,
        target_fasta,
        chain_file
    )
    BCFTOOLS_VARIANT_QC(BCFTOOLS_LIFTOVER.out.lifted)
    BCFTOOLS_INDEX(BCFTOOLS_VARIANT_QC.out.vcf)

    PLINK_MAKE_PFILE(BCFTOOLS_VARIANT_QC.out.vcf)
    PLINK_INDEP_PAIRWISE(PLINK_MAKE_PFILE.out.pfile, ld_exclude_bed)
    PLINK_PCA(PLINK_INDEP_PAIRWISE.out.pfile)
    PREPARE_COVARIATES(donors, PLINK_PCA.out)

    PREPARE_GENE_REGIONS(gene_gtf, target_fai)
    PREPARE_QTL_PHENOTYPES(adata, PREPARE_GENE_REGIONS.out.regions)

    ch_phenotypes = PREPARE_QTL_PHENOTYPES.out.phenotypes
        .flatten()
        .map { pheno ->
            def name = pheno.name.replaceFirst(/\.phenotype\.bed\.gz$/, '')
            tuple(name, pheno, pheno.resolveSibling("${name}.offset.tsv"))
        }

    JAXQTL_COMPUTE_PCS(ch_phenotypes, PREPARE_COVARIATES.out)
    JAXQTL_CIS_GENES(JAXQTL_COMPUTE_PCS.out, BCFTOOLS_INDEX.out.vcf)
    JAXQTL_NOMINAL(JAXQTL_COMPUTE_PCS.out, BCFTOOLS_INDEX.out.vcf)
}
