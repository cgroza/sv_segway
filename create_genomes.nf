params.vcf = "variants.vcf.gz"
params.vcf_lift = params.vcf
params.vcf2diploid = "vcf2diploid.jar"
params.fasta = "genome.fa"

vcf_ch = Channel.fromPath(params.vcf)
vcf_ch2 = Channel.fromPath(params.vcf)
vcf_lift_ch = Channel.fromPath(params.vcf_lift)
fasta_ch = Channel.fromPath(params.fasta)

process vcf_sample_names {
    input:
    file vcf_file from vcf_ch

    output:
    stdout vcf_samples_ch

    script:
    """
    bcftools query -l ${vcf_file}
    """
}

samples_ch = vcf_samples_ch.splitCsv().map{row -> "${row[0]}"}.combine(vcf_ch2).combine(vcf_lift_ch).combine(fasta_ch).view()

process vcf2diploid {
    time '4h'
    memory '10GB'
    cpus 2

    input:
    set val(sample_name), file(vcf_file), file(vcf_lift_file), file(fasta) from samples_ch

    output:
    file "${sample_name}_paternal" into paternal_genomes_ch
    file "${sample_name}_maternal" into maternal_genomes_ch

    script:
    """
    mkdir ${sample_name}_paternal
    mkdir ${sample_name}_maternal

    java -jar ${params.vcf2diploid} -id ${sample_name} -chr ${fasta} -vcf ${vcf_file}

    cat *paternal*.fa  > ${sample_name}_paternal/${sample_name}_paternal.fa
    mv paternal.chain ${sample_name}_paternal

    cat *maternal*.fa > ${sample_name}_maternal/${sample_name}_maternal.fa
    mv maternal.chain ${sample_name}_maternal

    CrossMap.py vcf ${sample_name}_paternal/paternal.chain ${vcf_lift_file} ${sample_name}_paternal/${sample_name}_paternal.fa ${sample_name}_paternal/${vcf_file.getName()} --compress
    CrossMap.py vcf ${sample_name}_maternal/maternal.chain ${vcf_lift_file} ${sample_name}_maternal/${sample_name}_maternal.fa ${sample_name}_maternal/${vcf_file.getName()} --compress

    bgzip ${sample_name}_paternal/${sample_name}_paternal.fa
    bgzip ${sample_name}_maternal/${sample_name}_maternal.fa

    rm *.map
    rm *.fa
    """
}

genomes_ch = paternal_genomes_ch.concat(maternal_genomes_ch)

process bwa_index_genome {
    time '12h'
    memory '16GB'
    cpus 8
    publishDir "genomes"

    input:
    file genome from genomes_ch

    output:
    file "${genome.getName()}"  into indexed_genomes_ch

    script:
    """
    bwa index ${genome}/${genome}.fa.gz
    """
}
