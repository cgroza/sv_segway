params.dataset = "dataset.csv"

reads_genome_ch = Channel.fromPath(params.dataset).splitCsv(header:true)map{row ->
    [file(row.fastq), file(row.genome)]
}.view()

process bwa_align {
    cpus 8
    memory '16GB'
    module 'mugqic/bwa:mugqic/samtools'
    time '24h'

    input:
    set file(reads), file(genome) from reads_genome_ch

    output:
    file("${reads.getSimpleName()}_${genome.getSimpleName()}.bam")

    publishDir 'aligned'

    script:
    """
    bwa mem -t 8 -p ${genome}/${genome.getSimpleName()}.fa.gz ${reads} | samtools sort -@8 -OBAM -o ${reads.getSimpleName()}_${genome.getSimpleName()}.bam
    """
}
