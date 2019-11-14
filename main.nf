nextflow.preview.dsl=2

/*
 * Define the default parameters
 */

params.ref_dir = "$baseDir/data/"
params.ref_name = "hg19_all_chromosomes.sorted.fa"
// params.readGroup_info = "@RG\\tID:000000000-BPC6F.1\\tPL:illumina\\tPU:None\\tLB:None\\tSM:103"
// params.reads = "$baseDir/test_data/test-datasets/testdata/1_S103_L001_R{1,2}_001.fastq.gz"
params.bam = "$baseDir/data/*.bam"
params.bundle_dir = "$baseDir/data/hg19/"
params.mills = "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
//params.omni = "$baseDir/data/hg19/1000G_omni2.5.hg19.sites.vcf.gz"

log.info """\
NEXTFLOW Exome variant calling
================================
reference dir     : $params.ref_dir
reference filename: $params.ref_name
bam file          : $params.bam
bundle dir        : $params.bundle_dir
known sites mills : $params.mills
"""
// reads             : $params.reads
// known sites omni  : $params.omni

include './processes' params(params)

workflow {
    input_ref = Channel.value( [file(params.ref_dir), params.ref_name] )
    // input_rg_info = params.readGroup_info
    input_ks_mills = Channel.value([file(params.bundle_dir), params.mills])
    // input_ks_omni = Channel.value(file(params.omni))
    // reads_ch = Channel.fromFilePairs(params.reads, flat: true)
    bam_ch = Channel.fromPath(params.bam)

    // Alignment_BWA(input_ref, input_rg_info, reads_ch)
    // Alignment_Samtools(Alignment_BWA.out)
    // MarkDuplicates(Alignment_Samtools.out)
    // BamOrder(input_ref, MarkDuplicates.out)
    // IndexBAM(MarkDuplicates.out)
    IndexBAM(bam_ch)
    GenBQSR(input_ref, input_ks_mills, IndexBAM.out)
    ApplyBQSR(input_ref, GenBQSR.out)
    HapCaller(input_ref, ApplyBQSR.out)
    CombineGVCFs(input_ref, HapCaller.out.collect())
    GenotypeGVCFs(input_ref, CombineGVCFs.out)
    //Variant_recal_SNP(input_ref, input_ks_mills, GenotypeGVCFs.out)
    // Variant_recal_indel(input_ref, GenotypeGVCFs.out, input_ks_mills)
    // Apply_VQSR_snp(input_ref, GenotypeGVCFs.out, Variant_recal_SNP.out)
    // Apply_VQSR_indel()
    GQ_filter(input_ref, GenotypeGVCFs.out)
}
