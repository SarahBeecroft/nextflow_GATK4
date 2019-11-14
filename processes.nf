process 'Alignment_BWA' {
    tag "$read_id"
    label 'BWA'

    input:
        tuple path(ref_dir), ref_filename
        val rg_info
        tuple read_id, path(read1), path(read2)

    output:
        path("${read_id}.bwa")

    script:
    """
    bwa mem -M -t $task.cpus -R "$rg_info" $ref_dir/$ref_filename $read1 $read2 > ${read_id}.bwa
    """
}

process 'Alignment_Samtools' {
    tag "$read_id_bwa"
    label 'Samtools'

    input:
        path read_id_bwa

    output:
        file("${read_id_bwa.baseName}.bam")

    script:
    """
    samtools view -bSu $read_id_bwa | samtools sort - > ${read_id_bwa.baseName}.bam
    """
}

process 'MarkDuplicates' {
    tag "$read_id_bam"
    label 'Picard'

    input:
        path read_id_bam

    output:
        tuple file("${read_id_bam.baseName}.markedDups.bam"), \
              file("${read_id_bam.baseName}.metrics.txt")

    script:
    java_mem = "-Xmx" + task.memory.toGiga() + "G"
    """
    java $java_mem -jar -Djava.io.tmpdir=/tmp /usr/picard/picard.jar \
        MarkDuplicates \
        I=$read_id_bam \
        O=${read_id_bam.baseName}.markedDups.bam \
        M=${read_id_bam.baseName}.metrics.txt
    """
}

process 'BamOrder' {
    tag "$read_id_bam"
    label 'Picard'

    input:
        tuple path(ref_dir), ref_filename
        tuple path(read_id_bam), path(_metrics)

    output:
        file("${read_id_bam.baseName}.sorted.bam")

    script:
    java_mem = "-Xmx" + task.memory.toGiga() + "G"
    """
    java $java_mem -jar -Djava.io.tmpdir=/tmp /usr/picard/picard.jar \
        ReorderSam \
        I=$read_id_bam \
        R=$ref_dir/$ref_filename \
        O=${read_id_bam.baseName}.sorted.bam
    """
}

process 'IndexBAM' {
    tag "$read_id_bam"
    label 'Samtools'

    input:
        // tuple path(read_id_bam), path(_metrics)
        path(read_id_bam)

    output:
        tuple path(read_id_bam), file("${read_id_bam}.bai")

    script:
    """
    samtools index $read_id_bam
    """
}

process 'GenBQSR' {
    tag "$read_id_bam"
    label 'GATK'

    input:
        tuple path(ref_dir), ref_filename
        tuple path(mills_dir), ks_mills
        // path ks_omni
        tuple path(read_id_bam), path(read_id_bai)

    output:
        tuple path(read_id_bam), path(read_id_bai), file("${read_id_bam.baseName}.recal")

    script:
    java_mem = "-Xmx" + task.memory.toGiga() + "G"
    """
    gatk --java-options "$java_mem -Djava.io.tmpdir=/tmp" BaseRecalibrator \
        -R $ref_dir/$ref_filename \
        -I $read_id_bam \
        --known-sites $mills_dir/$ks_mills \
        -O ${read_id_bam.baseName}.recal
    """
}

process 'ApplyBQSR' {
    tag "$read_id_bam"
    label 'GATK'

    input:
        tuple path(ref_dir), ref_filename
        tuple path(read_id_bam), path(read_id_bai), path(bqsr_recal)

    output:
        tuple path("${read_id_bam.baseName}.BQSR.bam"), path(read_id_bai)

    script:
    java_mem = "-Xmx" + task.memory.toGiga() + "G"
    """
    gatk --java-options "$java_mem -Djava.io.tmpdir=/tmp" ApplyBQSR \
        -R $ref_dir/$ref_filename \
        -I $read_id_bam \
        -bqsr $bqsr_recal \
        -O ${read_id_bam.baseName}.BQSR.bam
    """
}

process 'HapCaller' {
    tag "$read_id_bam"
    label 'GATK'

    input:
        tuple path(ref_dir), ref_filename
        tuple path(read_id_bam), path(read_id_bai)

    output:
        file("${read_id_bam.baseName}.hapCall")

    script:
    java_mem = "-Xmx" + task.memory.toGiga() + "G"
    """
    gatk --java-options "$java_mem -Djava.io.tmpdir=/tmp" HaplotypeCaller \
        -R $ref_dir/$ref_filename \
        -I $read_id_bam \
        -O ${read_id_bam.baseName}.hapCall \
        -ERC GVCF
    """
}

process 'CombineGVCFs' {
    label 'GATK'

    input:
        tuple path(ref_dir), ref_filename
        path gatheredGVCFs

    output:
        file("MergedGVCF.g.vcf")

    script:
    java_mem = "-Xmx" + task.memory.toGiga() + "G"
    // gatk needs '--variant ' before each filename
    variants = gatheredGVCFs.collect({ "--variant " + it }).join(' ')
    """
    gatk --java-options $java_mem CombineGVCFs \
        -R $ref_dir/$ref_filename \
        $variants \
        -O MergedGVCF.g.vcf
    """
}

process 'GenotypeGVCFs' {
    label 'GATK'

    input:
        tuple path(ref_dir), ref_filename
        path(merged_gvcf)

    output:
        file("jointGenotyped.raw.g.vcf")

    script:
    java_mem = "-Xmx" + task.memory.toGiga() + "G"

    """
    gatk --java-options "$java_mem -Djava.io.tmpdir=/tmp" \
        GenotypeGVCFs \
        -R $ref_dir/$ref_filename \
        -V $merged_gvcf \
        -O jointGenotyped.raw.g.vcf
    """
}

// process 'Variant_recal_SNP' {
//     label 'GATK'
//     // need to add in correct references ie dbsnp
//     input:
//         tuple path(ref_dir), ref_filename
//         tuple path(mills_dir), ks_mills
//         path(genotyped_vcf)

//     output:
//        tuple file("cohort_SNP.recal"), file("cohort_SNP.tranches")

//     script:
//     java_mem = "-Xmx" + task.memory.toGiga() + "G"
//     // --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
//     // --resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI \
//     // --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
//     """
//     gatk --java-options "$java_mem -Djava.io.tmpdir=/tmp" \
//         VariantRecalibrator \
//         -R $ref_dir/$ref_filename \
//         -V $genotyped_vcf \
//         --resource:1000G,known=false,training=true,truth=false,prior=10.0 $mills_dir/$ks_mills \
//         --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $mills_dir/$ks_mills \
//         -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
//         -mode SNP \
//         -O cohort_SNP.recal \
//         --tranches-file cohort_SNP.tranches
//     """
// }

// process 'Variant_recal_indel' {
//     label 'GATK'
// // need to add in correct references ie dbsnp
//     input:
//         tuple path(ref_dir), ref_filename
//         file("jointGenotyped.raw.g.vcf.gz")
//         tuple path(mills_dir), ks_mills

//     output:
//         tuple file("cohort_indel.recal"), file("cohort_indel.tranches")

//     script:
//     java_mem = "-Xmx" + task.memory.toGiga() + "G"

//     """
//         gatk --java-options "$java_mem -Djava.io.tmpdir=/tmp" \
//             VariantRecalibrator \
//             -R $ref_dir/$ref_filename \
//             --max-gaussians 4 \
//             --resource:mills,known=false,training=true,truth=true,prior=12 $MILLS \
//             -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
//             -mode INDEL \
//             -V jointGenotyped.raw.g.vcf.gz \
//             -O cohort_indel.recal \
//             --tranches-file cohort_indel.tranches

//     """
// }

// process 'Apply_VQSR_snp' {
//     label 'GATK'
// // need to add in correct references ie dbsnp
//     input:
//         tuple path(ref_dir), ref_filename
//         file("jointGenotyped.raw.g.vcf.gz")
//         tuple file("cohort_SNP.recal"), file("cohort_SNP.tranches")

//     output:
//         file("jointGenotyped.recalSNP.vcf")

//     script:
//     java_mem = "-Xmx" + task.memory.toGiga() + "G"

//     """

//         gatk --java-options "$java_mem -Djava.io.tmpdir=/tmp" \
//             ApplyVQSR \
//             -R $ref_dir/$ref_filename \
//             -V jointGenotyped.raw.g.vcf.gz \
//             -O jointGenotyped.recalSNP.vcf \
//             --recal-file cohort_SNP.recal \
//             --tranches-file cohort_SNP.tranches \
//             -mode SNP \
//             --ts-filter-level 99.0

//     """
// }

// process 'Apply_VQSR_indel = {
// ' {
//     label 'GATK'
// // need to add in correct references ie dbsnp
//     input:
//         tuple path(ref_dir), ref_filename
//         file("jointGenotyped.raw.g.vcf.gz")
//         tuple file("cohort_SNP.recal"), file("cohort_SNP.tranches")

//     output:
//         file("jointGenotyped.recalSNP.vcf")

//     script:
//     java_mem = "-Xmx" + task.memory.toGiga() + "G"

//     """

//          gatk --java-options "$java_mem -Djava.io.tmpdir=/tmp" ApplyVQSR
//      -R $REF
//      -V cohort.jointGenotyped.recalSNP.vcf
//      --tranches-file cohort.indel.tranches
//      --recal-file cohort.indel.recal
//      -O cohort.jointGenotyped.recalSNPIndel.vcf
//      --ts-filter-level 99.0
//      -mode INDEL

//     """
// }
// Apply_VQSR_indel = {
// exec """

// """
// }

process 'GQ_filter' {
    label 'GATK'

    input:
        tuple path(ref_dir), ref_filename
        path(genotyped_vcf)

    output:
        file("filtered.vcf")

    script:
    java_mem = "-Xmx" + task.memory.toGiga() + "G"
    """
    gatk --java-options "$java_mem -Djava.io.tmpdir=/tmp" \
        VariantFiltration \
        -R $ref_dir/$ref_filename \
        -V $genotyped_vcf \
        --genotype-filter-expression "GQ<20" \
        --genotype-filter-name "lowGQ" \
        -O filtered.vcf
    """
}
