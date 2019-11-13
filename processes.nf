process 'Alignment_BWA' {
    tag "$read_id"

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

    input:
        tuple path(read_id_bam), path(_metrics)

    output:
        tuple path(read_id_bam), file("${read_id_bam}.bai")

    script:
    """
    samtools index $read_id_bam
    """
}

process 'GenBQSR' {
    tag "$read_id_bam"

    input:
        tuple path(ref_dir), ref_filename
        tuple path(mills_dir), ks_mills
        // path ks_omni
        tuple path(read_id_bam), path(read_id_bai)

    output:
        file("${read_id_bam.baseName}.recal")

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

    input:
        tuple path(ref_dir), ref_filename
        path ks_mills
        // path ks_omni
        tuple path(read_id_bam), path(_metrics)

    output:
        file("${read_id_bam.baseName}.recal")

    script:
    java_mem = "-Xmx" + task.memory.toGiga() + "G"
    """
    gatk --java-options "-Xmx$JAVA_MEM -Djava.io.tmpdir=/SCRATCH/tmp" ApplyBQSR \
    -R $ref_dir/$ref_filename \
    -I $read_id_bam \
    -bqsr $input.recal \
    -O $output
    """
}