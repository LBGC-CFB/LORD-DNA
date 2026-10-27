rule minimap2:
    input: f"{config['dirs']['outdir']}/dorado_demux/{{sample}}.fastq.gz"
    #input: f"{config['dirs']['fastqdir']}/{{sample}}.fastq.gz"
    output:
        sam = temp(f"{config['dirs']['outdir']}/alignment/minimap2/{{sample}}.sam")
    params:
        ref = config['refs']['genome']
    threads:
        config['threads']
    singularity:
        config["singularity"]["minimap2"]
    shell:
        """
        minimap2 \
            -ax map-ont \
            -t {threads} \
            -Y \
            {params.ref} \
            -o {output.sam} \
            {input}
        """


rule samtools_index:
    input: f"{config['dirs']['outdir']}/alignment/minimap2/{{sample}}.sam"
    output:
        bam = f"{config['dirs']['outdir']}/alignment/minimap2/{{sample}}.bam",
        bai = f"{config['dirs']['outdir']}/alignment/minimap2/{{sample}}.bam.bai"
    params:
        bed=config['refs']['bed'],
    threads:
        config['threads']
    singularity: 
        config["singularity"]["samtools"]
    shell:
        """
        samtools view -@ {threads} -b {input} | \
        samtools sort -@ {threads} -o {output.bam} && \
        samtools index -@ {threads} {output.bam}
        """
