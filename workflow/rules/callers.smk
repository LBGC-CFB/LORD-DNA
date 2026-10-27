rule sniffles:
    input: 
        f"{config['dirs']['outdir']}/alignment/{{aligner}}/{{sample}}.bam"
    output:
        f"{config['dirs']['outdir']}/sniffles/{{sample}}_{{aligner}}.vcf.gz"
    params:
        genome = config['refs']['genome']
    threads: config['threads']
    singularity: 
        config["singularity"]["sniffles"]
    shell:
        """
        sniffles \
            --input {input} \
            --reference {params.genome} \
            --long-dup-coverage 1.15 \
            --vcf {output} \
            --threads {threads} \
            --sample-id {wildcards.sample}
        """

rule cutesv:
    input:
        bam = f"{config['dirs']['outdir']}/alignment/{{aligner}}/{{sample}}.bam",
        bai = f"{config['dirs']['outdir']}/alignment/{{aligner}}/{{sample}}.bam.bai",
        ref = config["refs"]["genome"],
        fai = config["refs"]["genome"] + ".fai",
        bed = config["refs"]["bed"],
    output:
        vcf = f"{config['dirs']['outdir']}/cutesv/{{aligner}}/{{sample}}/cutesv.vcf",
        vcfgz = f"{config['dirs']['outdir']}/cutesv/{{aligner}}/{{sample}}/cutesv.vcf.gz",
        tbi = f"{config['dirs']['outdir']}/cutesv/{{aligner}}/{{sample}}/cutesv.vcf.gz.tbi",
    params:
        outdir = f"{config['dirs']['outdir']}/cutesv/{{aligner}}/{{sample}}",
        container = config["singularity"]["cutesv"],
        bcftools = config["singularity"]["bcftools"],
        ins_bias = config["tools"]["cutesv"]["max_cluster_bias_INS"],
        ins_ratio = config["tools"]["cutesv"]["diff_ratio_merging_INS"],
        del_bias = config["tools"]["cutesv"]["max_cluster_bias_DEL"],
        del_ratio = config["tools"]["cutesv"]["diff_ratio_merging_DEL"],
        min_support = config["tools"]["cutesv"]["min_support"],
        min_size = config["tools"]["cutesv"]["min_size"],
        max_size = config["tools"]["cutesv"]["max_size"],
        genotype = "--genotype" if config["tools"]["cutesv"].get("genotype", False) else "",
    log:
        f"{config['dirs']['outdir']}/cutesv/{{aligner}}/{{sample}}/cutesv.log"
    threads: config["threads"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir:q}
        workdir=$(mktemp -d {params.outdir:q}/work.XXXXXX)
        trap 'rm -rf -- "$workdir"' EXIT
        mkdir -p "$workdir/signatures"

        singularity exec {params.container:q} cuteSV \
            --max_cluster_bias_INS {params.ins_bias} \
            --diff_ratio_merging_INS {params.ins_ratio} \
            --max_cluster_bias_DEL {params.del_bias} \
            --diff_ratio_merging_DEL {params.del_ratio} \
            --threads {threads} --sample {wildcards.sample:q} \
            --min_support {params.min_support} \
            --min_size {params.min_size} --max_size {params.max_size} \
            -include_bed {input.bed:q} {params.genotype} \
            {input.bam:q} {input.ref:q} "$workdir/cutesv.raw.vcf" "$workdir/signatures" \
            > {log:q} 2>&1

        singularity exec {params.bcftools:q} bcftools sort \
            -T "$workdir/sort" -Oz -o {output.vcfgz:q} "$workdir/cutesv.raw.vcf" \
            >> {log:q} 2>&1
        singularity exec {params.bcftools:q} bcftools index -t -f {output.vcfgz:q} \
            >> {log:q} 2>&1
        singularity exec {params.bcftools:q} bcftools view -Ov -o {output.vcf:q} {output.vcfgz:q} \
            >> {log:q} 2>&1
        """

rule clair3: 
    input: f"{config['dirs']['outdir']}/alignment/{{aligner}}/{{sample}}.bam"
    output:
        f"{config['dirs']['outdir']}/clair3/{{aligner}}/{{sample}}/merge_output.vcf.gz"
    params:
        bed = config['refs']['bed'],
        genome = config['refs']['genome'],
        model = config['models']['clair3'],
        sample_dir = directory(f"{config['dirs']['outdir']}/clair3/{{aligner}}/{{sample}}")
    threads: config['threads']
    singularity: 
        config["singularity"]["clair3"]
    shell:
        """
        run_clair3.sh \
            --bed_fn={params.bed} \
            -b {input} \
            -f {params.genome} \
            -m {params.model} \
            -t {threads} \
            -p ont \
            -o {params.sample_dir} \
        """

rule deepvariant:
    input:
        f"{config['dirs']['outdir']}/alignment/{{aligner}}/{{sample}}.bam"
    output:
        f"{config['dirs']['outdir']}/deepvariant/{{aligner}}/{{sample}}/output.vcf.gz"
    params:
        bed=config['refs']['bed'],
        genome = config['refs']['genome'],
        output_dir = directory(f"{config['dirs']['outdir']}/deepvariant/{{aligner}}/{{sample}}"),
        model = config['models']['deepvariant'],
        container = config["singularity"]["deepvariant_gpu"],
        min_qual=config.get("deepvariant_min_qual", 1),
    threads: config['threads']
    resources:
        gpu_slot=1
    shell:
        """
        RAW_VCF="{params.output_dir}/output.raw.vcf.gz"

        singularity exec --nv {params.container} \
        /opt/deepvariant/bin/run_deepvariant \
            --model_type={params.model} \
            --ref={params.genome} \
            --reads={input} \
            --regions={params.bed} \
            --output_vcf="$RAW_VCF" \
            --output_gvcf={params.output_dir}/output.g.vcf.gz \
            --num_shards={threads}

        # 2) Filter out noise (QUAL=0)
        singularity exec --nv {params.container} bcftools view -i 'QUAL>={params.min_qual}' -Oz -o {output} "$RAW_VCF"

        singularity exec --nv {params.container} tabix -p vcf {output}
        """

rule pepper_deepvariant:
    input:
        f"{config['dirs']['outdir']}/alignment/{{aligner}}/{{sample}}.bam"
    output:
        f"{config['dirs']['outdir']}/pepper_deepvariant/{{aligner}}/{{sample}}/output.vcf.gz"
    params:
        genome = config['refs']['genome'],
        output_dir = directory(f"{config['dirs']['outdir']}/pepper_deepvariant/{{aligner}}/{{sample}}"),
        model = config['models']['pepper_deepvariant'],
        container = config["singularity"]["pepper-margin-deepvariant_gpu"],
    threads: config['threads']
    shell:
        """
        singularity exec --nv {params.container} \
        run_pepper_margin_deepvariant call_variant \
            --{params.model} \
            -f {params.genome} \
            -b {input} \
            -o {params.output_dir} \
            -p output \
            -g \
            -t {threads}
        """


rule spip:
    input:
        vcf = snv_vcf_gz
    output:
        csv = f"{config['dirs']['outdir']}/spip/{{caller}}/{{sample}}_{{aligner}}.csv"
    params:
        threads = config['threads'],
        maxlines = config['tools']['spip']['maxlines'],
        transcriptome = config['tools']['spip']['transcriptome'],
        container = config["singularity"]["spip"]
    shell:
        r"""
        set +e
        mkdir -p "$(dirname "{output.csv}")"

        tmp_vcf="{output.csv}.tmp.vcf"
        tmp_out="{output.csv}.tmp"

        gunzip -c "{input.vcf}" > "$tmp_vcf"

        singularity exec --bind "{params.transcriptome}:/opt/SPiP/RefFiles/transcriptome_hg38.RData" "{params.container}" \
          Rscript /opt/SPiP/SPiPv2.1_main.r \
            --transcriptome /opt/SPiP/RefFiles/transcriptome_hg38.RData \
            -I "$tmp_vcf" \
            -g hg38 \
            -O "$tmp_out" \
            -t "{params.threads}" \
            -l "{params.maxlines}"

        rc=$?

        if [ $rc -eq 0 ] && [ -s "$tmp_out" ]; then
            mv "$tmp_out" "{output.csv}"
        else
            echo "[WARN] SPiP failed for sample={wildcards.sample} caller={wildcards.caller} aligner={wildcards.aligner} (exit=$rc)" >&2

            printf 'CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tvarID\tInterpretation\tInterConfident\tSPiPscore\tstrand\tgNomen\tvarType\tntChange\tExonInfo\texonSize\ttranscript\tgene\tNearestSS\tDistSS\tRegType\tSPiCEproba\tSPiCEinter_2thr\tdeltaMES\tBP\tmutInPBarea\tdeltaESRscore\tposCryptMut\tsstypeCryptMut\tprobaCryptMut\tclassProbaCryptMut\tnearestSStoCrypt\tnearestPosSStoCrypt\tnearestDistSStoCrypt\tposCryptWT\tprobaCryptWT\tclassProbaCryptWT\tposSSPhysio\tprobaSSPhysio\tclassProbaSSPhysio\tprobaSSPhysioMut\tclassProbaSSPhysioMut\n' > "{output.csv}"
            rm -f "$tmp_out"
        fi

        rm -f "$tmp_vcf"
        """
    