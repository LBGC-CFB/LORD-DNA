rule greenvaran:
    input:
        sv_vcf_gz
    output:
        f"{config['dirs']['outdir']}/greenvaran/{{svcaller}}/{{sample}}_{{aligner}}.vcf"
    params:
        db = config["databases"]["greenvaran"],
        dbschema = config["models"]["greendb"],
    threads: config["threads"]
    singularity:
        config["singularity"]["greenvaran"]
    shell:
        r"""
        mkdir -p $(dirname {output})

        /usr/GREEN-VARAN-1.1/workflow/bin/greenvaran sv \
            -i {input} \
            -o {output} \
            -d {params.db} \
            -s {params.dbschema} \
            -p 2

        sample_name=$(basename {output} | cut -d'_' -f1)
        sed -i "s/SAMPLE/$sample_name/g" {output}
        """
        
rule AnnotSV:
    input:
        f"{config['dirs']['outdir']}/greenvaran/{{svcaller}}/{{sample}}_{{aligner}}_final.vcf"
    output:
        f"{config['dirs']['outdir']}/AnnotSV/{{svcaller}}/{{sample}}_{{aligner}}.tsv"
    params:
        annotdir = config["databases"]["annotsv"],
        genome_version = config["refs"]["genome_version"],
    threads: config["threads"]
    singularity:
        config["singularity"]["annotsv"]
    resources:
        gpu_slot=1
    shell:
        r"""
        mkdir -p $(dirname {output})

        AnnotSV \
            -annotationsDir {params.annotdir} \
            -genomeBuild {params.genome_version} \
            -SVinputFile {input} \
            -outputFile {output}

        if [ ! -s {output} ]; then
            echo "AnnotSV output missing or empty. Converting input VCF to TSV format instead."
            awk -F '\t' 'BEGIN {{OFS="\t"}} /^#/ {{print}} /^[^#]/ {{print $1, $2, $4, $5, $8}}' {input} > {output}
        fi
        """

rule add_vaf_vcf_sv:
    input:
        f"{config['dirs']['outdir']}/greenvaran/{{svcaller}}/{{sample}}_{{aligner}}.vcf"
    output:
        f"{config['dirs']['outdir']}/greenvaran/{{svcaller}}/{{sample}}_{{aligner}}_final.vcf"
    params:
        container = config["singularity"]["python"],
        greendb = config["databases"]["greendb"],
        script = f"{workflow.basedir}/add_vaf_vcf.py",
    threads: config["threads"]
    shell:
        r"""
        singularity exec --nv {params.container} bash -c 'export TMPDIR=/mnt/csf/tmp && python3 {params.script} {input} {params.greendb} {output}'
        """

rule annotsv_to_excel:
    input:
        lambda wc: expand(
            f"{config['dirs']['outdir']}/AnnotSV/{wc.svcaller}/{{sample}}_{{aligner}}.tsv",
            sample=config["samples"],
            aligner=ALIGNERS,
        )
    output:
        filtered = f"{config['dirs']['outdir']}/filtered_final_sv_{{svcaller}}.xlsx",
        unfiltered = f"{config['dirs']['outdir']}/unfiltered_final_sv_{{svcaller}}.xlsx"
    params:
        container = config["singularity"]["python"],
        db = config["databases"]["LR_SV"],
        prom = config["refs"]["bed_prom"],
        occu = config["databases"]["LR_occu"],
        occu_sv = config["databases"]["LR_occu_sv"],
        script = f"{workflow.basedir}/annotSV_to_excel.py",
    threads: config["threads"]
    shell:
        r"""
        # openpyxl reopens temporary files; use local scratch instead of the SSHFS share.
        singularity exec --nv {params.container:q} bash -c '
          set -euo pipefail
          TMPDIR=$(mktemp -d /tmp/lr-excel.XXXXXX)
          export TMPDIR
          trap "rm -rf -- \"$TMPDIR\"" EXIT
          python3 "$@"
        ' bash {params.script:q} {input:q} {output.filtered:q} {output.unfiltered:q} {params.db:q} {params.prom:q} {params.occu:q} {params.occu_sv:q}
        """


rule greenvaran_snv:
    input:
        snv_vcf_gz
    output:
        f"{config['dirs']['outdir']}/greenvaran_snv/{{caller}}/{{sample}}_{{aligner}}.vcf"
    params:
        db = config['databases']['greenvaran'],
        dbschema = config['models']['greendb'],
        config = config['tools']['greenvaran']['conf'],
    threads: config['threads']
    singularity: 
        config["singularity"]["greenvaran"]
    shell:
        """
        mkdir -p $(dirname {output})

        /usr/GREEN-VARAN-1.1/workflow/bin/greenvaran smallvars \
            -i {input} \
            -o {output} \
            -d {params.db} \
            -s {params.dbschema} \
            --config {params.config}

        sample_name=$(basename {output} | cut -d'_' -f1)
        sed -i "s/SAMPLE/$sample_name/g" {output}
        """

rule SnpEff:
    input: f"{config['dirs']['outdir']}/pepper_deepvariant/{{aligner}}/{{sample}}/output.vcf.gz",
    output:
        f"{config['dirs']['outdir']}/SnpEff_SnpSift/{{sample}}_{{aligner}}.vcf"
    params:
        bed = config['refs']['bed']
    threads: config['threads']
    singularity: 
        config["singularity"]["snpeff"]
    shell:
        """
        java -Xmx8g -jar /opt/custflow/epi2meuser/conda/share/snpeff-5.1-2/snpEff.jar \
            -nodownload \
            -fi {params.bed} \
            GRCh38.p13 \
            {input} > \
            {output} \
        """

rule SnpSift_Clinvar:
    input: f"{config['dirs']['outdir']}/SnpEff_SnpSift/{{sample}}_{{aligner}}.vcf"
    output:
        f"{config['dirs']['outdir']}/SnpEff_SnpSift/{{sample}}_{{aligner}}_SnpSift_Clinvar.vcf"
    params:
        db = config['databases']['clinvar'],
    threads: config['threads']
    singularity:
        config["singularity"]["snpeff"]
    shell:
        """
        java -jar /opt/custflow/epi2meuser/conda/share/snpsift-5.1-0/SnpSift.jar \
            annotate \
            -v {params.db} \
            {input} > \
            {output}
        """

rule SnpSift_dbNSFP:
    input: f"{config['dirs']['outdir']}/SnpEff_SnpSift/{{sample}}_{{aligner}}_SnpSift_Clinvar.vcf"
    output:
        f"{config['dirs']['outdir']}/SnpEff_SnpSift/{{sample}}_{{aligner}}_SnpSift_Clinvar_dbNSFP.vcf"
    params:
        db = config['databases']['dbNSFP'],
    threads: config['threads']
    singularity: 
        config["singularity"]["snpeff"]
    shell:
        """
        java -jar /opt/custflow/epi2meuser/conda/share/snpsift-5.1-0/SnpSift.jar \
            dbnsfp \
            -v -db {params.db} \
            -n \
            -f Geuvadis_eQTL_target_gene \
            {input} > \
            {output}
        """

rule VEP_annotation:
    input: 
        f"{config['dirs']['outdir']}/greenvaran_snv/{{caller}}/{{sample}}_{{aligner}}.vcf"
    output: 
        f"{config['dirs']['outdir']}/vep/{{caller}}/{{sample}}_{{aligner}}_VEP_annotated.vcf"
    params:
        alphamissense = config['tools']['vep']['alphamissense'],
        assembly = config['refs']['genome_version'],
        container = config["singularity"]["vep"],
        genome = config['refs']['genome'],
        maxentscan = config['tools']['vep']['maxentscan'],
        revel = config['tools']['vep']['revel'],
        species = "homo_sapiens",
        spliceaisnv = config['tools']['vep']['spliceaisnv'],
        spliceaiindel = config['tools']['vep']['spliceaiindel'],
        vepdir = config['dirs']['vep'],
    threads: config['threads']
    shell:
        """
        mkdir -p $(dirname {output})

        singularity run --bind ~/vep_plugins:/plugins {params.container} \
            --verbose \
            --assembly {params.assembly} \
            -i {input} \
            --format vcf \
            -o {output} \
            --force_overwrite \
            --fork 4 \
            --cache \
            --dir {params.vepdir} \
            --offline \
            --fasta {params.genome} \
            --refseq \
            --plugin SpliceAI,snv={params.spliceaisnv},indel={params.spliceaiindel} \
            --plugin AlphaMissense,file={params.alphamissense} \
            --plugin REVEL,{params.revel} \
            --plugin MaxEntScan,{params.maxentscan} \
            --use_transcript_ref \
            --vcf \
            --sift b \
            --polyphen b \
            --regulatory \
            --vcf_info_field ANN \
            --hgvs \
            --hgvsg \
            --mane \
            --check_existing \
            --af \
            --af_gnomade \
            --af_gnomadg \
            --pick \
            --pick_order mane_select,canonical,appris,tsl,biotype,ccds,length
        """

rule split_vcf_spliceAI:
    input:
        f"{config['dirs']['outdir']}/vep/{{caller}}/{{sample}}_{{aligner}}_VEP_annotated.vcf"
    output:
        with_spliceAI=f"{config['dirs']['outdir']}/vep/{{caller}}/{{sample}}_{{aligner}}_VEP_spliceai.vcf",
        without_spliceAI=f"{config['dirs']['outdir']}/vep/{{caller}}/{{sample}}_{{aligner}}_VEP_nosplice.vcf"
    shell:
        r"""
        spliceAI_pattern='[-0-9]+\|[-0-9]+\|[-0-9]+\|[-0-9]+\|[0-9.]+\|[0-9.]+\|[0-9.]+\|[0-9.]+\|'

        rm -f {output.with_spliceAI} {output.without_spliceAI}

        with_spliceAI_count=0
        without_spliceAI_count=0

        while IFS= read -r line; do
            if [[ "$line" =~ ^# ]]; then
                echo "$line" >> {output.with_spliceAI}
                echo "$line" >> {output.without_spliceAI}
            elif [[ "$line" =~ $spliceAI_pattern ]]; then
                echo "$line" >> {output.with_spliceAI}
                with_spliceAI_count=$(expr $with_spliceAI_count + 1)
            else
                echo "$line" >> {output.without_spliceAI}
                without_spliceAI_count=$(expr $without_spliceAI_count + 1)
            fi
        done < {input}

        echo "VCF splitting complete!"
        echo "Lines with spliceAI annotation: $with_spliceAI_count (saved in: {output.with_spliceAI})"
        echo "Lines without spliceAI annotation: $without_spliceAI_count (saved in: {output.without_spliceAI})"
        """


rule spliceai_merge:
    input:
        with_spliceAI=f"{config['dirs']['outdir']}/vep/{{caller}}/{{sample}}_{{aligner}}_VEP_spliceai.vcf",
        without_spliceAI=f"{config['dirs']['outdir']}/vep/{{caller}}/{{sample}}_{{aligner}}_VEP_nosplice.vcf"
    output:
        inter=f"{config['dirs']['outdir']}/spliceai/{{caller}}/{{sample}}_{{aligner}}.vcf",
        final=f"{config['dirs']['outdir']}/spliceai/{{caller}}/{{sample}}_{{aligner}}_total.vcf"
    params:
        ref = config['refs']['genome'],
        container = config["singularity"]["spliceai"]
    threads: config['threads']
    resources:
        gpu_slot=1
    shell:
        r"""
        mkdir -p $(dirname "{output.final}")

        safe_vcf="{output.final}.safe.vcf"
        unsafe_body="{output.final}.unsafe.body.vcf"

        rm -f "$safe_vcf" "$unsafe_body" "{output.inter}" "{output.final}"

        # SAFE = biallelic AND (SNV or simple indel)
        awk -F'\t' 'BEGIN{{OFS="\t"}} \
          /^#/ {{print > safe; next}} \
          {{ref=$4; alt=$5; \
            if(index(alt,",")>0) {{print > unsafe; next}} \
            if(length(ref)==1 || length(alt)==1) {{print > safe}} \
            else {{print > unsafe}} }}' \
          safe="$safe_vcf" unsafe="$unsafe_body" "{input.without_spliceAI}"

        # Run SpliceAI on SAFE variants only
        if [ "$(grep -vc '^#' "$safe_vcf")" -gt 0 ]; then
          singularity exec --nv "{params.container}" spliceai \
            -I "$safe_vcf" \
            -O "{output.inter}" \
            -R "{params.ref}" \
            -A grch38
        else
          # No safe variants -> keep header-only file as inter
          cp "$safe_vcf" "{output.inter}"
        fi

        # Final output:
        # - header from spliceai output
        # - variants already annotated by VEP plugin split (with_spliceAI)
        # - variants scored by spliceai CLI (inter body)
        # - unsafe variants (kept unscored)
        grep "^#" "{output.inter}" > "{output.final}"
        grep -v "^#" "{input.with_spliceAI}" >> "{output.final}" || true
        grep -v "^#" "{output.inter}" >> "{output.final}" || true
        cat "$unsafe_body" >> "{output.final}" || true

        # Sort final
        awk '/^#/ {{print; next}} {{print | "sort -k1,1V -k2,2n"}}' "{output.final}" > "{output.final}.sorted"
        mv "{output.final}.sorted" "{output.final}"

        rm -f "$safe_vcf" "$unsafe_body"
        echo "Merge complete! Output saved in: {output.final}"
        """


rule add_vaf_vcf:
    input:
        vcf=f"{config['dirs']['outdir']}/spliceai/{{caller}}/{{sample}}_{{aligner}}_total.vcf",
        spip=f"{config['dirs']['outdir']}/spip/{{caller}}/{{sample}}_{{aligner}}.csv",
    output:
        f"{config['dirs']['outdir']}/{{sample}}_{{aligner}}_{{caller}}_final.vcf"
    params:
        container = config['singularity']['python'],
        greendb = config['databases']['greendb'],
        script = f"{workflow.basedir}/add_vaf_vcf.py",
    threads: config['threads']
    shell:
        """
        singularity exec --nv {params.container} bash -c 'export TMPDIR=/mnt/csf/tmp && python3 {params.script} {input.vcf} {params.greendb} {output} {input.spip}'
        """


rule vcf_to_excel:
    input:
        lambda wc: expand(
            f"{config['dirs']['outdir']}/{{sample}}_{{aligner}}_{wc.caller}_final.vcf",
            sample=config["samples"],
            aligner=ALIGNERS,
        )
    output:
        unfiltered=f"{config['dirs']['outdir']}/unfiltered_final_snv_{{caller}}.xlsx",
        filtered=f"{config['dirs']['outdir']}/filtered_final_snv_{{caller}}.xlsx",
    params:
        container = config['singularity']['python'],
        db=config['databases']['LR_SNV'],
        prom=config['refs']['bed_prom'],
        occu=config['databases']['LR_occu'],
        script = f"{workflow.basedir}/vcf_to_excel.py",
    threads: config['threads']
    shell:
        r"""
        # openpyxl reopens temporary files; use local scratch instead of the SSHFS share.
        singularity exec --nv {params.container:q} bash -c '
          set -euo pipefail
          TMPDIR=$(mktemp -d /tmp/lr-excel.XXXXXX)
          export TMPDIR
          trap "rm -rf -- \"$TMPDIR\"" EXIT
          python3 "$@"
        ' bash {params.script:q} {input:q} {output.filtered:q} {output.unfiltered:q} {params.db:q} {params.prom:q} {params.occu:q}
        """
