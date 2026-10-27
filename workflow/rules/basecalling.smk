import json


rule dorado_basecall:
    input:
        raw_data = config['dirs']['indir']
    output:
        basecalled_reads = temp(f"{config['dirs']['outdir']}/dorado_basecalled_reads.fastq")
    params:
        outdir = config['dirs']['outdir'],
        container = config['singularity']['dorado'],
        kit = config['tools']['dorado']['kit'],
        model = config['tools']['dorado']['model'],
    threads: config['threads']
    resources:
        gpu_slot=1
    shell:
        """
        mkdir -p {params.outdir:q}
        singularity exec --nv {params.container:q} \
            dorado basecaller {params.model:q} {input.raw_data:q} \
            -x cuda:all --recursive --emit-fastq --no-trim --kit-name {params.kit:q} \
            > {output.basecalled_reads:q}
        """


rule dorado_demux:
    input:
        reads = rules.dorado_basecall.output.basecalled_reads
    output:
        demux_dir = temp(directory(f"{config['dirs']['outdir']}/dorado_demux_raw"))
    params:
        container = config['singularity']['dorado'],
        kit = config['tools']['dorado']['kit'],
    threads: config['threads']
    shell:
        """
        mkdir -p {output.demux_dir:q}
        singularity exec {params.container:q} \
            dorado demux {input.reads:q} --output-dir {output.demux_dir:q} \
            --emit-fastq --kit-name {params.kit:q} -t {threads}
        """


rule dorado_postprocess:
    input:
        demux_dir = rules.dorado_demux.output.demux_dir,
        script = f"{workflow.basedir}/dorado_postprocess.py",
    output:
        expand(f"{config['dirs']['outdir']}/dorado_demux/{{sample}}.fastq.gz", sample=config['samples'])
    params:
        outdir = f"{config['dirs']['outdir']}/dorado_demux",
        mapping = json.dumps(config.get("barcode_to_sample", {})),
        samples = json.dumps(config['samples']),
    threads: config['threads']
    shell:
        """
        python3 {input.script:q} {input.demux_dir:q} {params.outdir:q} {params.mapping:q} {params.samples:q} --threads {threads}
        """
