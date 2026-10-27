import json

rule nanostat:
    input: 
        f"{config['dirs']['outdir']}/alignment/minimap2/{{sample}}.bam"
    output:
        f"{config['dirs']['outdir']}/multiqc/{{sample}}_nanostat.txt"
    params:
        outdir = directory(f"{config['dirs']['outdir']}/multiqc"),
        sname = f"{{sample}}_nanostat.txt"
    threads: config['threads']
    singularity: 
        config["singularity"]["nanostat"]
    shell:
        """
        NanoStat \
            --bam {input} \
            --outdir {params.outdir} \
            -n {params.sname} \
            -t {threads}
        """


rule samtools_stat:
    input: 
        f"{config['dirs']['outdir']}/alignment/minimap2/{{sample}}.bam"
    output:
        stat_ontarget = f"{config['dirs']['outdir']}/multiqc/samtools_stat/{{sample}}_ontarget.stat",
        stat_offtarget = f"{config['dirs']['outdir']}/multiqc/samtools_stat/{{sample}}_offtarget.stat"
    params:
        bed_ontarget = config['refs']['bed_ontarget'],
        bed_offtarget = config['refs']['bed_offtarget']
    threads: config['threads']
    singularity: 
        config["singularity"]["samtools"]
    shell:
        """
        samtools stats \
            -@ {threads} \
            {input} \
            -t {params.bed_ontarget} \
                > {output.stat_ontarget}
        
        samtools stats \
            -@ {threads} \
            {input} \
            -t {params.bed_offtarget} \
                > {output.stat_offtarget}
        """


def get_dico_json(section_name, description, xlab, ylab):
    dico_json = {
    "parent_id": "Adaptive sampling statistics",
    "parent_name": "Adaptive sampling statistics",
    "parent_description": "Various statistics of adaptive sampling results. Using `samtools stats` on bedfile used during the adaptive sampling (i.e. on target region) and on the bedtool complement of this bedfile (off target regions).",
    "plot_type": "linegraph",
    "section_name": section_name,
    "description": description,
    "pconfig": {
        "title": "",
        "xlab": xlab,
        "ylab": ylab,
        "categories": bool("true"),
        "colors": {}
        }
    }
    return dico_json


rule generate_cov_rl_files:
    input:
        ontarget_files = expand(f"{config['dirs']['outdir']}/multiqc/samtools_stat/{{sample}}_ontarget.stat", sample=config['samples']),
        offtarget_files = expand(f"{config['dirs']['outdir']}/multiqc/samtools_stat/{{sample}}_offtarget.stat", sample=config['samples'])
    output:
        coverages = f"{config['dirs']['outdir']}/multiqc/coverage_FlowcellA_on_off_ratio_target_mqc.json",
        read_lengths = f"{config['dirs']['outdir']}/multiqc/read_lengths_FlowcellA_on_off_ratio_target_mqc.json"
    params:
        outdir = f"{config['dirs']['outdir']}/multiqc/"
    run:
        all_onfiles = list(sorted(input.ontarget_files))
        flowcell_files = [all_onfiles[:4], all_onfiles[4:]]

        flowcount = 0
        for ontarget_files in flowcell_files:
            flowcount +=1
            tag = "FlowcellA" if flowcount ==1 else "FlowcellB"

            dico_json_RL = get_dico_json(f"Read lengths On/Off target - {tag}", "Linegraph of read lengths of on and off target regions per sample.", "Range of read lengths (bp)", "Nb of reads (ratio %)"); dico_json_RL["data"] = {}
            dico_json_COV = get_dico_json(f"Coverage On/Off target - {tag}", "Linegraph of coverage of on and off target regions per sample.", "Coverage", "Nb of reads (ratio %)"); dico_json_COV["data"] = {}
            dico_json_COV["pconfig"]["xmax"] = 100
            lst_colors = ["#ff1744", "#aeea00", "#1a237e", "#00bcd4", "#00c853", "#ffab00", "#e64a19", "#6d4c41"]

            if ontarget_files:
                for file in ontarget_files:
                    ix = ontarget_files.index(file)
                    sample = file.split("_ontarget.stat")[0]
                    sample_files = [sample+"_ontarget.stat", sample+"_offtarget.stat"]

                    #compute read lengths and coverage stat
                    lst_range = [i for i in range(0, 10500, 500)]
                    dico_RL = dico_RL_ratio = {f"{sample}_On-Target": {r: 0 for r in lst_range}, f"{sample}_Off-Target": {r: 0 for r in lst_range}}
                    dico_cov = dico_cov_ratio = {f"{sample}_On-Target": {}, f"{sample}_Off-Target": {"0": 0}}
                    c = 0; ON_RL_sum = 0; OFF_RL_sum = 0; ON_COV_sum = 0; OFF_COV_sum = 0
                    for stat_file in sample_files:
                        c += 1
                        reg = "On-Target" if c == 1 else "Off-Target"
                        dico_json_RL["pconfig"]["colors"][f"{sample}_{reg}"] = dico_json_COV["pconfig"]["colors"][f"{sample}_{reg}"] = lst_colors[ix]
                        with open(stat_file, "r") as filin:
                            for line in filin:
                                line = line.split()
                                if line[0] == "RL":
                                    for i in lst_range:
                                        ind = lst_range.index(i)
                                        try:
                                            if int(line[1]) in range(lst_range[ind], lst_range[ind+1]):
                                                dico_RL[f"{sample}_{reg}"][lst_range[ind]] += int(line[2])
                                                ON_RL_sum += int(line[2]) if c == 1 else 0
                                                OFF_RL_sum += int(line[2]) if c == 2 else 0
                                                break
                                        except IndexError:
                                            if int(line[1]) > lst_range[-1]:
                                                dico_RL[f"{sample}_{reg}"][lst_range[-1]] += int(line[2])
                                                ON_RL_sum += int(line[2]) if c == 1 else 0
                                                OFF_RL_sum += int(line[2]) if c == 2 else 0
                                elif line[0] == "COV":
                                    dico_cov[f"{sample}_{reg}"][line[2]] = int(line[3])
                                    ON_COV_sum += int(line[3]) if c == 1 else 0
                                    OFF_COV_sum += int(line[3]) if c == 2 else 0
                    
                    for region, dic_reg in dico_RL.items():
                        for r, count in dic_reg.items():
                            if region.split('/')[-1].split('_')[1] == "On-Target":
                                dico_RL_ratio[region][r] = (count/ON_RL_sum)*100
                            elif region.split('/')[-1].split('_')[1] == "Off-Target":
                                dico_RL_ratio[region][r] = (count/OFF_RL_sum)*100

                    dico_json_RL["data"].update(dico_RL_ratio)
                    
                    for region, dic_reg in dico_cov.items():
                        for r, count in dic_reg.items():
                            if region.split('/')[-1].split('_')[1] == "On-Target":
                                dico_cov_ratio[region][r] = (count/ON_COV_sum)*100
                            elif region.split('/')[-1].split('_')[1] == "Off-Target":
                                dico_cov_ratio[region][r] = (count/OFF_COV_sum)*100
                    
                    dico_json_COV["data"].update(dico_cov_ratio)

                with open(f"{params.outdir}coverage_{tag}_on_off_ratio_target_mqc.json", "w") as filout:
                    json.dump(dico_json_COV, filout, indent=6)
                
                with open(f"{params.outdir}read_lengths_{tag}_on_off_ratio_target_mqc.json", "w") as filout:
                    json.dump(dico_json_RL, filout, indent=6)

rule samtools_depth:
    input:
        bam = f"{config['dirs']['outdir']}/alignment/minimap2/{{sample}}.bam",
        bai = f"{config['dirs']['outdir']}/alignment/minimap2/{{sample}}.bam.bai",
        bed = config['refs']['bed_exons'],
    output:
        depth = f"{config['dirs']['outdir']}/multiqc/{{sample}}.depth",
    singularity:
        config["singularity"]["samtools"]
    shell:
        """
        samtools depth -b {input.bed} {input.bam} > {output.depth}
        """

rule coverage_per_gene:
    input:
        bed_exons = config['refs']['bed_exons'],
        bed_locus = config['refs']['bed_locus'],
        depths = expand(f"{config['dirs']['outdir']}/multiqc/{{sample}}.depth", sample=config['samples']),
    output:
        all30_ex = f"{config['dirs']['outdir']}/multiqc/exons_30X_mqc.tsv",
        all30_loc = f"{config['dirs']['outdir']}/multiqc/locus_30X_mqc.tsv",
    params:
        depths_dir = f"{config['dirs']['outdir']}/multiqc",
        script = f"{workflow.basedir}/coverage.py",
    singularity:
        config['singularity']['python']
    shell:
        """
        python3 {params.script} {input.bed_exons} {params.depths_dir} exons
        python3 {params.script} {input.bed_locus} {params.depths_dir} locus
        rm {params.depths_dir}/*.depth
        """

rule generate_multiQC_report:
    input:
        f"{config['dirs']['outdir']}/multiqc/read_lengths_FlowcellA_on_off_ratio_target_mqc.json",
        expand(f"{config['dirs']['outdir']}/multiqc/{{sample}}_nanostat.txt", sample=config['samples']),
        f"{config['dirs']['outdir']}/multiqc/exons_30X_mqc.tsv",
    output:
        f"{config['dirs']['outdir']}/multiqc/DNASeq_LongRead_Adaptive_Sampling_multiqc_report.html"
    params:
        dirr=directory(f"{config['dirs']['outdir']}/multiqc/"),
        configpath= f"{workflow.basedir}/../config/mqc_config.yaml",
    threads: config['threads']
    singularity: 
        config["singularity"]["multiqc"]
    shell:
        """
        export PYTHONNOUSERSITE=1;\
        multiqc \
            --force \
            --export \
            --interactive \
            --verbose \
            --outdir {params.dirr} \
            --title DNASeq_LongRead_Adaptive_Sampling \
            --config {params.configpath} \
            {params.dirr}
        """
