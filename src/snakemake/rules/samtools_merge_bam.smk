rule samtools_merge_samples_sst:
    """
    Created:
        2018-01-24 18:03:14
    Aim:
        This rule use Samtools merge samples as define in salva_runs.tsv by sample_merge_list value.
    Test:
        out/samtools/merge_samples_sst/CapStarr_155_170_178_204/GRCm38/mT_DHS_PGK_input.bam
    """
    input:
        #tsv=TSV, Annoying to get analysis redone each time TSV is updated.
        bam=input_bam_samtools_merge_samples_sst
    output:
        bam="out/samtools/merge_samples_sst/{exp_name}/{assembly}/{sample_name}.bam"
    wildcard_constraints:
        exp_name="[a-zA-Z0-9_-]+",
        sample_name="[a-zA-Z0-9_-]+",
    threads:
        1
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools merge {output.bam} {input.bam}"

rule ln_sst_exp:
    """
    Created:
        2018-01-29 14:49:25
    Aim:
        Put samples from different runs but for same exp_name into a common directory for easier downstream analysis.
    Test:
        out/ln/sst_exp/CapStarr_155_170_178_204/GRCm38/mT_DHS_PGK_P5424_stimulated_rep3.bam
        out/ln/sst_exp/CapStarr_155_170_178_204/GRCm38/mT_DHS_PGK_input.bam
        out/ln/sst_exp/ChIP/mm9/P5424_WT_K27ac.bam
        out/ln/sst_exp/CapStarr_107/hg19/K562_Capstarr_hpromoter_rep1.bam
    """
    input:
        #tsv=TSV,
        bam=input_bam_ln_sst_exp
    output:
        bam="out/ln/sst_exp/{exp_name}/{assembly}/{sample_name}.bam"
    shell:
        "ln -f {input.bam} {output.bam}"

#def input_bam_ln_sst_pool(wildcards):
#    """
#    Created:
#        2018-02-11 11:54:46
#    Aim:
#        Pool all samples from Salva projects into the same directory to allow the use of SAMPLE VS CONTROL rules for different runs without having to define an 'exp_name' in the spreadsheet.
#    """
#    assembly = wildcards['assembly']
#    sample_name = wildcards['sample_name']
#
#    for row in csv.DictReader(open(TSV),delimiter='\t'):
#        if row['sample_name'] == wildcards['sample_name']:
#            if row['exp_name'] != wildcards['exp_name']:
#                raise Exception("This sample can't be in this experiment based on metadata.")
#            if ',' in row['sample_merge_list']:
#                path = "out/samtools/merge_samples_sst/" + exp_name + "/" + assembly + "/" + sample_name + ".bam"
#            else:
#                path = "out/samtools/sort/samtools/view_bSh/bowtie2/" + row['se_or_pe'] + "_" + assembly +"/sickle/" + row['se_or_pe'] + "_-t_sanger_-q_20/" + "gunzip/merge_lanes_nextseq500_" + row['se_or_pe'] + "_raw/" + row['tgml_fastq_prefix'] + ".bam"
#
#    return(path)
#
#
#
#rule ln_sst_pool:
#    """
#    Created:
#        2018-02-11 11:56:28
#    Aim:
#        Pool all samples from Salva projects into the same directory to allow the use of SAMPLE VS CONTROL rules for different runs without having to define an 'exp_name' in the spreadsheet.
#    Test:
#    """
#    input:
#        tsv=TSV,
#        bam=input_bam_ln_sst_pool
#    output:
#        bam="out/ln/sst_pool/{assembly}/{sample_name}.bam"
#    shell:
#        """
#        ln -f {input.bam} {output.bam}
#        """
