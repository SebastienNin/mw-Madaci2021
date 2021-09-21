localrules: config_targets, config_qc_targets, config_blc2fastq_targets
rule config_targets:
    input: mwconf['targets']

rule config_qc_targets:
    input: mwconf['qc_targets']

rule config_blc2fastq_targets:
    input: mwconf['bcl2fastq_targets']

#localrules: config_targets_analysis
#rule config_targets_analysis:
#    """
#    Entry point:
#        ./mw-thymus/src/bash/process_analysis_summary.sh
#    """
#    input:
#        mwconf['targets_analysis']

rule ln_sst_table:
    """
    Created:
        2017-12-26 15:17:50
    Modified:
        2018-02-13 09:18:43
    Aim:
        Produces all files correctly described in sequencing_summary and map them to sst.
    Test:
        sm out/wget/salva_runs.xlsx -f; sm out/python/xlsx2tsv/wget/salva_runs.tsv; smqb ln_sst_table
        rsync -aruzPv /gpfs/projects/spicuglia/mw/ Charbonnier@10.1.1.26:/mnt/SynologySpicuglia/guillaume/users/gcharbonnier/mw --include-from src/rsync/include-from/sst.txt
    """
    input:
        tsv=TSV,
        #tsv="src/snakemake/tables/salva_runs.tsv",
        dependencies=input_ln_sst_table
    output:
        tsv="out/ln/sst_table.tsv"
    params:
        tmp=TMP,
        tsv=TSV
    run:
        #shell("rm -rf sst")
        mapper_mw_to_sst(usage='run',tsv=TSV)
        #shell("chgrp -R thymus sst")
        shell("rm -rf {params.tmp}")

#rule ln_sst_table_test:
#    """
#    Created:
#        2019-02-27 10:47:23
#    Aim:
#        Testing rule. Does the same as ln_sst_table but only for a subset of samples.
#        Faster to debug than generating the DAG for all samples.
#    Test:
#    """
#    input:
#        tsv=TSV_TEST,
#        dependencies=input_ln_sst_table
#    output:
#        tsv="out/ln/sst_table.tsv"
#    params:
#        tmp=TMP,
#        tsv=TSV
#    run:
#        mapper_mw_to_sst(usage='run',tsv=TSV)
#        shell("rm -rf {params.tmp}")

#rule salva_2018_01_09:
#    """
#    Created:
#        2018-01-09 12:57:22
#    """
#    input:
#        "out/ChromHMM/LearnModel_H3K27ac-in-thymic-populations_numstates-20_assembly-hg38/emissions_20.png",
#        expand("out/deepTools/bamCoverage_binSize-20_normalizeUsingRPKM_extendReads-200/inp/bam/hg38/H3K27ac/thymus/{sample}.bw", sample=THYMUS_SAMPLES)#,
#        #expand("out/")

#rule salva_2018_07_30:
#    """
#    Created:
#        2018-07-30 13:26:48
#    """
#    input:
#        expand("out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/Quantile_normalized_VS_SP8_bw_hg38/{mark}_{stage}.bw", mark=SIX_MAIN_HISTONE_MARKS, stage=ALT_MHSC_THYMIC_STAGES),
#        expand("out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/{sample}.bed", sample=THYMIC_STAGES),
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI.CD4+_CD25+_CD127-_Treg_62.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD4+_CD25-_CD45RA+_Naive_62.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD4+_CD25-_CD45RO+_Memory_62.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD4+_CD25-_IL17+_PMA-Ionomcyin_stimulated_Th17_62.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD4+_CD25int_CD127+_Tmem_332.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD4+_CD25int_CD127+_Tmem_62.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD4+_CD25-_Th_332.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD4_Memory_100_7.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD4_Memory_100_8.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD4_Naive_100_7.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD4_Naive_101_8.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD8_Memory_100_7.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD8_Naive_100_7.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_CD8_Naive_100_8.bed",
#        "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_Mobilized_CD34_1480.bed",
#       "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_Mobilized_CD34_1536.bed",
#       "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/BI_Mobilized_CD34_1549.bed",
#       "out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/ChromHMM_hg38/ROADMAP/UW_CD14.bed"
#   shell:
#        """
#        cp -rv out/ln/alias/to_copy_to_Blueprint_data_processed_for_Salva/* /gpfs/projects/spicuglia/Blueprint_processed_data/
#        """

"""
tmp <- c(
"MONO_CB_1_subsample_size-3439_replicate-1",
"MONO_CB_1_subsample_size-3439_replicate-2",
"MONO_CB_1_subsample_size-3439_replicate-3",
"MONO_CB_2_subsample_size-3439_replicate-1",
"MONO_CB_2_subsample_size-3439_replicate-2",
"MONO_CB_2_subsample_size-3439_replicate-3",
"MONO_PB_subsample_size-3439_replicate-1",
"MONO_PB_subsample_size-3439_replicate-2",
"MONO_PB_subsample_size-3439_replicate-3",
"MONO_ENCODE_1_subsample_size-3439_replicate-1",
"MONO_ENCODE_1_subsample_size-3439_replicate-2",
"MONO_ENCODE_1_subsample_size-3439_replicate-3",
"MONO_ENCODE_2_subsample_size-3439_replicate-1",
"MONO_ENCODE_2_subsample_size-3439_replicate-2",
"MONO_ENCODE_2_subsample_size-3439_replicate-3",
"MONO_ROADMAP_subsample_size-3439_replicate-1",
"MONO_ROADMAP_subsample_size-3439_replicate-2",
"MONO_ROADMAP_subsample_size-3439_replicate-3",
"IM_1_subsample_size-3439_replicate-1",
"IM_1_subsample_size-3439_replicate-2",
"IM_1_subsample_size-3439_replicate-3",
"IM_2_subsample_size-3439_replicate-1",
"IM_2_subsample_size-3439_replicate-2",
"IM_2_subsample_size-3439_replicate-3",
"IM_3_subsample_size-3439_replicate-1",
"IM_3_subsample_size-3439_replicate-2",
"IM_3_subsample_size-3439_replicate-3",
"AAM_subsample_size-3439_replicate-1",
"AAM_subsample_size-3439_replicate-2",
"AAM_subsample_size-3439_replicate-3")
"""


rule salva_2018_10_31:
    """
    """
    input:
        histone_marks=expand("out/deepTools/plotHeatmap_--sortRegions_keep{perGroup}_--boxAroundHeatmaps_no_-colorList_blueCyanYellowOrangeRed_-whatToShow_phc_-xAxisLabel_peak-center_-refPointLabel_0/deepTools/computeMatrix_reference-point_--referencePoint_center_-b_5000_-a_5000_-bs_200_--sortRegions_keep_-R_{regions}_-S_hg38-{mark}-thymus-merged-wiq.pdf", perGroup=["","_--perGroup"], mark=["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3","H3K9me3","H3K36me3"], regions=["hg38-dClust-rowFeature-no-rmsk-mxy-no-donor-effect-distal","hg38-dClust-rowFeature-all-cpg-hypo-meth-call-no-rmsk-mxy-no-donor-effect-distal"]),
        bs=expand("out/deepTools{v2}/plotHeatmap_--sortRegions_no{perGroup}_--missingDataColor_1_-colorList_blueCyanYellowOrangeRed_-whatToShow_phc_-xAxisLabel_peak-center_-refPointLabel_0/deepTools{v2}/computeMatrix_reference-point_--referencePoint_center_-b_{ba}_-a_{ba}{bs}_--sortRegions_keep_-R_{regions}_-S_hg38-BS-thymus-110-91.pdf",
            perGroup=["","_--perGroup"],
            bs=[""],
            ba=["5000"],
            v2=["2"],
            regions=["hg38-dClust-rowFeature-no-rmsk-mxy-no-donor-effect-distal","hg38-dClust-rowFeature-all-cpg-hypo-meth-call-no-rmsk-mxy-no-donor-effect-distal","hg38-bs-hypometh-thymus-common-all"])
        #expand("out/deepTools/plotHeatmap_sortRegions-descend_sortUsing-region_length_averageTypeSummaryPlot-mean_missingDataColor-1_colorList-blueCyanYellowOrangeRed_heatmapHeight-28_heatmapWidth-1_whatToShow-phc_xAxisLabel-peak-center_refPointLabel-0_boxAroundHeatmaps-no/deepTools/computeMatrix_reference-point_--referencePoint_center_-b_{ba}_-a_{ba}_-bs_200_--sortRegions_keep_-R_{regions}_-S_hg38-BS-thymus-110-91.pdf",
        #    regions=["hg38-dClust-rowFeature-no-rmsk-mxy-no-donor-effect-distal","hg38-dClust-rowFeature-all-cpg-hypo-meth-call-no-rmsk-mxy-no-donor-effect-distal","hg38-bs-hypometh-thymus-common-all"],
        #    ba=["2000","5000"])
    output:
        "out/ln/salva_2018_10_31/done"
    shell:
        """

        for FILE in {input};
        do
            BASENAME=$(basename $FILE)
            MD5=$(echo $(dirname $FILE) | md5sum | cut -f1 -d ' ')
            ln -f $FILE out/ln/salva_2018_10_31/${{MD5}}_${{BASENAME}}
            ln -f $FILE ../guillaume/plots/003_histones_bs_heatmaps_by_windows/${{MD5}}_${{BASENAME}}
        done
        touch {output}
        """


rule salva_2017_04_26:
    """
    Created:
        2017-05-04 13:54:13
    Mail:
        Voici toutes les infos pour l’analyse des données de l’atelier.
        Nous aurons 2 expériences de capstarr-seq faites avec la librairie Human Promoters
        Fichier d’annotation: CRM_hPromoters (ci-joint): contient 4 classes (hPro; Random; Positive; PosEpromoters)
        Genome: HG19
        Extension: 314
        Input: from run107
        Filtre: pas de filtre
        Critère de sélection: inflection point

        Fichiers de sortie:
        Fichier PDF ci-joint (plutôt png)
        Fichier de résultats « allData_all.tab » si possible sans les colonnes « groupes FDR » mis avec la colonne « Gene name »
        Fichier .bw normalisé
        Tu peux faire un test avec l'une de données K562
    Note:
        fichier de crms ici: inp/bed/crms/hProm_posEprom.bed
    """
    input:
        "out/capstarrseq/merge_all_data_hProm_posEprom/awk/extend_reads_314/bedtools/bamtobed/samtools/sam_to_bam/bowtie2/se_hg19/sickle/se_-t_sanger_-q_20/gunzip/merge_lanes_nextseq500_single_end/ln/rename_run107_tgml/K562_rep1_over_Input.allData.tsv",
        "out/poppler/pdftoppm_png_singlefile/capstarrseq/grouping_crms/capstarrseq/fold_change_hProm_posEprom/awk/extend_reads_314/bedtools/bamtobed/samtools/sam_to_bam/bowtie2/se_hg19/sickle/se_-t_sanger_-q_20/gunzip/merge_lanes_nextseq500_single_end/ln/rename_run107_tgml/K562_rep1_over_Input.inflexionPointGroups.png",
        expand("out/deepTools/bamCoverage_binSize-20_minMappingQuality-0_normalizeUsingRPKM/samtools/sam_to_bam/bowtie2/se_hg19/sickle/se_-t_sanger_-q_20/gunzip/merge_lanes_nextseq500_single_end/ln/rename_run107_tgml/{samples}.bw", samples=["K562_rep1","Input"]),
        "out/deepTools/bamCompare_scaleFactorsMethod-readCount_ratio-log2_binSize-20_smoothLength-200_minMappingQuality-0/samtools/sam_to_bam/bowtie2/se_hg19/sickle/se_-t_sanger_-q_20/gunzip/merge_lanes_nextseq500_single_end/ln/rename_run107_tgml/K562_rep1_vs_Input.bw"
    params:
        shr_dir="shr/salva_2017_04_26"
    shell:
        """
        mkdir -p {params.shr_dir}
        for FILE in {input}
        do
            BASENAME=`basename $FILE`
            ln --force $FILE {params.shr_dir}/$BASENAME
        done
        """

rule salva_2017_04_10:
    """
    Created:
        2017-04-10 14:35:18
    Mail:
        All Samples for Run178

        mT-DHS-PGK-p5424_rep3
        Extension: 314
        Mouse mm9
        input: input_plasmid_library (mT_DHS_PGK_input) from run 150
        CRM : mouse DHS (comme pour le run 150)
        No filter on input FPKM

        Il faut aligner chez l’homme (H38) ou souris (mm9) et générer les bigwig normalisés par le nombre de reads sans enlever les « artefacts PCR » avec ou pas ratio input
        Human
        h_SE_SCP1_input
        h_SE_SCP1_Jurkat_rep1
        h_SE_pRag2_input
        h_SE_pRag2_Jurkat_rep1
        h_SE_pRag2_Jurkat_rep2

        Souris
        BAC-mT-DHS-PGK-p5424_rep1
        BAC-mT-DHS-PGK-p5424_rep2
        BAC-mT-DHS-PGK-input

    Note:
        ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz

    Test:
        out/19_mTDHS/awk/extend_reads_314/bedtools/bamtobed/samtools/sam_to_bam/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_20/gunzip/ln/rename_run178_tgml/BAC_mT_DHS_PGK_p5424_rep1_over_BAC_mT_DHS_PGK_input.allData.tsv
   """
    input:
        # BW
        ## Raw
        expand("out/deepTools/bamCoverage_binSize-20_minMappingQuality-0_normalizeUsingRPKM/samtools/sam_to_bam/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_20/gunzip/ln/rename_run178_tgml/{mouse_samples}.bw", mouse_samples=["BAC_mT_DHS_PGK_p5424_rep1","BAC_mT_DHS_PGK_p5424_rep2","BAC_mT_DHS_PGK_input"]),
        expand("out/deepTools/bamCoverage_binSize-20_minMappingQuality-0_normalizeUsingRPKM/samtools/sam_to_bam/bowtie2/pe_GRCh38/sickle/pe_-t_sanger_-q_20/gunzip/ln/rename_run178_tgml/{human_samples}.bw", human_samples=["h_SE_SCP1_input","h_SE_SCP1_Jurkat_rep1","h_SE_pRag2_input","h_SE_pRag2_Jurkat_rep1","h_SE_pRag2_Jurkat_rep2"]),
        # ADDD EXTEND READS RAW AND RATIO FOR mT_DHS_PGK.
        ## Ratios
        ### Standard smoothed
        expand("out/deepTools/bamCompare_scaleFactorsMethod-readCount_ratio-log2_binSize-20_smoothLength-200_minMappingQuality-0/samtools/sam_to_bam/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_20/gunzip/ln/rename_run178_tgml/{id_bam1}_vs_{id_bam2}.bw", id_bam1=["BAC_mT_DHS_PGK_p5424_rep1","BAC_mT_DHS_PGK_p5424_rep2"], id_bam2=["BAC_mT_DHS_PGK_input"]),
        expand("out/deepTools/bamCompare_scaleFactorsMethod-readCount_ratio-log2_binSize-20_smoothLength-200_minMappingQuality-0/samtools/sam_to_bam/bowtie2/pe_GRCh38/sickle/pe_-t_sanger_-q_20/gunzip/ln/rename_run178_tgml/{id_bam1}_vs_{id_bam2}.bw", id_bam1=["h_SE_SCP1_Jurkat_rep1"], id_bam2=["h_SE_SCP1_input"]),
        expand("out/deepTools/bamCompare_scaleFactorsMethod-readCount_ratio-log2_binSize-20_smoothLength-200_minMappingQuality-0/samtools/sam_to_bam/bowtie2/pe_GRCh38/sickle/pe_-t_sanger_-q_20/gunzip/ln/rename_run178_tgml/{id_bam1}_vs_{id_bam2}.bw", id_bam1=["h_SE_pRag2_Jurkat_rep1","h_SE_pRag2_Jurkat_rep2"], id_bam2=["h_SE_pRag2_input"]),
        ### blacklist filtered smoothed
        expand("out/deepTools/bamCompare_scaleFactorsMethod-readCount_ratio-log2_binSize-20_blackListFileName-no-super-enhancers-hg38_smoothLength-200_minMappingQuality-0/samtools/sam_to_bam/bowtie2/pe_GRCh38/sickle/pe_-t_sanger_-q_20/gunzip/ln/rename_run178_tgml/{id_bam1}_vs_{id_bam2}.bw", id_bam1=["h_SE_SCP1_Jurkat_rep1"], id_bam2=["h_SE_SCP1_input"]),
        expand("out/deepTools/bamCompare_scaleFactorsMethod-readCount_ratio-log2_binSize-20_blackListFileName-no-super-enhancers-hg38_smoothLength-200_minMappingQuality-0/samtools/sam_to_bam/bowtie2/pe_GRCh38/sickle/pe_-t_sanger_-q_20/gunzip/ln/rename_run178_tgml/{id_bam1}_vs_{id_bam2}.bw", id_bam1=["h_SE_pRag2_Jurkat_rep1","h_SE_pRag2_Jurkat_rep2"], id_bam2=["h_SE_pRag2_input"]),
        # Capstarrseq analysis
        #"out/19_mTDHS/awk/extend_reads_314/bedtools/bamtobed/samtools/sam_to_bam/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_20/gunzip/ln/rename_run178_tgml/BAC_mT_DHS_PGK_p5424_rep1_over_BAC_mT_DHS_PGK_input.allData.tsv",
        # Even if the run 178 is paired-end, this sample has to be considered se to be compared to the se input from run 155.
        "out/19_mTDHS/awk/extend_reads_314/bedtools/bamtobed/samtools/sam_to_bam/bowtie2/se_mm9/sickle/se_-t_sanger_-q_20/gunzip/ln/rename_run178_tgml/mT_DHS_PGK_p5424_rep3_1_over_mT_DHS_PGK_input.allData.tsv",
        # Here is just the reanalysis of the rep1 from run155 to compare with the results from fq_to_bam and see if the two workflows produce the same results.
        "out/19_mTDHS/awk/extend_reads_314/bedtools/bamtobed/samtools/sam_to_bam/bowtie2/se_mm9/sickle/se_-t_sanger_-q_20/gunzip/input/fastq/run155/CapSTAR-seq_mSilencers_mT_DHS_PGK_P5424_rep1_over_CapSTAR-seq_mSilencers_mT_DHS_PGK_input.allData.tsv",
    params:
        shr_dir="shr/salva_2017_04_10"
    shell:
        """
        mkdir -p {params.shr_dir}
        for FILE in {input}
        do
            BASENAME=`basename $FILE`
            ln --force $FILE {params.shr_dir}/$BASENAME
        done

        """

CONDITIONS_RUN168_RNASEQ=[
    "T4",
    "T5",
    "T6",
    "T8",
    "T9",
    "T10",
    "T11",
    "T12",
    "CD34"]

CONDITIONS_RUN151=[
    "POMMIER",
    "SCHWEIZER",
    "CHEVALIER",
    "BOUDJELTHIA",
    "MARINHO",
    "EVRARD",
    "LEFRANCOIS",
    "CROZAT",
    "HALLOU",
    "SCHOLZ",
    "MONTIER",
    "BAGNOST"]

FILTER_RNA_STRAND=["filterRNAstrand-forward_","filterRNAstrand-reverse_",""]

rule salva_2017_12_05:
    """
    Created:
        2017-12-05 16:51:40
    Aim:
        First step to uniformized directory for Salva's projects.
    """
    input:
        expand(
            "out/deepTools/bamCoverage_{filterRNAstrand}binSize-20_minMappingQuality-0_normalizeUsingRPKM/star/se_GRCh38/sickle/se_-t_sanger_-q_20/gunzip/merge_lanes_nextseq500_single_end/ln/alias/fastq/run168/RNA-seq/{condition}.bw",
            condition=CONDITIONS_RUN168_RNASEQ,
            filterRNAstrand=FILTER_RNA_STRAND),
        expand(
            "out/deepTools/bamCoverage_{filterRNAstrand}binSize-20_minMappingQuality-0_normalizeUsingRPKM/star/pe_GRCh38/sickle/pe_-t_sanger_-q_20/gunzip/to-stdout/ln/alias/fastq/run151/{condition}.bw",
            condition=CONDITIONS_RUN151,
            filterRNAstrand=FILTER_RNA_STRAND)
