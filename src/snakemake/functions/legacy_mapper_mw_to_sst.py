def mapper_mw_to_sst(usage,tsv):
    """
    Created:
        2018-01-01 22:29:32
    Aim:
        Define links between mw files and sst ones.
        usage is either 'input' when the function is used as rule input to define dependecies, or 'run' when function is used inside the rule to create the hardlinks.
    """
    paths = []
    d_paths_in_to_ln = {}

    # I should choose between pandas and csv modules to handle the table... I try both for the moment.
    df = pandas.read_table(tsv, header=0, comment='#', skiprows=[1])
    #print(df.columns.values.tolist())

    # Checking sample_name are unique
    seen_sample_name = set()
    uniq_sample_name = []

    for sample_name in df.sample_name:
        if sample_name in seen_sample_name:
            print('duplicated sample_name: ' + str(sample_name))
        else:
            uniq_sample_name.append(sample_name)
            seen_sample_name.add(sample_name)


    d_metadata = csv.DictReader(open(tsv),delimiter='\t')
    data = []
    for row in d_metadata:
        data.append(row)

    #data = {row["sample_name"]: row["tgml_fastq_prefix"] for row in d_metadata}
    #print(data)

    for row in csv.DictReader(open(tsv),delimiter='\t'):
    #for row in d_metadata:

        # Only look for samples that have available files
        if re.match('^inp/fastq/run*', row['tgml_fastq_prefix']):
            #print("sample_name : " + row['sample_name'])

            if row['type'] == '':
                TYPE = 'Unknown'
            else:
                TYPE = row['type']

            if row['merged_or_unmerged'] == 'unmerged':
                ln_base_prefix = "sst/" + TYPE + "/run" + row['run']
                #in_base_suffix = "gunzip/merge_lanes_nextseq500_" + row['se_or_pe'] + "_raw/" + row['tgml_fastq_prefix']
                in_base_suffix = "gunzip/merge_lanes_nextseq500_" + row['se_or_pe'] + "_raw/ln/updir/mw/" + row['tgml_fastq_prefix']
                in_fq_prefix = "out/gzip/to-stdout/" + in_base_suffix

                #print('unmerged')
                if row['se_or_pe'] == 'se':
                    path_fq_in = in_fq_prefix + ".fastq.gz"
                    path_fq_ln = ln_base_prefix + "/fastq/" + row['sample_name'] + ".fastq.gz"
                    d_paths_in_to_ln[path_fq_in] = path_fq_ln
                if row['se_or_pe'] == 'pe':
                    path_fq_R1_in = in_fq_prefix + "_1.fastq.gz"
                    path_fq_R1_ln = ln_base_prefix + "/fastq/" + row['sample_name'] + "_R1.fastq.gz"
                    d_paths_in_to_ln[path_fq_R1_in] = path_fq_R1_ln

                    path_fq_R2_in = in_fq_prefix + "_2.fastq.gz"
                    #print(path_fq_R2_in)
                    path_fq_R2_ln = ln_base_prefix + "/fastq/" + row['sample_name'] + "_R2.fastq.gz"
                    d_paths_in_to_ln[path_fq_R2_in] = path_fq_R2_ln


                # path_bam_in="out/gunzip/merge_lanes_nextseq500_se_raw/" + row['tgml_path_prefix'] + ".fastq"
                if row['specie'] == 'human':
                    # Add other assemblies here:
                    assembly_list = ["GRCh38","hg19"]
                    gsize = "hs"
                elif row['specie'] == 'mouse':
                    assembly_list = ["GRCm38","mm9"]
                    gsize = "mm"

                # I do not want to throw an error if specie is not referred and just skip sample instead.
                if 'assembly_list' in locals():

                    for assembly in assembly_list:
                        ln_processed_prefix = ln_base_prefix + "/" + assembly

                        if row['type'] == 'RNA':
                            in_aligned_suffix = "star/" + row['se_or_pe'] + "_" + assembly +"/sickle/" + row['se_or_pe'] + "_-t_sanger_-q_20/" + in_base_suffix
                        else:
                            in_aligned_suffix = "samtools/sort/samtools/view_bSh/bowtie2/" + row['se_or_pe'] + "_" + assembly +"/sickle/" + row['se_or_pe'] + "_-t_sanger_-q_20/" + in_base_suffix

                        path_bam_in = "out/" + in_aligned_suffix + ".bam"
                        path_bam_ln = ln_processed_prefix + "/bam/" + row['sample_name'] + ".bam"
                        d_paths_in_to_ln[path_bam_in] = path_bam_ln

                        path_bai_in = path_bam_in + ".bai"
                        path_bai_ln = path_bam_ln + ".bai"
                        d_paths_in_to_ln[path_bai_in] = path_bai_ln

                        path_bw_in="out/deepTools/bamCoverage_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM/" + in_aligned_suffix + ".bw"
                        path_bw_ln = ln_processed_prefix + "/bw/" + row['sample_name'] + ".bw"
                        d_paths_in_to_ln[path_bw_in] = path_bw_ln

                        if row['type'] in ['ChIP','ATAC']:
                            # Peak-calling without control:
                            #"out/macs2/callpeak_broad_format-BAM_gsize-hs_no_control/ln/alias/experiments/hg38_H3K27ac_thymus/input_peaks.broadPeak"
                            path_broadPeak_in = "out/macs2/callpeak_broad_format-AUTO_gsize-" + gsize + "_no_control/" + in_aligned_suffix + "_peaks.broadPeak" 
                            path_broadPeak_ln = ln_processed_prefix + "/peaks/broad/" + row['sample_name'] + ".bed"
                            d_paths_in_to_ln[path_broadPeak_in] = path_broadPeak_ln

                            path_narrowPeak_in = "out/macs2/callpeak_format-AUTO_gsize-" + gsize + "_no_control/" + in_aligned_suffix + "_peaks.narrowPeak" 
                            path_narrowPeak_ln = ln_processed_prefix + "/peaks/narrow/" + row['sample_name'] + ".bed"
                            d_paths_in_to_ln[path_narrowPeak_in] = path_narrowPeak_ln
                            
                            #if (row['control_name'] not in ["",'irrelevant']) and (row['type'] in ['ChIP']):
                            #print('debug if row')
                            if row['control_name'] not in ["",'irrelevant']:
                                print('TODO: Add to meta.rules to process peaks with control by input for sample: ' +  row['sample_name'])
                                # out/ln/sst_exp/ChIP/mm9/P5424_WT_K27ac.bam
                                path_broadPeak_in = "out/macs2/callpeak_broad_format-AUTO_gsize-" + gsize + "/ln/sst_exp/" + row['type'] + "/" + assembly + "/" + row['sample_name'] + "_over_" + row['control_name'] + "_peaks.broadPeak"
                                path_broadPeak_ln = ln_processed_prefix + "/peaks/broad/" + row['sample_name'] + "_VS_" + row['control_name'] + ".bed"
                                d_paths_in_to_ln[path_broadPeak_in] = path_broadPeak_ln


                        # TODO: Check first here that control_name is a valid one and store its row['tgml_fastq_prefix'].
                        # Note: Maybe the easiest to deal with such case is to keep all samples in a pooled directory so all of them can be used in pattern like "CHIP VS CONTROL".
                        #if row['control_name'] != "":
                        #    print('dev control for sample: ' + data['control_name'])
                        #    EXTENDREADS = ""
                        #    path_bw_ratio_in="out/deepTools/bamCompare_scaleFactorsMethod-SES_ratio-log2_binSize-20_extendReads-" + EXTENDREADS + "/" + in_aligned_suffix + "_vs_samtools/sort/samtools/view_bSh/bowtie2/" + row['se_or_pe'] + "_" + assembly +"/sickle/" + row['se_or_pe'] + "_-t_sanger_-q_20/" + "gunzip/merge_lanes_nextseq500_" + row['se_or_pe'] + "_raw/" + row['tgml_fastq_prefix']
                        #    {id_bam2}.bw"

                        if row['type'] == 'RNA':
                            path_bw_forward_in="out/deepTools/bamCoverage_--filterRNAstrand_forward_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM/" + in_aligned_suffix + ".bw"
                            path_bw_forward_ln = ln_processed_prefix + "/bw/forward/" + row['sample_name'] + ".bw"
                            d_paths_in_to_ln[path_bw_forward_in] = path_bw_forward_ln
        
                            path_bw_reverse_in="out/deepTools/bamCoverage_--filterRNAstrand_reverse_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM/" + in_aligned_suffix + ".bw"
                            path_bw_reverse_ln= ln_processed_prefix + "/bw/reverse/" + row['sample_name'] + ".bw"
                            d_paths_in_to_ln[path_bw_reverse_in] = path_bw_reverse_ln

                        #if row['type'] == 'ChIP':
                            #path_peaks
        

        # Dealing with samples and merge of samples gathered into an experiment:
        # Currently working on it for Capstarrseq:
        if row['exp_name'] != '':
            print('sample_name "' + row['sample_name'] + '" in exp_name "' + row['exp_name'] +'"')
            if row['type'] == 'CapStarr' and row['specie'] == 'mouse':
                path_bam_in = "out/ln/sst_exp/" + row['exp_name'] + "/mm9/" + row['sample_name'] + ".bam"
                path_bam_ln = "sst/experiments/" + row['exp_name'] + "/mm9/" + row['sample_name'] + ".bam"
                d_paths_in_to_ln[path_bam_in] = path_bam_ln
                
                path_bai_in = path_bam_in + ".bai"
                path_bai_ln = path_bam_ln + ".bai"
                d_paths_in_to_ln[path_bai_in] = path_bai_ln
        
                path_bw_in="out/deepTools/bamCoverage_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM_--extendReads_314/" + in_aligned_suffix + ".bw"
                path_bw_ln = "sst/experiments/" + row['exp_name'] + "/mm9/" + row['sample_name'] + ".bw"
                d_paths_in_to_ln[path_bw_in] = path_bw_ln
 
                if row['control_name'] not in ['','irrelevant']:
                    path_tsv_in = 'out/capstarrseq/merge_all_data_mTDHS/awk/extend_reads_314/awk/keep_first_mate_for_pe_bedtools_bamtobed/bedtools/bamtobed/ln/sst_exp/' + row['exp_name'] + "/mm9/" + row['sample_name'] + "_over_" + row['control_name'] + ".allData.tsv"
                    path_tsv_ln = "sst/experiments/" + row['exp_name'] + "/mm9/" + row['sample_name'] + "_over_" + row['control_name'] + ".allData.tsv"
                    d_paths_in_to_ln[path_tsv_in] = path_tsv_ln

                    path_pdf_in = 'out/capstarrseq/grouping_crms/capstarrseq/fold_change_mTDHS/awk/extend_reads_314/awk/keep_first_mate_for_pe_bedtools_bamtobed/bedtools/bamtobed/ln/sst_exp/' + row['exp_name'] + "/mm9/" + row['sample_name'] + "_over_" + row['control_name'] + ".inflexionPointGroups.pdf"
                    path_pdf_ln = "sst/experiments/" + row['exp_name'] + "/mm9/" + row['sample_name'] + "_over_" + row['control_name'] + ".inflexionPointGroups.pdf"
                    d_paths_in_to_ln[path_pdf_in] = path_pdf_ln

    # This file will log links done between mw and sst
    if usage == 'run':
        file = open("out/ln/sst_table.tsv","w")
        for path_in in d_paths_in_to_ln.keys():
            path_ln = d_paths_in_to_ln[path_in]
            shell("mkdir -p `dirname {path_ln}`; ln -f {path_in} {path_ln}")
            file.write(path_in + "\t" + path_ln)
        file.close()
    
    return list(d_paths_in_to_ln.keys())

def input_ln_sst_table(wildcards):
    """
    Created:
        2017-10-03 10:42:05
    Aim:
        Map files to aliases used by the rule ln_alias.
    """
    #alias_id = wildcards['alias_id']
    #paths = mapper_mw_to_sst('input', tsv=input['tsv'])
    #paths = mapper_mw_to_sst('input', tsv='src/snakemake/tables/salva_runs.tsv')
    paths = mapper_mw_to_sst('input', tsv=TSV)
    return paths


## function to create an ID table for bam to merge:
#for row in csv.DictReader(open(TSV),delimiter='\t'):
#    paths = list()
#    print('Creating IDs for samples to merge')
#    if row['sample_merge_list'] != "":
#        samples = row['sample_merge_list'].split(",")
#        for sample in samples:
#            for row2 in csv.DictReader(open(TSV),delimiter='\t'):
#                if row2['sample_name'] == sample:
#                paths.append(mwconf[sample]['path_bam_in'])
#    return(paths)

def input_bam_samtools_merge_samples_sst(wildcards):
    """
    Created:
        2018-01-01 22:29:32
    Aim:
    """
    paths = list()
    assembly = wildcards['assembly']

    for row in csv.DictReader(open(TSV),delimiter='\t'):

        # Only look for samples that have available files
        #if ',' in row['sample_merge_list']:
        #if row['exp_name'] == wildcards['exp_name']:···
        if row['sample_name'] == wildcards['sample_name']:
            if row['exp_name'] != wildcards['exp_name']:
                raise Exception("This sample can't be in this experiment based on metadata.")
            #print("sample_name : " + row['sample_name'])
            samples = row['sample_merge_list'].split(",")
            for sample in samples:
                #print('sample: ' + sample)
                for row2 in csv.DictReader(open(TSV),delimiter='\t'):
                    if row2['sample_name'] == sample:
                        path_bam_in = "out/samtools/sort/samtools/view_bSh/bowtie2/" + row2['se_or_pe'] + "_" + assembly +"/sickle/" + row2['se_or_pe'] + "_-t_sanger_-q_20/" + "gunzip/merge_lanes_nextseq500_" + row2['se_or_pe'] + "_raw/ln/updir/mw/" + row2['tgml_fastq_prefix'] + ".bam"
                        print('bam: ' + path_bam_in)
                        paths.append(path_bam_in)

    return(paths)

def input_bam_ln_sst_exp(wildcards):
    """
    Created:
        2018-01-29 14:54:54
    Aim:
    """
    assembly = wildcards['assembly']
    exp_name = wildcards['exp_name']
    sample_name = wildcards['sample_name']

    for row in csv.DictReader(open(TSV),delimiter='\t'):
        if row['sample_name'] == wildcards['sample_name']:
            if exp_name not in [row['exp_name'],'ChIP']:
                # 'ChIP' added here because I want all ChIP-Seq samples to be put together in a 'ChIP' type experiment, mainly because it is easier then to deal with peak callin with control samples.
                raise Exception("This sample can't be in this experiment based on metadata.")
            if ',' in row['sample_merge_list']:
                path = "out/samtools/merge_samples_sst/" + exp_name + "/" + assembly + "/" + sample_name + ".bam"
            else:
                path = "out/samtools/sort/samtools/view_bSh/bowtie2/" + row['se_or_pe'] + "_" + assembly +"/sickle/" + row['se_or_pe'] + "_-t_sanger_-q_20/" + "gunzip/merge_lanes_nextseq500_" + row['se_or_pe'] + "_raw/ln/updir/mw/" + row['tgml_fastq_prefix'] + ".bam"

    return(path)


