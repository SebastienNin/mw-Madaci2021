cd /gpfs/tgml/nin/mw-Madaci2021

mkdir -p out/cite-seq_count/Run_342/

cp out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_--expected_cells_10000_-o_Results_HTO/cat/merge-nexsteq500-pe/S003960_HTO_Hashtag_R342/tags.csv out/cite-seq_count/Run_342/tags.csv

cp out/cellranger/count_--sample_mRNA_Hashtag_R342_GRCh38-2020-A/cellranger/mkfastq/gpfs/tgml/reads/NS500_output_from_sept_2020/Run_342_NS500-249_11-12-2020_LM/Run_342/outs/filtered_feature_bc_matrix/barcodes.tsv.gz out/cite-seq_count/Run_342/barcodes.tsv.gz
gunzip out/cite-seq_count/Run_342/barcodes.tsv.gz

CELLNB=`wc -l out/cite-seq_count/Run_342/barcodes.tsv | cut -d " " -f 1`
sed 's/-1//g' out/cite-seq_count/Run_342/barcodes.tsv > out/cite-seq_count/Run_342/barcodes_gsub.tsv
echo $CELLNB
source activate cite-seq-count
CITE-seq-Count -cbf 1 -cbl 16 -umif 17 -umil 28 --whitelist out/cite-seq_count/Run_342/barcodes_gsub.tsv --expected_cells $CELLNB -o out/cite-seq_count/Run_342/ -R1 out/cat/merge-nexsteq500-pe/S003960_HTO_Hashtag_R342_1.fastq.gz -R2 out/cat/merge-nexsteq500-pe/S003960_HTO_Hashtag_R342_2.fastq.gz -t out/cite-seq_count/Run_342/tags.csv

