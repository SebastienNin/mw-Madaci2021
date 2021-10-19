cd /gpfs/tgml/nin/mw-Madaci2021

cat out/cat/merge-nexsteq500-pe/S003718_19_3_HTO_Hashtag_R306_1.fastq.gz out/cat/merge-nexsteq500-pe/S003718_19_1_HTO_Hashtag_R309_1.fastq.gz > out/cat/merge-nexsteq500-pe/HTO_R306_R309_1.fastq.gz
cat out/cat/merge-nexsteq500-pe/S003718_19_3_HTO_Hashtag_R306_2.fastq.gz out/cat/merge-nexsteq500-pe/S003718_19_1_HTO_Hashtag_R309_2.fastq.gz > out/cat/merge-nexsteq500-pe/HTO_R306_R309_2.fastq.gz

mkdir -p out/cite-seq_count/Run_306_Run_309/

cp out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_--expected_cells_10000_-o_Results_HTO/cat/merge-nexsteq500-pe/S003960_HTO_Hashtag_R342/tags.csv out/cite-seq_count/Run_306_Run_309/tags.csv

cp out/cellranger/count_mRNA_R306_R309/outs/filtered_feature_bc_matrix/barcodes.tsv.gz out/cite-seq_count/Run_306_Run_309/barcodes.tsv.gz
gunzip out/cite-seq_count/Run_306_Run_309/barcodes.tsv.gz

CELLNB=`wc -l out/cite-seq_count/Run_306_Run_309/barcodes.tsv | cut -d " " -f 1`
sed 's/-1//g' out/cite-seq_count/Run_306_Run_309/barcodes.tsv > out/cite-seq_count/Run_306_Run_309/barcodes_gsub.tsv

source activate cite-seq-count
#CITE-seq-Count -cbf 1 -cbl 16 -umif 17 -umil 28 --expected_cells 20000 -o out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_--expected_cells_20000_-o_Results_HTO/cat/merge-nexsteq500-pe/HTO_R306_R309/ -R1 out/cat/merge-nexsteq500-pe/HTO_R306_R309_1.fastq.gz -R2 out/cat/merge-nexsteq500-pe/HTO_R306_R309_2.fastq.gz -t out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_--expected_cells_20000_-o_Results_HTO/cat/merge-nexsteq500-pe/HTO_R306_R309/tags.csv
CITE-seq-Count -cbf 1 -cbl 16 -umif 17 -umil 28 -cells $CELLNB -o out/cite-seq_count/Run_306_Run_309/ -R1 out/cat/merge-nexsteq500-pe/HTO_R306_R309_1.fastq.gz -R2 out/cat/merge-nexsteq500-pe/HTO_R306_R309_2.fastq.gz -t out/cite-seq_count/Run_306_Run_309/tags.csv -wl out/cite-seq_count/Run_306_Run_309/barcodes_gsub.tsv

