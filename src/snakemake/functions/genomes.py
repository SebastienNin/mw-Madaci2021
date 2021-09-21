# Method to define genomes from a genomes sheet in Sequencing_summary.xlsx

#if 'ids' not in config:
#    mwconf['ids'] = {}
#
#genomes = pandas.read_excel("../mw-sst/Sequencing_summary.xlsx", sheet_name="genomes")
#
#for index, row in genomes.iterrows():
#    ASSEMBLY = str(row['assembly'])
#    ANNOTATION = str(row['annotation'])
#    PROCESS = str(row['process'])
#
#    if PROCESS == 'yes':
#        GENOME = ASSEMBLY + '-' + ANNOTATION
#
#        mwconf['ids']['fa-genome-'+ ASSEMBLY] = str([str(row['fasta'])])
#        print(mwconf['ids']['fa-genome-'+ ASSEMBLY])
#        mwconf['ids']['gtf-'+ GENOME] = str(row['gtf'])
#        mwconf['ids']['star-idx-' + GENOME] = "out/star/build_index/fa-genome-" + GENOME + "_gtf-" + GENOME + "/Genome" 
#        #mwconf['ids']['bowtie2-idx-' + GENOME] = ["out/bowtie2-build/" + str(row['fasta']) + "." + str(parts) + ".bt2" for parts in ["1","2","3","4","rev.1","rev.2"]]
#        mwconf['ids']['bowtie2-idx-' + ASSEMBLY] = ["out/bowtie2-build/fa-genome-" + ASSEMBLY + "." + str(parts) + ".bt2" for parts in ["1","2","3","4","rev.1","rev.2"]]
#        mwconf['ids']['bwa-idx-' + ASSEMBLY] = ["out/bwa/index/fa-genome-" + ASSEMBLY + "." + str(parts) for parts in ["amb","ann","bwt","pac","sa"]]
#
#        #out/bowtie2-build/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.1.bt2
#
#

#genomes.loc[((genomes['vernacular_name'] == specie) | (genomes['specie_name'] == specie) | (genomes['assembly_id'] == specie)) & (genomes['process'] == 'no') ]
#out/star/build_index/fa-GRCh38_chr22_gtf-GRCh38_chr22/Genome

