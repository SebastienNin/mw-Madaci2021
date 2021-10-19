# Based on rocker - https://github.com/rocker-org/rocker-versioned
FROM rocker/tidyverse:4.1.1

MAINTAINER SNin (sebastien.nin@inserm.fr)

# ##################
# OTHER DEPENDENCIES
# ##################

#### Python
# pip
RUN apt-get update && apt-get install -y \
    python3-pip \
    && pip3 install --upgrade --user pip

# Umap processing library (python)
RUN python3 -m pip install umap-learn

# Install libglpk40 (required for igraph)
RUN apt-get update && apt-get install -y \
  build-essential \
  libglpk40

# ##########
# R PACKAGES 
# ##########

#### Figures & layout management
# ggplot2
RUN Rscript -e 'install.packages( "ggplot2")'
RUN Rscript -e 'install.packages( "cowplot")'        # plot_grid, themes, ...
RUN Rscript -e 'install.packages( "ggpubr")'         # add_summary, geom_signif, ...
RUN Rscript -e 'install.packages( "ggrepel")'        # geom_text_repel, geom_label_repel
RUN Rscript -e 'install.packages( "gridExtra")'      # grid.arrange, ...
RUN Rscript -e 'BiocManager::install( "patchwork")'  # +/ operators for ggplots

# plotly
RUN Rscript -e 'install.packages( "plotly")'

# general
RUN Rscript -e 'install.packages( "gplots")'         # heatmap.2
RUN Rscript -e 'install.packages( "heatmaply")'      # heatmaply (interactive)
RUN Rscript -e 'BiocManager::install( "iheatmapr")'  # iheatmap (interactive, uses plotly), dependencies OK with BiocManager
RUN Rscript -e 'install.packages( "pheatmap")'       # pheatmap
RUN Rscript -e 'install.packages( "data.table")'
RUN Rscript -e 'install.packages( "factoextra")'

#### Reporting
RUN Rscript -e 'install.packages( "DT")'             # datatable
RUN Rscript -e 'install.packages( "pander")'         # pander
RUN Rscript -e 'install.packages( "bookdown")'       # bookdown

#### General
RUN Rscript -e 'install.packages( "funr")'           # get_script_path
RUN Rscript -e 'install.packages( "reshape2")'        # melt


#### Technology specific
RUN Rscript -e 'BiocManager::install( "Seurat")'     # Dependencies OK with BiocManager (https://github.com/satijalab/seurat/issues/2409)
RUN Rscript -e 'install.packages( "umap")'

#### Add package for better heatmap
RUN Rscript -e 'BiocManager::install( "ComplexHeatmap")' 

RUN Rscript -e 'BiocManager::install( c("Biobase", "limma"))'

#### To automatically generate a R documentation for each function
RUN Rscript -e 'install.packages( "roxygen2")'

#### To test each function
RUN Rscript -e 'install.packages( "testthat")'

#### Add remark for syntax, style,... of R script
RUN Rscript -e 'install.packages( "lintr")'

#### Add remark for syntax, style,... of R script
RUN Rscript -e 'install.packages( "styler")'


#### Custom
#RUN Rscript -e 'install.packages( "")'

# #################################################
# DECLARE EXPECTED PORTS TO BE EXPOSED
# #################################################

# Rstudio listens on port 8787 by default
# Container typically started with `-p 8787:8787`
EXPOSE 8787





