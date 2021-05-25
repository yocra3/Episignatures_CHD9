# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_12

# Install required Bioconductor packages
RUN R -e 'BiocManager::install(c("minfi", "Gviz", "a4Base", "tidyverse", "limma", "e1071", "sva", "cowplot", "dendextend", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))'

# Install github packages (last commit when creating docker)
RUN R -e 'devtools::source_url("https://raw.githubusercontent.com/isglobal-brge/EASIER/HEAD/installer.R")' ## EASIER - for Enrichments
RUN R -e 'devtools::install_github("perishky/meffil")' ## meffil
