##### cellHTS2 install #####

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("cellHTS2")
library(cellHTS2)

# namespace ‘vctrs’ 0.4.2 is already loaded, but >= 0.5.0 is required 
# ^ if this error shows up, uninstall package and re-install

##### Importing data #####

setwd("C:/Users/misak/Desktop/ICR/siRNAscreen/siRNA_screen_validation_scripts/datafiles") 
  # siRNAscreen is on GitHub repository

# create cellHTS object calld "x"
x <- readPlateList(
  filename="Platelist.txt",
  name="SUM149viab",
  path="./"
)

# add description and plate layout to x
x <- configure(
  x,
  descripFile="Description.txt",
  confFile="Plateconfig.txt",
  path="./"
)

# add information about genes targeted by siRNA
x <- annotate(
  x,
  geneIDFile="Genelibrary.txt",
  path="./"
) #  if there is an error message, close R and reload cellHTS2

##### Normalization & z-score #####

# Perform normalisation of each plate in the screen. 
# Log transform and center data to the median intensity of each plate
xn <- normalizePlates(
  x,
  scale="multiplicative",
  log=TRUE,
  method="median",
  varianceAdjust = "none",
  posControls="siPLK1",
  negControls="siCON1|siCON2|allstar"
)

# Scale the well intensities to the median absolute
# deviation across the library to produce Z-scores
xsc <- scoreReplicates(
  xn,
  method="zscore",
  sign="+"
)

# Write a table with z scores for all the replicates before summarizing
genes <- geneAnno(xsc)
plates <-plate(xsc)
wells <- well(xsc)
scores <- Data(xsc)
all_z <- data.frame(
  gene=genes,
  plate=plates,
  well=wells,
  zscore=scores
)

write.table(
  all_z,
  "all_zscores.txt",
  col.names=TRUE,
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)


# Summarise the Z-scores as medians
xsc <- summarizeReplicates(
  xsc,
  summary="median"
)


##### Create report #####

# configure the HTML report created by writeReport()
setSettings(
  list(
    plateList=list(
      reproducibility=list(
        include=TRUE,
        map=TRUE
      ),
      intensities=list(
        include=TRUE,
        map=TRUE)
    ),
    screenSummary=list(
      scores=list(
        range=c(-20, 10),
        map=TRUE
      )
    )
  )
)

# Produce a detailed HTML report in a subdirectory
writeReport(
  raw=x,
  normalized=xn,
  scored=xsc,
  outdir="./report",
  force=TRUE,
  posControls="siplk1",
  negControls="sicon1|sicon2|allstar", # with the original dataset, putting posControls and negControls in lowercase works
  mainScriptFile="C:/Users/misak/Desktop/ICR/siRNA/get_zscore.R"
)

# For convenience, write Z-scores to a file that can be
# joined with further screens for analysis 
genes <- geneAnno(xsc)
plates <-plate(xsc)
wells <- well(xsc)
scores <- Data(xsc)[,1,1]
combinedz <- data.frame(
  gene=genes,
  plate=plates,
  well=wells,
  zscore=scores
)


write.table(
  combinedz,
  "summarized_zscores.txt",
  col.names=TRUE,
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)



# # normalized values
# genes <- geneAnno(xn)
# plates <-plate(xn)
# wells <- well(xn)
# values <- Data(xn)[,1,1]
# combinednorm <- data.frame(
#   gene=genes,
#   plate=plates,
#   well=wells,
#   normalized=values
# )
# write.table(
#   combinednorm,
#   "normalized_values.txt",
#   col.names=TRUE,
#   sep="\t",
#   quote=FALSE,
#   row.names=FALSE
# )
