library(here)

# Directori arrel
rootDir <- here()

# Primer nivell
dataDir <- here(rootDir, "data")
analysisDir <- here(rootDir, "analysis")
RscriptsDir <- here(rootDir, "Rscripts")

# Segon nivell
rawdataDir <- here(rootDir, "data", "rawdata")
cleandataDir <- here(rootDir, "data","cleandata")