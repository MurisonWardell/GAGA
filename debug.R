## Testing that gaga works

setwd("C:/Users/cwardell/Desktop")

## Load library
library(GAGA)

## Load synthetic data set
#data("gaga_synthetic_data","gaga_synthetic_data_annotation")
data("new_gaga_synthetic_data")

## Execute gaga - good data, don't mess with it...
#x=gaga(gaga_synthetic_data, gaga_synthetic_data_annotation, number_of_clones=6, iterations=20, contamination=0)
#y=gaga(gaga_synthetic_data, gaga_synthetic_data_annotation, number_of_clones=6, iterations=10, contamination=0)
z=gaga(new_gaga_synthetic_data, number_of_clones=6, iterations=10, contamination=0)


gagaReport(gaga_synthetic_data,x,outType="complete")
gagaReport(gaga_synthetic_data,z,outType="fitness")
gagaReport(gaga_synthetic_data,z,outType="heatmap")
gagaReport(gaga_synthetic_data,z,outType="proportion")
gagaReport(gaga_synthetic_data,z,outType="phylogeny")

setwd("N:/MORGAN/ChrisW/documents/projects/131213_gaga/GAGA")
