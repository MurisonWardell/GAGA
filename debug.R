## Testing that gaga works

setwd("C:/Users/cwardell/Desktop/temp")

## Load library
library(GAGA)

## Load synthetic data set
#data("gaga_synthetic_data","gaga_synthetic_data_annotation")
data("hidden_gaga")
data("gaga_simple_data")
data("gaga_synthetic_data")
data("gaga_synthetic_data_jittered")
data("BYB1_G07_pruned")

## Execute gaga - good data, don't mess with it...
#x=gaga(gaga_synthetic_data, gaga_synthetic_data_annotation, number_of_clones=6, iterations=20, contamination=0)
#y=gaga(gaga_synthetic_data, gaga_synthetic_data_annotation, number_of_clones=6, iterations=10, contamination=0)
gdata=gaga(gaga_synthetic_data, number_of_clones=6, iterations=3000)
gdataj=gaga(gaga_synthetic_data_jittered, number_of_clones=6, iterations=3000)
gagaYeast=gaga(BYB1_G07_pruned, number_of_clones=6, iterations=1000)
gagaYeastC=gaga(BYB1_G07_pruned, number_of_clones=6, iterations=1000,contamination=1)

## Correct solution for 
simpleDataSolution=gaga(gaga_simple_data, number_of_clones=4, nroot=0,iterations=3000)
simpleDataSolution=gaga(gaga_simple_data, number_of_clones=4, nroot=1,iterations=500) # will hit correct solution faster


gagaReport(gaga_synthetic_data,gdata,output_file_prefix="gaga_synthetic_data")
gagaReport(gaga_synthetic_data_jittered,gdataj,output_file_prefix="gaga_synthetic_data_jittered")
gagaReport(BYB1_G07_pruned,gagaYeast,output_file_prefix="gagaYeast")
gagaReport(BYB1_G07_pruned,gagaYeastC,output_file_prefix="gagaYeastC")

#x=BYB1_G07_pruned
#y=gagaYeast
#y=gagaYeastC
gagaReport(x,y,outType="fitness")
gagaReport(x,y,outType="heatmap")
gagaReport(x,y,outType="proportion")
gagaReport(x,y,outType="phylogeny")
gagaReport(x,y)


gagaReport(gaga_synthetic_data,x,outType="complete")
gagaReport(gaga_synthetic_data,gdata,outType="fitness")
gagaReport(gaga_synthetic_data,z,outType="heatmap")
gagaReport(gaga_synthetic_data,z,outType="proportion")
gagaReport(gaga_synthetic_data,z,outType="phylogeny")

setwd("N:/MORGAN/ChrisW/documents/projects/131213_gaga/GAGA")


## The horrifying building vignettes section:
library(devtools)
build_vignettes()


