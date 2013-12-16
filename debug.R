## Testing that gaga works

## Load library
library(GAGA)

## Load synthetic data set
data("gaga_synthetic_data","gaga_synthetic_data_annotation")

## Execute gaga - good data, don't mess with it...
#x=gaga(gaga_synthetic_data, gaga_synthetic_data_annotation, number_of_clones=6, iterations=20, contamination=0)


gagaReport(gaga_synthetic_data,x)


