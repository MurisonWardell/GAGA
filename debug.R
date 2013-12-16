## Testing that gaga works
library(GAGA)

#gaga("data/synthetic_data_90.txt", "xxx", "data/synthetic_data_annotation.txt", 6, 6, iterations=3, contamination=0)

data("gaga_synthetic_data","gaga_synthetic_data_annotation")

x=gaga(gaga_synthetic_data, gaga_synthetic_data_annotation, 6, 6, iterations=3, contamination=0)



