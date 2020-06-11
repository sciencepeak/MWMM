#!/bin/bash

time Rscript 03.GO_similarity.01.calculation.R BLCA &
time Rscript 03.GO_similarity.01.calculation.R BRCA &
time Rscript 03.GO_similarity.01.calculation.R COAD &
time Rscript 03.GO_similarity.01.calculation.R ESCA &
time Rscript 03.GO_similarity.01.calculation.R HNSC &
time Rscript 03.GO_similarity.01.calculation.R KICH &
time Rscript 03.GO_similarity.01.calculation.R KIRC &
time Rscript 03.GO_similarity.01.calculation.R KIRP &
time Rscript 03.GO_similarity.01.calculation.R LIHC &
time Rscript 03.GO_similarity.01.calculation.R LUAD &
time Rscript 03.GO_similarity.01.calculation.R LUSC &
time Rscript 03.GO_similarity.01.calculation.R PRAD &
time Rscript 03.GO_similarity.01.calculation.R STAD &
time Rscript 03.GO_similarity.01.calculation.R THCA &
time Rscript 03.GO_similarity.01.calculation.R UCEC &

wait
echo "all are finished"
