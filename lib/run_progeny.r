library(progeny)
library(dplyr)
library(tibble)
library(ggplot2)

message("STARTING PROGENY PROCESS")

message("Reading input arguments:")
args = commandArgs(trailingOnly = TRUE)
progeny_file = args[1]
organism = args[2] # "Human"
zscores = as.logical(args[3]) # T
top = as.numeric(args[4]) # 100

message(progeny_file)
message(organism)
message(zscores)
message(top)

message("Creating output file")

file_csv = paste0("progeny_scores_", organism, "_", top, ".csv")
message(file_csv)

message("Reading input file")

progeny_data <- as.matrix(read.csv(progeny_file, row.names = 1))

message("Running Progeny") 

PathwayActivity_counts <- progeny::progeny(progeny_data, 
                                           scale = TRUE, 
                                           organism = organism, 
                                           top = top, 
                                           perm = 10000, 
                                           z_scores = zscores)

PathwayActivity_counts <- as.data.frame(t(PathwayActivity_counts))

message("Writing output file")

write.csv(PathwayActivity_counts, file_csv, quote = F)

message("FINISHING PROGENY PROCESS")
