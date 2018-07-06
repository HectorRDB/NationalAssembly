# Load data
library(tidyverse)
voting_record <- read_csv("data/voting_record.csv")
meta <- voting_record %>% select(circo, dept, name, surname, identifiant,
                                 iden, chamber, NbVote, NbYes, NbNo, NbAbst)
voting_record <- voting_record %>% select(-one_of(colnames(meta))) %>%
                                   t(.) 
scrutins <- read_csv("data/scrutins.csv")
colnames(voting_record) <- rownames(meta) <- meta$identifiant
rownames(voting_record) <- rownames(scrutins) <- scrutins$Number
library(Biobase)
Assay <- ExpressionSet(assayData = voting_record, 
                       phenoData = AnnotatedDataFrame(meta) ,
                       featureData = AnnotatedDataFrame(scrutins))
