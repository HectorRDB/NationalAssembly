Understanding RNA-seq analysis through an example: voting behavior at the French parliament
================
Hector Roux de Bézieux
May 2018

-   [Introduction](#introduction)
-   [Context](#context)
-   [Loading the data and building an ExpressionSet object](#loading-the-data-and-building-an-expressionset-object)
-   [Data Cleaning](#data-cleaning)
    -   [Filtering the cells (delegate)](#filtering-the-cells-delegate)
    -   [Filtering the genes (votes)](#filtering-the-genes-votes)
-   [Normalization and correcting for zero inflation](#normalization-and-correcting-for-zero-inflation)
-   [Recovering / discovering biological types](#recovering-discovering-biological-types)
    -   [Clustering](#clustering)
    -   [Finding differentially expressed genes](#finding-differentially-expressed-genes)
-   [R session](#r-session)
-   [Thanks](#thanks)

Introduction
============

Single-cell RNA-seq (scRNA-seq) is a very potent biological tool used for many applications. A far from exhaustive list would include identifying new cell types, finding differentially expressed (DE) genes, and discovering lineages among cells. However, the usual framework might seem a little daunting for beginners and, while many well-crafted tutorials exists, they all share the same idea: use a biological dataset as an example. Here, we want to use a dataset that would be more understandable to a broader public to explain the usual steps in scRNA-seq analysis. Only some basics on data analysis are needed to understand this tutorial. Knowledge of R helps to understand the code but is not necessary to follow along.

Context
=======

Our dataset is the voting record of the French delegates of the 14<sup>*t**h*</sup> legislature, going from May 2012 to May 2017. Delegates are characterized by their names and surnames, their department and circonscription, their political group. For each of the 644 votes (which can be bills, amendments, choosing the prime minister, ...), we also have the voting behavior of each delegate (voted yes, no, abstain or did not took part in the vote) and the reason of the vote. We are interested in investigating the link between voting behavior and political affiliations.

Loading the data and building an ExpressionSet object
=====================================================

The usual scRNA-seq data is a matrix of *J* genes by *n* biological cells. Each cell of the matrix represents the number of copies of messenger RNA of a given gene in a given biological cell.This number serves as a proxy for the number of proteins from this mRNA and therefore represent how active a gene is in a given cell. We also eventually have some information about cells

Here, we instead have a *J* bills by *n* delegates matrix. Each cell represent the vote of the delegate for that low (−1 if against, 1 is for, 0 if not voting or abstained). We store the data in an *ExpressionSet* object, where the *phenodata* is the information about the delegates (the cells) and the *featuredata* is the information about the votes (the genes).

``` r
# Create the phenodata matrix
voting_record <- read_csv("data/voting_record.csv")
meta <- voting_record %>% select(circo, dept, name, surname, identifiant, iden,
                                 group1, group2, group3, NbVote, NbYes, NbNo, NbAbst)
voting_record <- voting_record %>% select(-one_of(colnames(meta))) %>%
                                   t(.)

# We count abstaining and not voting as the same thing for simplification
voting_record[is.na(voting_record)] <- "NA"
voting_record[voting_record == "For"] <- "1"
voting_record[voting_record == "Against"] <- "-1"
voting_record[voting_record == "Abstain"] <- "0"
voting_record[voting_record == "NotVoting"] <- "0"
voting_record[voting_record == "NA"] <- NA
voting_record <- apply(voting_record, 2, as.numeric)

# Build the ExpressionSet object
votes <- read_csv("data/votes.csv")
colnames(voting_record) <- rownames(meta) <- meta$identifiant
rownames(voting_record) <- rownames(votes) <- votes$Number
Assay <- ExpressionSet(assayData = voting_record, 
                       phenoData = AnnotatedDataFrame(meta) ,
                       featureData = AnnotatedDataFrame(votes))
```

Data Cleaning
=============

Filtering the cells (delegate)
------------------------------

In scRNA-seq settings, we want to filter the cells with too few total reads, or with too many reads of poor quality (that do not match to the genome of reference). Here, we instead filter the delegates with too few voting possibilities. This may happen either because a delegate had to leave office before the end of the session, or was appointed to the government, and then its replacement is not around long enough.

We therefore filter delegates whose number of votes is more than 2 SD below the mean (so equal or less to 40).

``` r
Nb_Votes <- data.frame(Nb_Votes = colSums(!is.na(exprs(Assay))),
                       delegate = rownames(phenoData(Assay)))

ggplot(Nb_Votes, aes(x = Nb_Votes, fill = Nb_Votes > 40)) +
      stat_bin(breaks = seq(0, max(Nb_Votes$Nb_Votes), 20)) +
      labs(x = "Number of votes") +
      guides(fill = guide_legend(title = "Could vote\nmore than\n40 times")) +
      scale_fill_manual(values = brewer.pal(3, "Set1")[c(1, 3)],
                        breaks = c(TRUE, FALSE)) +
  theme_classic()
```

<img src="Single-cell_files/figure-markdown_github/filtering cells-1.png"  style="display: block; margin: auto;" />

``` r
# Filter
Assay <- Assay[, Nb_Votes$Nb_Votes > 40]
```

We then simplify the data by replacing NAs by 0 (as abstentions). We can then check whether some delegates rarely voted. Again, we remove people whose number of votes is more than 2 SD below the mean number of votes (so delegates who actually voted less than 38 votes). Among those, we of course find the president of the parliament, who does not take part in the votes and therefore never voted.

``` r
# Remplace NAs by zeros
E <- exprs(Assay)
E[is.na(E)] <- 0
pd <- phenoData(Assay)
fd <- featureData(Assay)
colnames(E) <- rownames(pd) <- pd@data$identifiant
rownames(E) <- rownames(fd) <- fd@data$Number
Assay <- ExpressionSet(assayData = E, phenoData = pd, featureData = fd)

# Filter
Nb_Votes <- data.frame(Nb_Votes = colSums(exprs(Assay) != 0),
                       delegate = rownames(phenoData(Assay)))
ggplot(Nb_Votes, aes(x = Nb_Votes, fill = Nb_Votes > 38)) +
      stat_bin(breaks = seq(0, max(Nb_Votes$Nb_Votes), 19)) +
      labs(x = "Number of votes") +
  guides(fill = guide_legend(title = "Actually voted\nmore than\n38 times")) +
  scale_fill_manual(values = brewer.pal(3, "Set1")[c(1, 3)],
                        breaks = c(TRUE, FALSE)) +
  theme_classic()
```

<img src="Single-cell_files/figure-markdown_github/change NA to zeros-1.png"  style="display: block; margin: auto;" />

``` r
Assay <- Assay[, Nb_Votes$Nb_Votes > 38]
```

Filtering the genes (votes)
---------------------------

In scRNA-seq settings, we want to filter out genes that are expressed in a small number of cells. This is done for two reasons. The first is computational. Most datasets would have dozens of thousands of genes, over dozens or hundreds of cells. Reducing the number of genes helps speed up calculations. Secondly, we are usually interested in a specific tissue or process. Most genes are not being expressed in the samples of interest and just add noise to the analysis. Hence filtration. Simple heuristics in the form of "at least *i* counts in *j* cells" are usually quite efficient. The aim is too be quite broad. We want to filter as many genes as possible while retaining most of the reads. In usual settings, you would only keep ~20% of the genes at most but that would still capture more than 90% of the reads.

In our context, we only have 644 votes so the computational goal does not really hold. But we can filter out votes where many delegates did not vote. Those are probably non-political technical votes. We therefore use the filter rules "at least *j* delegates cast a vote on the bill".

``` r
# Compute the number of votes and number of vote casted as a function of j
Nb_Votes <- rowSums(exprs(Assay) != 0)
voting_behavior <- data.frame(n_votes_cast = rep(0, nrow(Assay)),
                             n_votes = rep(0, nrow(Assay)),
                             j = 1:nrow(Assay))
for(j in 1:nrow(Assay)){
  Nb_Votes <- Nb_Votes[Nb_Votes >= j]
  voting_behavior$n_votes[j] <- length(Nb_Votes)
  voting_behavior$n_votes_cast[j] <- sum(Nb_Votes)
}

voting_behavior <- voting_behavior %>% mutate(n_votes_cast = n_votes_cast /
                                                              max(n_votes_cast),
                                              n_votes = n_votes / max(n_votes))

ggplot(voting_behavior,
       aes(x = n_votes, y = n_votes_cast, col = j)) +
  geom_point() + theme_classic() + 
  labs(x = "% of bills kept",
       y = "% of votes kept",
       col = "Cutoff value") +
  geom_label(label = "Final\nCutoff",
            data = voting_behavior %>% filter(j == 130),
            nudge_x = -.15, nudge_y = .1,
            col = "black") +
  geom_curve(data = voting_behavior %>% filter(j == 130), curvature = 0,
            aes(x = n_votes -.095, xend = n_votes,
                y = n_votes_cast + .1, yend = n_votes_cast),
            col = "black", size = .9, 
            arrow = arrow(length = unit(0.03, "npc"))) +
  geom_point(col = "red", data = voting_behavior %>% filter(j == 130))
```

<img src="Single-cell_files/figure-markdown_github/filtering genes-1.png"  style="display: block; margin: auto;" />

We can see a clear inflexion point in the curve, that we pick as the final cutoff value. Our final filtering rule for bills (genes) is "At least 130 delegates cast at vote."

``` r
# Filter
Assay <- Assay[rowSums(exprs(Assay) != 0) >= 130, ]
```

Normalization and correcting for zero inflation
===============================================

In scRNA-seq settings, after filtering comes a normalization step. scRNA-sequencing has a lot of technical issues and this technical noise can be stronger than any biological signal. Normalization tries to adjust for this. Here, we have no technical noise so we do not need to adjust for it.

This is also the time to adjust for dropouts. Indeed, because of technical difficulties, some genes may seem non-expressed (i.e their count is zero) while they in fact are. Zero-inflation methods try to adjust for that by categorizing the zeros into real and dropouts and, in the case of dropouts, try to input the missing values.

To better understand dropouts, let us go back to our example. In our voting record matrix, zeros represent abstention (i.e real zeros) and not being able to votes (i.e dropouts, we do not know what the delegate would have done if he could have voted). Adjusting for dropouts means identifying which zeros are true zeros (true abstentions) and then, for dropouts, estimating the real value. Since the technique that would be used to correct for dropouts in our political example would be different to biostatistical techniques used in scRNA-seq settings, we will not cover it here.

Recovering / discovering biological types
=========================================

In most scRNA-seq settings, we have no information about the biological type. We therefore want to cluster the cells and then use marker genes to identify known-cell types. Alternatively, we can use clustering to discover new cell types and identify marker genes for those new types.

In our example, let us imagine for a moment that we do not know the political group to which each delegate belong. Can we recover those groups? This will be done in two steps: discover stable clusters of delegates (cells) then identify vote (genes) that help distinguish those clusters.

Clustering
----------

We use a broad clustering framework based on the *ClusterExperiment* package, that test many clustering algorithms and use all the results to find stable clusters. As a result, some cells (delegates) are left unclustered.

``` r
# Build the object
se <- SummarizedExperiment(assays = exprs(Assay))

# Run all clustering methods
ce <- clusterMany(x = se, k = seq(5, 50, 5),
                  clusterFunction = c("pam", "hierarchicalK"),
                  reduceMethod = c("PCA"), 
                  ks = 2:20, minSizes = 5,
                  run=TRUE, transFun = function(x){x})
# Find stable clusters
ce <- combineMany(ce, minSize = 12, proportion = 0.6,
                  clusterLabel = "combineMany")

# Merge clusters
ce <- makeDendrogram(ce)
ce <- mergeClusters(ce, mergeMethod="adjP", plotInfo = "adjP",
                  cutoff = 0.35, clusterLabel = "Clusters", plot = F)
```

After running the algorithm, we find 6 clusters that we can visualize on the first 2 principal components.

``` r
# Plot clusters
pca <- prcomp(t(exprs(Assay)))
Clusters <- data.frame(clusters = ce@clusterMatrix[,1],
                       x = pca$x[,1], y = pca$x[,2],
                       group = phenoData(Assay)@data$group1)
ggplot(Clusters, aes(x = x, y = y, col = factor(clusters))) +
  geom_point(size = 2, alpha = 0.7) +
  labs(col = "Clusters", x = "PC1", y = "PC2") +
  scale_color_manual(values = brewer.pal(7, "Set2"),
                    breaks = c(-1, 1:6),
                    labels = c("Unclustered", paste("Cluster", 1:6))) +
  theme_classic() + 
  theme(legend.title = element_text(hjust = 0.7))
```

<img src="Single-cell_files/figure-markdown_github/plot-1.png"  style="display: block; margin: auto;" />

We can see how those clusters compare with the actual political groups

``` r
ggplot(Clusters, aes(x = x, y = y,
       col = factor(group, levels = c("Independent", "LR", "SER", "GDR",
                                      "Greens", "RRPD", "UDI")))) +
  geom_point(size = 2, alpha = 0.7) +
  labs(col = "Clusters", x = "PC1", y = "PC2") +
  scale_color_manual(values = brewer.pal(7, "Set2")) +
  theme_classic() + 
  theme(legend.title = element_text(hjust = 0.7))
```

<img src="Single-cell_files/figure-markdown_github/plot groups-1.png"  style="display: block; margin: auto;" />

We can see that we correctly identified 5 main groups: SER (left, governing party), as well as the two main opposition parties, LR (right) and UDI (center-right). The far left party GDR is also correclty identify. However, we cannot identify the RRPD (center-left), which was part of the governing coalition, nor the greens (which can be identified if we choose a smaller minimum cluster size but we then start to see to many clusters).

We can morever see that we partitionshed the governing coalition in three. This situation has an analog in biological settings. We cluster our data in an unsupervised manner even though we know some rough biological types. Then we can compare our clusters and the types. If, like it does here, we have some good matches, we can be confident in our clustering techniques. Any new partition we discover might therefore be interesting and can be investigated.

Finding differentially expressed genes
--------------------------------------

We can now look at the votes and ask: which votes set a given group appart from the others? This can be done in two ways. Either looking at the clusters we found to try and identify them (like we will do for clusters 2, 4 & 5), or look at the political labels that we know and find the votes that separate them. In the first case, we care about the characterizing the clusters, in the other about finding bills that create a lot of political differentiation.

In scRNA-seq settings, the first case would be looking at clusters and use known-marker genes, gene ontology or pathway analysis to identify them. The other case would be identifying knwon-marker genes, or finding DE genes between treatment and control, or many other possible scenario.

In this example, we will look at the second case: which votes distinguish political groups? We can see it in two regards in our example. What votes distinguish each group from the governing coalitions (SER-RRPD-Greens)? What votes distinguish each group from all others?

``` r
groups <- phenoData(Assay)$group1
groups[groups %in%  c("SER", "Greens", "RRPD")] <- "Gov"
groups <- as.numeric(factor(groups), levels = c("Independent", "LR", "UDI",
                                                "GRD", "GOV"))

contrasts <- getBestFeatures(exprs(Assay), groups,
                             contrastType = "Pairs",
                             number = Inf) %>%
             filter(str_detect(Contrast, "-Cl05")) %>%
             group_by(Contrast) %>%
             arrange(P.Value) %>%
             top_n(-5, P.Value) %>%
             select(Contrast, Feature) %>%
             mutate(Feature = as.numeric(Feature)) %>%
             left_join(Assay@featureData@data,
                       by = c("Feature" = "Number")) %>%
  select(-Feature) %>%
  group_by(Contrast) %>%
  mutate(groupid = row_number()) %>%
  spread(Contrast, Title) %>%
  select(-groupid)

colnames(contrasts) <- c("Independent", "LR (right)", "UDI (center-right)",
                         "GRD (far-right)")
kable(contrasts)
```

| Independent                                                                                                                                                                                                                    | LR (right)                                                                                                     | UDI (center-right)                                                                                                                                               | GRD (far-right)                                                                                                                                                  |
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| la motion de censure déposée en application de l'article 49, alinéa 2 de la Constitution par M. Christian Jacob et 144 membres de l'Assemblée.                                                                                 | l'ensemble de la première partie du projet de loi de finances pour 2014.                                       | sur l'ensemble du projet de loi organique relatif à la programmation et à la gouvernance des finances publiques.                                                 | l'ensemble du projet de loi relatif à la lutte contre la fraude fiscale et la grande délinquance économique et financière.                                       |
| l'ensemble du projet de loi relatif à la prévention de la récidive et à l'individualisation des peines.                                                                                                                        | l'ensemble du projet de loi de financement de la sécurité sociale pour 2014.                                   | la proposition de loi organique relative à la nomination du président du conseil d'administration de l'Agence française pour la biodiversité (première lecture). | l'ensemble du projet de loi relatif à la lutte contre la fraude fiscale et la grande délinquance économique et financière.                                       |
| l'ensemble du projet de loi relatif à la mobilisation du foncier public en faveur du logement et au renforcement des obligations de production de logement social (première lecture).                                          | l'ensemble du projet de loi d'orientation et de programmation pour la refondation de l'école de la République. | la proposition de loi organique relative à la nomination des dirigeants de la SNCF.                                                                              | la proposition de loi organique relative à la nomination du président du conseil d'administration de l'Agence française pour la biodiversité (première lecture). |
| l'ensemble de la proposition de loi tendant à modifier la loi n° 2011-814 du 7 juillet 2011 relative à la bioéthique en autorisant sous certaines conditions la recherche sur l'embryon et les cellules souches embryonnaires. | la première partie du projet de loi de finances pour 2013.                                                     | l'ensemble du projet de loi pour une République numérique (première lecture).                                                                                    | l'ensemble du projet de loi relatif à la réforme de l'asile (première lecture).                                                                                  |
| l'ensemble de la proposition de loi constitutionnelle visant à rendre constitutionnel le principe d'indisponibilité du corps humain (première lecture).                                                                        | l'ensemble du projet de loi de financement de la sécurité sociale pour 2014 (nouvelle lecture).                | l'ensemble du projet de loi relatif à la réforme de l'asile (première lecture).                                                                                  | le projet de loi organique relatif au procureur de la République financier (lecture définitive).                                                                 |

R session
=========

``` r
sessionInfo()
```

    ## R version 3.5.0 (2018-04-23)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.5
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] bindrcpp_0.2.2              stringr_1.3.1              
    ##  [3] RColorBrewer_1.1-2          clusterExperiment_2.0.2    
    ##  [5] SingleCellExperiment_1.2.0  SummarizedExperiment_1.10.1
    ##  [7] DelayedArray_0.6.1          BiocParallel_1.14.2        
    ##  [9] GenomicRanges_1.32.3        GenomeInfoDb_1.16.0        
    ## [11] IRanges_2.14.10             S4Vectors_0.18.3           
    ## [13] readr_1.1.1                 ggplot2_3.0.0              
    ## [15] tidyr_0.8.1                 dplyr_0.7.6                
    ## [17] matrixStats_0.53.1          Biobase_2.40.0             
    ## [19] BiocGenerics_0.26.0         knitr_1.20                 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-137           bitops_1.0-6           progress_1.2.0        
    ##  [4] httr_1.3.1             doParallel_1.0.11      rprojroot_1.3-2       
    ##  [7] prabclus_2.2-6         tools_3.5.0            backports_1.1.2       
    ## [10] R6_2.2.2               HDF5Array_1.8.1        lazyeval_0.2.1        
    ## [13] colorspace_1.3-2       ade4_1.7-11            trimcluster_0.1-2     
    ## [16] nnet_7.3-12            withr_2.1.2            prettyunits_1.0.2     
    ## [19] tidyselect_0.2.4       gridExtra_2.3          compiler_3.5.0        
    ## [22] xml2_1.2.0             pkgmaker_0.27          labeling_0.3          
    ## [25] diptest_0.75-7         scales_0.5.0           DEoptimR_1.0-8        
    ## [28] mvtnorm_1.0-8          robustbase_0.93-1      NMF_0.21.0            
    ## [31] digest_0.6.15          rmarkdown_1.10         XVector_0.20.0        
    ## [34] pkgconfig_2.0.1        htmltools_0.3.6        bibtex_0.4.2          
    ## [37] highr_0.7              limma_3.36.2           rlang_0.2.1           
    ## [40] howmany_0.3-1          bindr_0.1.1            mclust_5.4.1          
    ## [43] dendextend_1.8.0       RCurl_1.95-4.10        magrittr_1.5          
    ## [46] modeltools_0.2-21      GenomeInfoDbData_1.1.0 Matrix_1.2-14         
    ## [49] Rhdf5lib_1.2.1         Rcpp_0.12.17           munsell_0.5.0         
    ## [52] ape_5.1                viridis_0.5.1          stringi_1.2.3         
    ## [55] whisker_0.3-2          yaml_2.1.19            MASS_7.3-50           
    ## [58] zlibbioc_1.26.0        rhdf5_2.24.0           flexmix_2.3-14        
    ## [61] plyr_1.8.4             grid_3.5.0             crayon_1.3.4          
    ## [64] rncl_0.8.2             lattice_0.20-35        splines_3.5.0         
    ## [67] hms_0.4.2              pillar_1.2.3           uuid_0.1-2            
    ## [70] fpc_2.1-11             rngtools_1.3.1         reshape2_1.4.3        
    ## [73] codetools_0.2-15       XML_3.98-1.11          glue_1.2.0            
    ## [76] evaluate_0.10.1        RNeXML_2.1.1           data.table_1.11.4     
    ## [79] foreach_1.4.4          locfdr_1.1-8           gtable_0.2.0          
    ## [82] purrr_0.2.5            kernlab_0.9-26         assertthat_0.2.0      
    ## [85] gridBase_0.4-7         phylobase_0.8.4        xtable_1.8-2          
    ## [88] RSpectra_0.13-1        class_7.3-14           viridisLite_0.3.0     
    ## [91] tibble_1.4.2           iterators_1.0.9        registry_0.5          
    ## [94] cluster_2.0.7-1

Thanks
======

Special thanks to Vincent Viers for the initial inspiration of this project.

The code used for scraping the data is based on the blog post <https://freakonometrics.hypotheses.org/50973>.
