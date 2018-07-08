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
-   [Thanks](#thanks)

Introduction
============

Single-cell RNA-seq (scRNA-seq) is a very potent biological tool used for many applications. A far from exhaustive list would include identifying new cell types, finding differentially expressed (DE) genes, and discovering lineages among cells. However, the usual framework might seem a little daunting for beginners and, while many well-crafted tutorials exists, they all share the same idea: use a biological dataset as an example. Here, we want to use a dataset that would be more understandable to a broader public to explain the usual steps in scRNA-seq analysis.

Context
=======

Our dataset is the voting record of the French delegates of the 14<sup>*t**h*</sup> legislature, going from May 2012 to May 2017. Delegates are characterized by their names and surnames, their department and circonscription, their political group. For each of the 644 votes (which can be bills, amendments, choosing the prime minister, ...), we also have the voting behavior of each delegate (voted yes, no, abstain or did not took part in the vote) and the reason of the vote. We are interested in investigating the link between voting behavior and political affiliations.

Loading the data and building an ExpressionSet object
=====================================================

The usual scRNA-seq data is a matrix of *J* genes by *n* cells. Each cell represent the number of copies of the mRNA of a given gene in a given cell. Here, we instead have a *J* bills by *n* delegates matrix. Each cell represent the vote of the delegate for that low (−1 if against, 1 is for, 0 if not voting or abstained). We store the data in an *ExpressionSet* object, where the *phenodata* is the information about the delegates (the cells) and the *featuredata* is the information about the votes (the genes).

``` r
voting_record <- read_csv("data/voting_record.csv")
meta <- voting_record %>% select(circo, dept, name, surname, identifiant,
                                 iden, group, NbVote, NbYes, NbNo, NbAbst)
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

In scRNA-seq settings, we want to filter the cells with too few total reads, or reads of too poor quality (that do not match to the genome of reference). Here, we instead filter the delegates with too few votes. This may happen either because a delegate had to leave office before the end of the session, or was appointed to the government, and then its replacement is not around long enough.

We therefore filter delegates whose number of votes is more than 2 SD below the mean (so equal or less to 40).

``` r
Nb_Votes <- data.frame(Nb_Votes = colSums(!is.na(exprs(Assay))),
                       delegate = rownames(phenoData(Assay)))

ggplot(Nb_Votes, aes(x = Nb_Votes, fill = Nb_Votes <= 40)) +
                 stat_bin(breaks = seq(0, max(Nb_Votes$Nb_Votes), 20)) +
                 labs(x = "Number of votes") +
            guides(fill = guide_legend(title = "Could vote\nless than\n40 times"))
```

<img src="Single-cell_files/figure-markdown_github/filtering cells-1.png"  style="display: block; margin: auto;" />

``` r
Assay <- Assay[, Nb_Votes$Nb_Votes > 40]
```

We then simplify the data by replacing NAs by 0 (as abstentions). We can then check whether some delegates voted only rarely. Again, we remove people whose number of votes is more than 2 SD below the mean number of votes (so delegates who voted less than 38 votes). Among those, we of course find the president of the parliament, who does not take part in the votes.

``` r
E <- exprs(Assay)
E[is.na(E)] <- 0
pd <- phenoData(Assay)
fd <- featureData(Assay)
colnames(E) <- rownames(pd) <- pd@data$identifiant
rownames(E) <- rownames(fd) <- fd@data$Number
Assay <- ExpressionSet(assayData = E, phenoData = pd, featureData = fd)

Nb_Votes <- data.frame(Nb_Votes = colSums(exprs(Assay) != 0),
                       delegate = rownames(phenoData(Assay)))
ggplot(Nb_Votes, aes(x = Nb_Votes, fill = Nb_Votes <= 38)) +
                 stat_bin(breaks = seq(0, max(Nb_Votes$Nb_Votes), 19)) +
                 labs(x = "Number of votes") +
            guides(fill = guide_legend(title = "Actually voted\nless than\n38 times"))
```

<img src="Single-cell_files/figure-markdown_github/change NA to zeros-1.png"  style="display: block; margin: auto;" />

``` r
Assay <- Assay[, Nb_Votes$Nb_Votes > 38]
```

Filtering the genes (votes)
---------------------------

In scRNA-seq settings, we want to filter out genes that are expressed in a small number of cells. This is done for two reasons. The first is computational. Most datasets would have dozens of thousands of genes, over dozens or hundreds of cells. Reducing the number of genes to keep speeds up calculations. Secondly, we are usually interested in a specific tissue or process. Most genes are not being expressed in the samples of interest and just add noise to the analysis. Hence filtration. Simple heuristics in the form of "at least *i* counts in *j* cells" are usually quite efficient. The aim is too be quite broad. We want to filter as many genes as possible while retaining most of the reads. In usual settings, you would only keep ~20% of the genes at most but that would still capture more than 90% of the reads.

In our context, we only have 644 votes so the computational goal does not really hold. But we can filter out votes where many delegates did not vote. Those are probably non-political technical votes. We therefore use the filter rules "at least X delegates cast a vote on the bill".

``` r
Nb_Votes <- rowSums(exprs(Assay) != 0)
voting_behavior<- data.frame(n_votes_cast = rep(0, nrow(Assay)),
                             n_votes = rep(0, nrow(Assay)),
                             X = 1:nrow(Assay))
for(X in 1:nrow(Assay)){
  Nb_Votes <- Nb_Votes[Nb_Votes >= X]
  voting_behavior$n_votes[X] <- length(Nb_Votes)
  voting_behavior$n_votes_cast[X] <- sum(Nb_Votes)
}

ggplot(voting_behavior,
       aes(x = n_votes, y = n_votes_cast, col = X)) +
  geom_point() +
  labs(x = "Number of bills kept",
       y = "Number of votes kept",
       col = "Cutoff value") +
  geom_point(col = "red", data = voting_behavior %>% filter(X == 130)) +
  geom_label(label = "Final\nCutoff",
            data = voting_behavior %>% filter(X == 130),
            nudge_x = -120, nudge_y = 1E4,
            col = "black") +
  geom_curve(data = voting_behavior %>% filter(X == 130), curvature = 0,
            aes(x = n_votes -70, xend = n_votes,
                y = n_votes_cast + 1E4, yend = n_votes_cast),
            col = "black",
            arrow = arrow(length = unit(0.03, "npc")))
```

<img src="Single-cell_files/figure-markdown_github/filtering genes-1.png"  style="display: block; margin: auto;" />

We can see a clear inflexion point in the curve, that we pick as the final cutoff value.

Normalization and correcting for zero inflation
===============================================

Recovering / discovering biological types
=========================================

In most scRNA-seq settings, we have no information about the biological type. We therefore want to cluster the cells and then use marker genes to identify known-cell types. Alternatively, we can use clustering to discover new cell types and identify marker genes for those new types.

In our example, let us imagine for a moment that we do not know the political group to which each delegate belong. Can we recover those groups? This will be done in two steps: discover stable clusters of delegates (cells) then identify vote (genes) that help distinguish those clusters. \#\# Clustering

``` r
pca <- prcomp(t(exprs(Assay)), center = T, scale = T)
ggplot(data.frame(x = pca$x[,1], y = pca$x[,2], col = phenoData(Assay)@data$group),
       aes(x = x, fill = col)) +
  facet_wrap(col~.) + 
  geom_density()
```

<img src="Single-cell_files/figure-markdown_github/clustering-1.png"  style="display: block; margin: auto;" />

Thanks
======

Special thanks to Vincent Viers for the initial inspiration of this project.

The code used for scraping the data is based on the blog post <https://freakonometrics.hypotheses.org/50973>.
