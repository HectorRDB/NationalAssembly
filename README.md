A short tutorial to help beginners in biostatistics and genomics understand an usual single-cell RNA-sequencing workflow
================
Hector Roux de Bézieux
May 2018

Aim
============

Single-cell RNA-sequencing (scRNA-seq) is a very potent biological tool used for many applications. A far from exhaustive list would include identifying new cell types, finding differentially expressed (DE) genes, and discovering lineages among cells. To get an overview of the principle, see [here](https://en.wikipedia.org/wiki/Single_cell_sequencing#Single-cell_RNA_sequencing_(scRNA-seq)).

However, the usual framework might seem a little daunting for beginners and, while many well-crafted tutorials exists, they all share the same idea: use a biological dataset as an example. Here, we want to use a dataset that would be more understandable to a broader public to explain the usual steps in scRNA-seq analysis

Description
============

The data consists of the record of votes of all delegates of the <a href="https://www.codecogs.com/eqnedit.php?latex=14^{th}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?14^{th}" title="14^{th}" /></a> legislature  (2012-2017) of the French parliament. Only some basics on data analysis are needed to understand this tutorial. Knowledge of R helps to understand the code but is not necessary to follow along.

To access the tutorial, see [here](https://hectorrdb.github.io/NationalAssembly/).

Thanks
======

Special thanks to [Vincent Viers](https://github.com/vviers) for the initial inspiration of this project.  The code used for scraping the data is based on the blog post https://freakonometrics.hypotheses.org/50973.
