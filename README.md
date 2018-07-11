Understanding RNA-seq analysis through an example: voting behavior at the French parliament
================
Hector Roux de Bézieux
May 2018

Introduction
============

Single-cell RNA-seq (scRNA-seq) is a very potent biological tool used for many applications. A far from exhaustive list would include identifying new cell types, finding differentially expressed (DE) genes, and discovering lineages among cells. However, the usual framework might seem a little daunting for beginners and, while many well-crafted tutorials exists, they all share the same idea: use a biological dataset as an example. Here, we want to use a dataset that would be more understandable to a broader public to explain the usual steps in scRNA-seq analysis. Only some basics on data analysis are needed to understand this tutorial. Knowledge of R helps to understand the code but is not necessary to follow along.

Thanks
======

Special thanks to Vincent Viers for the initial inspiration of this project.

The code used for scraping the data is based on the blog post <https://freakonometrics.hypotheses.org/50973>.
