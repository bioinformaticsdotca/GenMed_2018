---
layout: tutorial_page
permalink: /GenMed_2017_module6_lab
title: GenMed Lab 6
header1: Workshop Pages for Students
header2: Genomic Medicine 2017 Module 6 Lab
image: /site_images/CBW_population_icon.jpg
home: https://bioinformaticsdotca.github.io/genomic_medicine_2017
description: Epigenetic Profiling in Disease
author: Andrei Turinsky
modified: May 11th, 2017
---

# Module 6: Epigenetic Profiling in Disease 

by Andrei Turinsky, *PhD*

## Introduction

### Description of the lab
In this module's lab, we will use R packages to analyse epigenetic data related to disease.

* First, we will examine a Down syndrome dataset at the Gene Ontology Omnibus (GEO).
* Second, we will explore the R/Bioconductor minfi package.
* Third, we will retrieve the GEO dataset using and convert in into minfi data structure.
* Finally, we will use a custom-made R script to perform a series of exploratory data analysis tasks, such as clustering, classification, batch correction, and detection of differentially methylated sites and regions. 

### Local software that we will use
* A web browser
* RStudio


## Tutorial

### Download a Gene Expression Omnibus microarray annotation file  
* Create a new directory (folder) for Module 6.

* As preparation for future analysis steps, download a large file [GPL13534.txt](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL13534&form=text&view=full) from GEO (size 215 Mb) into the newly created directory.


### Explore the Gene Expression Omnibus dataset
* Open a web browser on your computer, and load the URL [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52588). It will present a repository with 87 DNA methylation samples obtained using Illumina HumanMethylation450 microarrays: 29 methylomes corresponding to Down sysndrome individuals as well as matching methylomes of unaffected siblings and mothers. 

* Examine the data. Note the Series Matrix (recommended) and SOFT data files available for download... but no need to download them now, the data will be retrieved by the R script later during the tutorial. 

* Follow the link [Analyze with GEO2R](https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE52588) and explore the features of the automated analysis pipeline provided by the GEO - including the automatically generated R script.

### Explore minfi tutorial  
* Open the web page for the [minfi tutorial](https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html) and examine the Contents. This will give you an idea what data analysis tasks can be accomplished using minfi.

### Data retrieval and exploratory analysis in R

* Place the provided [R script](https://github.com/bioinformaticsdotca/Genomic_Med_2017/blob/master/mod6/cbw_mod6.R) into the newly created directory for Module 6.

* Open the RStudio. Follow the menu File > New Project > Existing Directory, and choose your Module 6 directory to create a new R project.  
   
* Load the R script into the RStudio.

* For the remainder of the tutorial, run the R script line by line, e.g. by pressing Ctrl-Enter or (Cmd-Enter on Mac). Examine the results and the script comments... until you are done! 
 
