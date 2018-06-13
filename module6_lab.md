---
layout: tutorial_page
permalink: /GenMed_2018_module6_lab
title: GenMed Lab 6
header1: Workshop Pages for Students
header2: Genomic Medicine 2018 Module 6 Lab
image: /site_images/CBW_population_icon.jpg
home: https://bioinformaticsdotca.github.io/genmed_2018
description: Epigenetic Profiling in Disease
author: Andrei Turinsky
modified: May 11th, 2017
---

# Module 6: Epigenetic Profiling in Disease 

by Andrei Turinsky, *PhD*

## Introduction

### Description of the lab
In this module's lab, we will use R packages to analyse epigenetic data related to disease (Down syndrome dataset extracted from Gene Expression Omnibus). We will use a custom-made R script to perform a series of exploratory data analysis tasks, such as clustering, classification, detection of differentially methylated sites and regions, and batch correction. 

### Local software that we will use
* A web browser
* RStudio


## Tutorial

### Exploratory analysis of DNA methylation dataset in R

* Create a new directory (folder) for Module 6.

* Place the provided [R script](https://github.com/bioinformaticsdotca/GenMed_2018/tree/master/mod6/cbw-mod6-2018.R) into the newly created directory for Module 6.

* Place 3 data files [samplesheet.csv](https://github.com/bioinformaticsdotca/GenMed_2018/tree/master/mod6/samplesheet.csv), [betas.csv](https://github.com/bioinformaticsdotca/GenMed_2018/tree/master/mod6/betas.csv) and [annotations.csv](https://github.com/bioinformaticsdotca/GenMed_2018/tree/master/mod6/annotations.csv) into the newly created directory for Module 6.

* Open the RStudio. Follow the menu File > New Project > Existing Directory, and choose your Module 6 directory to create a new R project.  
   
* Load the R script into the RStudio.

* For the remainder of the tutorial, run the R script line by line, e.g. by pressing Ctrl-Enter or (Cmd-Enter on Mac), following the instructor. Examine the results and the script comments... until you are done! 
 
