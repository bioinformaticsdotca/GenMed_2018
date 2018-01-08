---
layout: workshop_main_2day
permalink: /genomic_medicine_2018
title: Genomic Medicine
header1: Workshop Pages for Students
header2: Bioinformatics of Genomic Medicine
image: /site_images/CBW_population_icon.jpg
---

# Welcome <a id="welcome"></a>

Welcome to Bioinformatics of Genomic Medicine.  

## Install these tools before the workshop:  

1) A robust text editor.   

* For Windows/PC - [notepad++](http://notepad-plus-plus.org/)  
* For Linux - [gEdit](http://projects.gnome.org/gedit/)  
* For Mac – [TextWrangler](http://www.barebones.com/products/textwrangler/download.html)

2) A file decompression tool.  

* For Windows/PC – [7zip](http://www.7-zip.org/).  
* For Linux – [gzip](http://www.gzip.org).   
* For Mac – already there.

3) A robust internet browser such as Firefox or Safari (Internet Explorer and Chrome are not recommended because of Java issues).

4) SSH client - Mac and Linux users already have a command line ssh program that can be run from the terminal. For Windows users, please download [PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html).  

5) SCP/SFTP client - We will be moving data from the servers to the student laptops for visualization. Mac and Linux users already have a command line scp and sftp program. For Windows users, please install [WinSCP](http://winscp.net/eng/download.php).

6) A PDF viewer (Adobe Acrobat or equivalent).

7) Latest version of [R](http://www.r-project.org/) and [RStudio](http://www.rstudio.com/).  If asked to also install git, select yes.

8) Install the BioConductor core packages. To do this, open RStudio and type the following at the > prompt:

```r
source("http://bioconductor.org/biocLite.R");
biocLite();
```

9) Install some workshop specific packages in RStudio. At the > prompt, type:

```r
install.packages("SNFtool");
install.packages("gplots");
install.packages("RColorBrewer");
install.packages("glmnet");
install.packages("ggplot2");
install.packages("reshape2");
install.packages("pROC");
install.packages("caret");
install.packages("matrixStats");
source("http://bioconductor.org/biocLite.R");
biocLite("minfi");
biocLite("bumphunter");
biocLite("limma");
biocLite("lumi");
biocLite("sva");

```  

## Read these papers before the workshop:  

[The Human Phenotype Ontology: A Tool for Annotating and Analyzing Human Hereditary Disease](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2668030/)  
[PhenoTips: patient phenotyping software for clinical and research use](https://www.ncbi.nlm.nih.gov/pubmed/23636887)  
[IHEC Data portal](http://www.cell.com/cell-systems/abstract/S2405-4712(16)30362-3)  

## Do these tutorials before the workshop:  

1) **R Preparation tutorials**: You are expected to have completed the following tutorials in **R** beforehand. The tutorial should be very accessible even if you have never used **R** before.

* The [CBW R tutorial](http://bioinformatics-ca.github.io/CBW_R_Tutorial/)
* [Quick and Dirty Guide to **R**](http://ww2.coastal.edu/kingw/statistics/R-tutorials/text/quick&dirty_R.txt)  
* The [R Tutorial](http://www.cyclismo.org/tutorial/R/) up to and including 5. Basic Plots
* The [R command cheat sheet](https://github.com/bioinformaticsdotca/bioinformaticsdotca.github.io/blob/master/resources/R_Short-refcard.pdf)

2) **UNIX Preparation tutorials**:  

* [UNIX Bootcamp](http://rik.smith-unna.com/command_line_bootcamp/?id=9xnbkx6eaof)
* [Unix Cheat sheet](http://www.rain.org/~mkummel/unix.html) 

3) [Sequencing Terminology](http://www.ncbi.nlm.nih.gov/projects/genome/glossary.shtml)

Please note that these instructions may change prior to the workshop.  
