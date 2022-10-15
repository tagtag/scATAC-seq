 # scATAC-seq

This is  a sample R code for
"Tensor Decomposition Discriminates Tissues Using scATAC-seq"
by Y-H. Taguchi  and Turki Turki

https://doi.org/10.1101/2022.08.04.502875

============================

For the people who would like to try this procedure by themselves, we put a R script, TD_scATAC_seq.R, which includes a function TD_scATAC_seq that can perform the analysis in the above paper.

Usage:

TD_scATAC_seq(x_files,L,nseed)

Input:

x_files:  A character vector whose length is equal to the number of data sets analyzed. Each component is the name of file (full path) to be loaded to be analyzed. Each file should be a R object that includes a sparse matrix named "x_all" whose number of columns and rows are equal to the number of single cells and binned genomic regions, respectively. The values stored in a sparse matrix is binned accessibility (for details, see the paper)

L: the number of component generated and used as input for UMAP. It is denoted as V_{\ell' j} in the paper. Default value is 10.

nseed: random seed for UMAP. Default value is 0.

Output: A UMAP object  whose columns and rows correspond  to single cells and embedding dimensions, respectively.
Suppose the number of single cells in individual x_alls  are N1, N2, .....  The first N1 columns correspond to single cells in the first x_all, the columns from N1+1 to N1+N2 correspond to those in the second x_all, and so on.

Example:
Suppose eight files, x_all, x_all_1, ....,  x_all_7, are placed in the current directory. Then x_files should be

x_files <- c("x_all","x_all_1","x_all_2","x_all_3","x_all_4","x_all_5","x_all_6","x_all_7")

Then execute

result <- TD_scATAC_seq(x_files)

This returns an UMAP object whose columns and rows correspond  to single cells and the ten embedding  dimensions, respectively.


PS We are not responsible to generate binned accessibility from the files available. In the case in our study, check sample.R in this repository.


