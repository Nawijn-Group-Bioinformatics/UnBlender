# UnBlender workflow

## Step 1:
 - The user chooses which sample type they wish to deconvolute, and which cell types are relevant to their research question.
 - After selection of the cell types of interest, reference gene expression profiles of these cell types are randomly sampled
from the HLCA, based on samples of the same type. A total of 300 cells x the number of cell types is sampled.
 - UnBlender also generates HLCA-derived pseudobulk samples of the selected sample type. Pseudobulk samples are created by summing
the gene counts across the cells in a scRNA-seq sample, thereby making the single cell sample resemble a bulk RNA-seq sample.
 - These ‘pseudobulks’ are then deconvoluted using using the reference data that was sampled for the chosen cell types of interest. 
Because we know the true cell type composition of the pseudobulk samples, we can evaluate whether this deconvolution strategy is feasible,
or in other words: whether we can trust the results (per cell type).

## Step 2: 
 - The user decides whether the deconvolution accuracy is high enough. 
 - If so, the same reference data is used to deconvolute the real bulk RNA-seq dataset of interest. 
 - If not, the user can go back to step 1 and adjust their deconvolution strategy: changing the cell type selection by e.g. choosing cell type labels of lower resolution (e.g. ‘macrophages’ instead of ‘alveolar macrophages’) or by merging categories that are too alike in terms of their gene expression profiles (e.g. merging ‘club cell’ and ‘goblet cell’ labels into ‘secretory’).

## Step 3:
- This optimisation may take a few tries. Once you’ve found a strategy that yields accurate results, you proceed to step 3 of the workflow, 
in which you can upload your own bulk RNA-seq dataset to deconvolute using the validated strategy. Alternatively, you can also download the reference data sampled in step 1, to perform cell type deconvolution on your own computer.



