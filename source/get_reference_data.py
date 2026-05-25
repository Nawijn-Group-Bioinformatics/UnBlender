### About this script ###
# Extract the reference data for the signature matrix:
# X (default 300) cells per cell type from the relevant sample type. 
# (Default: non-balanced, i.e. randomly sampled; balanced sampling 
# results in equal numbers of reference cells per cell type.)

# command line arguments: tissue n_cells path/to/HLCA.h5ad output_path/ config.yml balanced_ref

print('### Generating reference data')

### Load libraries
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
import sys
import warnings

# Set up & load HLCA data
sample_type = sys.argv[1]
n_cells = int(sys.argv[2])
output_path = sys.argv[4]
config_yml = sys.argv[5]
balanced_ref = sys.argv[6]
print("loading data")
adata = sc.read(sys.argv[3])

# Suppress FutureWarnings
warnings.filterwarnings("ignore", category=FutureWarning)


### Cell type selection (from config file)
print('loading cell type selection')
with open(config_yml, "r") as stream:

	try:
		temp = yaml.safe_load(stream)
		include_cell_types = temp['cell_types']

	except yaml.YAMLError as exc:
		print(exc)

# run check for typing errors
warning = False
for item in include_cell_types:
    
    if type(item[1]) == str: # if just the one cell type (simple)
        if item[1] not in adata.obs[item[0]].unique():
            print("Error: " + item[1] + " not found in " + item[0])
            warning = True
        
    else: # if list of cell types to be combined
        for cell_type in item[1]:
            if cell_type not in adata.obs[item[0]].unique():
                print("Error: " + cell_type + " not found in " + item[0])
                warning = True
                
if warning == False:
    print("... selected cell types all found!")
else:
	print("... refer to cell_type_names.tsv for available cell type names and levels")
	sys.exit("... shutting down.")


### Subset HLCA to relevant sample types
print('selecting sample type')
if sample_type == "parenchyma":
	adata = adata[adata.obs['anatomical_region_coarse'] == 'parenchyma']

elif sample_type == "bronchial_biopsy":
	adata = adata[adata.obs['sample_type'] == 'biopsy']
	adata = adata[adata.obs['anatomical_region_coarse'].isin(['airway','Intermediate Bronchi',"Trachea"])]

elif sample_type == "nasal_brush":
	adata = adata[adata.obs['sample_type'].isin(['brush', 'scraping'])]
	adata = adata[adata.obs['anatomical_region_coarse'].isin(['nose','Inferior turbinate'])]

elif sample_type == "bronchial_brush":
	adata = adata[adata.obs['sample_type'] == 'brush']
	adata = adata[adata.obs['anatomical_region_coarse'] == "Distal Bronchi"]

else:
	sys.exit("ERROR: sample_type not valid: " + sample_type)


### Subsampling the reference data
print('downsampling reference data')
if balanced_ref == 'balanced_true':

    # Subsample to n_cells cells per cell type
    count_subsets = 0
    for item in include_cell_types:
        
        # subset adata to relevant cell types
        if type(item[1]) == str:
            adata_subset = adata[adata.obs[item[0]] == item[1]].copy()
            adata_subset.obs['custom_label'] = item[1]
        else:
            adata_subset = adata[adata.obs[item[0]].isin(item[1])].copy()
            adata_subset.obs['custom_label'] = ' & '.join(item[1])
        
        # subsample to n_cells cells
        if adata_subset.obs.shape[0] > n_cells:
            sc.pp.subsample(adata_subset, n_obs=n_cells, random_state=None) # random seed: set to 0 for tests
        
        # merge with previous cell types' data
        if count_subsets == 0:
            adata_subsampled = adata_subset
        else:
            adata_subsampled = adata_subsampled.concatenate(adata_subset)
        
        count_subsets +=1

elif balanced_ref == 'balanced_false':

    # Subset to cell types in include_cell_types
    count_subsets = 0
    for item in include_cell_types:
        
        # subset adata to relevant cell types
        if type(item[1]) == str:
            adata_subset = adata[adata.obs[item[0]] == item[1]].copy()
            adata_subset.obs['custom_label'] = item[1]
        else:
            adata_subset = adata[adata.obs[item[0]].isin(item[1])].copy()
            adata_subset.obs['custom_label'] = ' & '.join(item[1])
                
        # merge with previous cell types' data
        if count_subsets == 0:
            adata_subsampled = adata_subset
        else:
            adata_subsampled = adata_subsampled.concatenate(adata_subset)
        
        count_subsets +=1

    # Then randomly subsample to n_cells * number_of_cell_types cells
    total_number_of_cells = n_cells * count_subsets
    sc.pp.subsample(adata_subsampled, n_obs=total_number_of_cells, random_state = None) # random seed: set to 0 for tests

else:
    sys.exit("ERROR: balanced_ref not valid: " + balanced_ref)


### Extract & write the reference data with HUGO or ENSG IDs # TODO combine the two code blocks
print('saving reference data')

## Write for HUGO
adata = adata_subsampled

# write scRNA-seq counts matrix with cols = samples and rows = genes
counts_layer = pd.DataFrame(adata.X.todense(), index=adata.obs.index, columns=adata.var.index)    
counts_layer = counts_layer.transpose()

# get labels as first row (calling it GeneSymbols to get the colname of the first col right)
# (not entirely sure if that's neccessary, but let's run with it..)
adata.obs['GeneSymbols'] = adata.obs['custom_label']
counts_withInfo = pd.DataFrame(adata.obs['GeneSymbols']).transpose()
counts_withInfo = counts_withInfo.append(counts_layer)

# write
counts_withInfo.to_csv(output_path + "/" + sample_type + "_subsampled_matrix_max" + str(n_cells) + "cells_HUGO.txt", 
                       header=False, sep='\t')


## Write for ENSG
adata = adata_subsampled
        
# write scRNA-seq counts matrix with cols = samples and rows = genes
counts_layer = pd.DataFrame(adata.X.todense(), index=adata.obs.index, columns=adata.var['gene_ids'])    
counts_layer = counts_layer.transpose()

# get labels as first row (calling it GeneSymbols to get the colname of the first col right)
# (not entirely sure if that's neccessary, but let's run with it..)
adata.obs['GeneSymbols'] = adata.obs['custom_label']
counts_withInfo = pd.DataFrame(adata.obs['GeneSymbols']).transpose()
counts_withInfo = counts_withInfo.append(counts_layer)

# write
counts_withInfo.to_csv(output_path + "/" + sample_type + "_subsampled_matrix_max" + str(n_cells) + "cells_ENSG.txt", 
                       header=False, sep='\t')
print('\n')
