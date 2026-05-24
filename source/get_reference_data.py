# Get reference data for the signature matrix:
# X cells per cell type included.

# command line arguments: tissue n_cells path/to/HLCA.h5ad output_path/ config.yml

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
import sys

sc.logging.print_header()

# set up & load HLCA data
sample_type = sys.argv[1] # tissue = "parenchyma"
n_cells = int(sys.argv[2]) #200
adata = sc.read(sys.argv[3]) # sc.read("../HLCA_v1.1_20220725.h5ad") 
output_path = sys.argv[4] #"../pseudobulk_analysis/"
config_yml = sys.argv[5]
balanced_ref = sys.argv[6]

print("CONFIG FILE")
print(config_yml)


### CELL TYPE SELECTION directly from config file
with open(config_yml, "r") as stream:

	try:
		temp = yaml.safe_load(stream)
		print(temp)
		include_cell_types = temp['cell_types']

	except yaml.YAMLError as exc:
		print(exc)


# run check for typing errors
warning = False
for item in include_cell_types:
    
    if type(item[1]) == str: # if just the one cell type (simple)
        if item[1] not in adata.obs[item[0]].unique():
            print("Error: " + item[1] + " not found in column " + item[0])
            warning = True
        
    else: # if list of cell types to be combined
        for cell_type in item[1]:
            if cell_type not in adata.obs[item[0]].unique():
                print("Error: " + cell_type + " not found in column " + item[0])
                warning = True
                
if warning == False:
    print("Chosen cell types all found!")
else:
	sys.exit("Shutting down.")


# Subset HLCA to relevant sample types
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

elif sample_type == "IPF_DEG_analysis": # these are also parenchymal samples
    adata = adata[adata.obs['anatomical_region_coarse'] == 'parenchyma']

else:
	sys.exit("ERROR: sample_type not valid: " + sample_type)

print('starting subsampling')

# Subsampling
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
            sc.pp.subsample(adata_subset, n_obs=n_cells, random_state=None) # RANDOM SEED, FOR ACTUAL TESTS MAKE SURE ITS 0.
        
        # merge with previous cell types' data
        if count_subsets == 0:
            adata_subsampled = adata_subset
        else:
            adata_subsampled = adata_subsampled.concatenate(adata_subset)
        
        count_subsets +=1

elif balanced_ref == 'balanced_false':

    # Subset to cell types in include_cell_types
    count_subsets = 0

    print(include_cell_types)

    for item in include_cell_types:

        print(item)
        
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
    sc.pp.subsample(adata_subsampled, n_obs=total_number_of_cells, random_state = None) # RANDOM SEED, FOR ACTUAL TESTS MAKE SURE ITS 0.


else:
    sys.exit("ERROR: balanced_ref not valid: " + balanced_ref)

print('subsampling done')

# Extract & write the reference data with HUGO or ENSG IDs # TODO combine the two code blocks

# Write for HUGO
adata = adata_subsampled

    
# write scRNA-seq counts matrix with cols = samples and rows = genes
counts_layer = pd.DataFrame(adata.X.todense(), index=adata.obs.index, 
                            columns=adata.var.index)    

counts_layer = counts_layer.transpose()


# get labels as first row (calling it GeneSymbols to get the colname of the first col right)
# (not entirely sure if that's neccessary, but let's run with it..)
adata.obs['GeneSymbols'] = adata.obs['custom_label']
counts_withInfo = pd.DataFrame(adata.obs['GeneSymbols']).transpose()

counts_withInfo = counts_withInfo.append(counts_layer)

# write
counts_withInfo.to_csv(output_path + "/" + sample_type + "_subsampled_matrix_max" + str(n_cells) + "cells_HUGO.txt", 
                       header=False, sep='\t')


# Write for ENSG
adata = adata_subsampled
        
# write scRNA-seq counts matrix with cols = samples and rows = genes
counts_layer = pd.DataFrame(adata.X.todense(), index=adata.obs.index, 
                            columns=adata.var['gene_ids'])    
counts_layer = counts_layer.transpose()

# get labels as first row (calling it GeneSymbols to get the colname of the first col right)
# (not entirely sure if that's neccessary, but let's run with it..)
adata.obs['GeneSymbols'] = adata.obs['custom_label']
counts_withInfo = pd.DataFrame(adata.obs['GeneSymbols']).transpose()
counts_withInfo = counts_withInfo.append(counts_layer)

# write
counts_withInfo.to_csv(output_path + "/" + sample_type + "_subsampled_matrix_max" + str(n_cells) + "cells_ENSG.txt", 
                       header=False, sep='\t')

