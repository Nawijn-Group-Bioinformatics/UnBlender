# Get pseudobulks for testing the signature matrix.

# command line arguments: sample_type config_yml path/to/HLCA_file output_path/


import numpy as np
import pandas as pd
import scanpy as sc
import yaml
import sys
import random as rdm

# Set up, load HLCA data
sample_type = sys.argv[1] # tissue = "parenchyma"
if sample_type == "IPF_DEG_analysis": # specifically to run the disease DEG removal analysis:  # TODO: remove IPF pseudobulks from public version of code
    adata = sc.read("HLCA_full_v1.1_hasCountsLayer_colsRenamed.h5ad") # going to hard code the separate pseudobulk source file (extended HLCA) 
else:
    adata = sc.read(sys.argv[3]) # sc.read("../HLCA_v1.1_20220725.h5ad") 
output_path = sys.argv[4] #"../pseudobulk_analysis/"
config_yml = sys.argv[2]

# Load cell type selection directly from config file
with open(config_yml, "r") as stream:
	try:
		temp = yaml.safe_load(stream)
		include_cell_types = temp['cell_types']

	except yaml.YAMLError as exc:
		print(exc)


# Function: counts per cell type as ground truth
def write_count_file(adata, outfile, include_cell_types):
    
    samples = adata.obs['sample'].unique()
    
    with open(outfile, 'w') as f:
        f.write("sample,custom_label,number\n") # header
        
        for i in include_cell_types:
            for s in samples:
                
                if type(i[1]) == str:
                    nr = adata.obs[(adata.obs['sample'] == s) & (adata.obs[i[0]] == i[1])].shape[0]
                    f.write(s + ',' + i[1] + ',' + str(nr) + '\n')
                else:
                    nr = adata.obs[(adata.obs['sample'] == s) & (adata.obs[i[0]].isin( i[1]))].shape[0]
                    f.write(s + ',' + ' & '.join(i[1]) + ',' + str(nr) + '\n')


# Function: aggregate pseudobulk samples
def pseudobulk(adata, groupby):
    
    # groupby is a list of names of adata.obs columns used to split cells into pseudobulk samples
    # e.g. groupby = ['donor', 'celltype'] --> each pseudobulk 'sample' is the aggregation of a given celltype 
    # for a given donor
    
    # add metadata column for pseudobulk 'sample'
    first_col = True
    for col in groupby: 
        if first_col == True:
            adata.obs['pseudobulk_sample'] = adata.obs[col].astype(str)
            first_col = False
        else:
            adata.obs['pseudobulk_sample'] = adata.obs['pseudobulk_sample'] + '_' + adata.obs[col].astype(str)
    
    # get subset of raw count matrix for each 'sample' --> sum into single pseudobulk column.
    first_label = True
    for label in adata.obs['pseudobulk_sample'].unique():
        bool_mask = adata.obs['pseudobulk_sample'] == label # boolean mask for cells of pseudobulk sample


        subset_counts = adata.X[bool_mask,] # subset count matrix for pseudobulk sample

        # sum values over all cells in 'sample' for each feature:
        subset_pseudobulk = pd.DataFrame(subset_counts.sum(axis=0).transpose(), columns=[label])
    
        # merge with the other samples
        if first_label == True:
            pseudobulk_matrix = subset_pseudobulk
            first_label = False
        else:
            pseudobulk_matrix = pd.concat([pseudobulk_matrix, subset_pseudobulk], axis=1)
        
        
    # add feature names
    pseudobulk_matrix.index = adata.var.index

    return pseudobulk_matrix


# Identify samples to be used
if sample_type == "parenchyma":
	adata = adata[adata.obs['anatomical_region_coarse'] == 'parenchyma']
	rdm.seed(0)
	samples_to_include = rdm.sample(list(adata.obs.groupby(['sample']).size().index[adata.obs.groupby(['sample']).size() > 2500]), 20)
	samples_to_include = sorted(samples_to_include)

elif sample_type == "bronchial_biopsy":
	adata = adata[adata.obs['sample_type'] == 'biopsy']
	adata = adata[adata.obs['anatomical_region_coarse'].isin(['airway','Intermediate Bronchi',"Trachea"])]
	# get samples with >2500 UMIs
	samples_to_include = adata.obs.groupby(['sample']).size().index[adata.obs.groupby(['sample']).size() > 2500]

elif sample_type == "nasal_brush":
	adata = adata[adata.obs['sample_type'].isin(['brush', 'scraping'])]
	adata = adata[adata.obs['anatomical_region_coarse'].isin(['nose','Inferior turbinate'])]
	# get samples with >2500 UMIs
	samples_to_include = adata.obs.groupby(['sample']).size().index[adata.obs.groupby(['sample']).size() > 2500]

elif sample_type == "bronchial_brush":
	adata = adata[adata.obs['sample_type'] == 'brush']
	adata = adata[adata.obs['anatomical_region_coarse'] == "Distal Bronchi"]
	samples_to_include = adata.obs.groupby(['sample']).size().index # using all of them

elif sample_type == "IPF_DEG_analysis": # specifically to run the disease DEG removal analysis
    adata = adata[adata.obs['disease'] == 'pulmonary fibrosis']
    adata = adata[adata.obs['study'] == 'Banovich_Kropski_2020']
    # get samples with >2000 UMIs (this leaves me 9 samples instead of 7 if it's >2500, and that was a very arbitrary cut-off anyway)
    samples_to_include = adata.obs.groupby(['sample']).size().index[adata.obs.groupby(['sample']).size() > 2000]

else:
	sys.exit("ERROR: sample_type not valid: " + sample_type)


# Subset HLCA to those samples 
adata = adata[adata.obs['sample'].isin(samples_to_include)]

# Add the custom_label annotation for selected cell types
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
        adata_subs = adata_subset
    else:
        adata_subs = adata_subs.concatenate(adata_subset)
    
    count_subsets +=1

adata = adata_subs


# Get counts per cell type
write_count_file(adata, output_path + '/cell_counts.csv', include_cell_types)

# Make pseudobulks
pseudobulks = pseudobulk(adata, ['sample'])
pseudobulks.to_csv(output_path + '/pseudobulks.csv', sep='\t')
