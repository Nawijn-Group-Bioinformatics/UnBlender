UnBlender allows respiratory scientists to perform cell type deconvolution in a validated manner. Not all cell type deconvolution approaches are feasible, especially when using highly granular (i.e. high resolution, very specific) cell type labels. UnBlender leverages the Human Lung Cell Atlas reference data [(Sikkema et al., 2023)](https://www.nature.com/articles/s41591-023-02327-2) to allow the user to specify a cell type label granularity that suits the needs of their research question, and then validates whether that approach is able to yield accurate results.

Why is this important? Because otherwise there is a very realistic risk of generating meaningless deconvolution results. For reliable results, it is essential to test whether the chosen combination of cell types *can* accurately be deconvoluted in the chosen sample type.

How does UnBlender evaluate deconvolution accuracy? In short: by deconvoluting pseudo-bulk samples with a known cell type composition. More information can be found in [TODO].

### Graphical user interface
UnBlender will soon be available as an easy to use graphical user interface (GUI), with a point-and-click menus to configure analysis parameters. The UnBlender GUI evaluates deconvolution accuracy of the specified approach, and can then be used to deconvolute your own transcriptomics dataset. The GUI can be used for deconvolution of nasal brush, bronchial brush, bronchial biopsy and parenchymal resection samples.

Output of the UnBlender GUI:
- A visual summary of expected deconvolution accuracy per cell type
- A table of deconvolution accuracy evaluation results
- Deconvolution analysis results on your transcriptomics dataset

Want to try it now? A beta version of the GUI with full functionality is available for testers, contact us at unblender.info@gmail.com or submit an issue to this repository.

### Command line interface
For those who prefer command line functionality, this GitHub repository provides a pipeline for the command line interface (CLI). The UnBlender CLI evaluates a specified deconvolution approach in the selected sample type, and returns the following:

- A summary of the expected deconvolution accuracy
- A signature matrix that can be used to deconvolute your transcriptomics dataset
- The reference data used to generate the signature matrix

Then, the user can take the signature matrix (or reference data) and perform cell type deconvolution analysis on their dataset. The signature matrix can be saved and shared for future analyses, and easily be incorporated into an existing deconvolution workflow using your tool of preference. 

#### How to use the UnBlender CLI:
1) Install CIBERSORTx (docker), and activate it.
2) Install & activate the conda environment supplied with the code.
3) Configure your pipeline run using a config.yaml file (example file is provided, see below for parameters)
4) [Download the .h5ad HLCA file (core).](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293)
4) In your command line interface, go to the directory that contains the Snakemake file.
5) Run the pipeline using the following command: 
	snakemake -c1 --configfile [config_file_name.yaml]

#### Configuration of parameters:
- config_filename: the absolute path + file name to the config.yaml file, e.g. /home/user/Documents/my_deconvolution_analysis/config.yaml
- sample_type: "parenchyma", "bronchial_brush", "nasal_brush", or "bronchial_biopsy"
- output_dir: absolute path to the output directory, e.g.  /home/user/Documents/my_deconvolution_analysis/output
- DEG_file: optional, the absolute path + file name to a file with any gene names you want to remove from the analysis (e.g. because they are differentially expressed in your bulk RNA-seq cohort). If you don't want to use this feature, give it value "no_filter". It is *not* recommended to remove large numbers of genes, as this will reduce deconvolution accuracy. However, in specific cases it may be worthwhile to remove a few. Format: .csv/.txt/.tsv file with a single gene identifier (HGNC gene names) on each new line.
- HLCA_file: path to the .h5ad HLCA file.
- email: the email address with which you run CIBERSORTx
- token: your private CIBERSORTx token
- cell_types: the cell types you wish to deconvolute, and which HLCA annotation level they correspond to. (See [HLCA supplementary table 4](https://www.nature.com/articles/s41591-023-02327-2#Sec57) for available cell types & levels.) Format: a nested YAML list. Sometimes it is useful to merge two cell types with similar gene expression profiles into one category, which is possible as shown in the example below for basal and secretory cells:

### Partial example cell type selection of a deconvolution setup into three categories:
```cell_types:
  - - "ann_level_2_clean"
    - "Lymphatic EC"

  - - "ann_level_3_clean"
    - - "Basal"
      - "Secretory"

  - - "ann_level_3_clean"
    - "Multiciliated lineage"
```

# Citing UnBlender
Have you used UnBlender (GUI or CLI) to perform a cell type deconvolution analysis? Please cite us as: [TODO]

# UnBlender is in beta version
UnBlender pipeline code is in beta version. If you encounter any bugs, please submit an issue to this repository. Or contact unblender.info@gmail.com.
(A known bug is that currently, use of certain cell types with uncommon non-alphanumeric characters in their name may cause problems running the pipeline.)
