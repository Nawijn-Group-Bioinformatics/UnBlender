# UnBlender: reliable cell type deconvolution

UnBlender allows respiratory scientists to perform cell type
deconvolution with a custom, validated approach. Not all cell type
deconvolution analyses will yield reliable results, especially when
using highly granular (i.e. high resolution, very specific) cell type
labels. UnBlender leverages the Human Lung Cell Atlas [(Sikkema et al.,
2023)](https://www.nature.com/articles/s41591-023-02327-2) to perform
deconvolution into a cell type label granularity that suits the need of
a research question, and validates whether that approach is able to
yield accurate results.

<p align="center">
<img src="./UnBlender_workflow.png" alt="The UnBlender workflow" width="550">
</p>
**Why is this important?** Because otherwise there is a very realistic
risk of generating meaningless deconvolution results. For reliable
results, it is essential to test whether the chosen combination of cell
types *can* accurately be deconvoluted in the chosen sample type.

**How does UnBlender evaluate deconvolution accuracy?** In short: by
deconvoluting pseudo-bulk samples with a known cell type composition.
More information can be found in [TODO].

**Which sample types are available?** At the moment, UnBlender can be
run to evaluate deconvolution strategies for nasal brush, bronchial
brush, bronchial biopsy and parenchymal resection samples.

------------------------------------------------------------------------

### Graphical user interface

UnBlender will soon be available as an easy to use graphical user
interface (GUI), with a point-and-click menu to configure analysis
parameters. The GUI evaluates the accuracy of the chosen deconvolution
approach, and - if sufficient - deconvolutes the user's transcriptomics
dataset.

Output: 

- A summary of the expected deconvolution accuracy per cell
type: figures and a table to download for inclusion in a manuscript
supplement 
- Deconvolution analysis results of your transcriptomics
dataset

Want to use it now? A beta version of the GUI with full functionality is
available for testers, contact us at
[unblender.info\@gmail.com](mailto:unblender.info@gmail.com){.email} or
submit an issue to this repository.

------------------------------------------------------------------------

### Command line interface

This repository provides a pipeline for the command line interface
(CLI), which evaluates a chosen deconvolution approach and returns:

-   A summary of the expected deconvolution accuracy per cell type:
    figures and a table to download for inclusion in a manuscript
    supplement
-   A signature matrix to use to deconvolute your transcriptomics
    dataset
-   The reference data used to generate the signature matrix

The user can take the signature matrix (or reference data) and perform a
deconvolution analysis on their dataset. The signature matrix can be
saved and shared for future analyses, and easily be incorporated into an
existing workflow using your favourite tool.

##### How to use the UnBlender CLI:

1)  Install CIBERSORTx (docker), and activate it.
2)  Install & activate the conda environment supplied with the code.
3)  Configure your pipeline run using a config.yaml file (example file
    is provided, see below for parameters)
4)  [Download the .h5ad HLCA file
    (core).](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293)
    into the UnBlender /source folder
5)  In your command line interface, go to the directory that contains
    the Snakemake file.
6)  Run the pipeline using the following command: snakemake -c1
    --configfile [config_file_name.yaml]

###### Configuration of parameters:

-   `config_filename`: the absolute path + file name to the config.yaml
    file, e.g. `/home/user/Documents/my_deconv_analysis/config.yaml`
-   `sample_type`: `"parenchyma"`, `"bronchial_brush"`, `"nasal_brush"`,
    or `"bronchial_biopsy"`
-   `output_dir`: absolute path to the output directory, e.g.
    `/home/user/Documents/my_deconv_analysis/output`
-   `DEG_file`: optional, the absolute path + file name to a file with
    any gene names to remove from the analysis (e.g. because they are
    differentially expressed in your data). If you don't want to use
    this feature, give it value `"no_filter"`. It is *not* recommended
    to remove large numbers of genes, as this will reduce deconvolution
    accuracy. However, in specific cases it may be worthwhile to remove
    a few. Format: .csv/.txt/.tsv file with a single gene identifier
    (HGNC gene names) on each new line.
-   `email`: the email address with which you run CIBERSORTx
-   `token`: your private CIBERSORTx token
-   `cell_types`: the cell types you wish to deconvolute, and which HLCA
    annotation level they correspond to. (See "cell_type_names.tsv" file
    for an overview of available cell types, listed per level.) Format:
    a nested YAML list. Sometimes it is useful to merge two cell types
    with similar gene expression profiles into one category, which is
    possible as shown in the example below for basal and secretory
    cells:

**Partial example cell type selection of a deconvolution setup into
three categories:**

``` cell_types:
  - - "ann_level_2_clean"
    - "Lymphatic EC"

  - - "ann_level_3_clean"
    - - "Basal"
      - "Secretory"

  - - "ann_level_3_clean"
    - "Multiciliated lineage"
```

------------------------------------------------------------------------

### Citing UnBlender

Have you used UnBlender (GUI or CLI) to perform a cell type
deconvolution analysis? Please cite UnBlender as: [TODO]

### Questions?

In case of questions, or if you encounter any bugs, please submit an
issue to this repository. (Or contact
[unblender.info\@gmail.com](mailto:unblender.info@gmail.com){.email}.

------------------------------------------------------------------------

###### To do:

Future CLI updates may include: - User manual on the UnBlender
workflow + tips & tricks (need) - Clean up the pipeline to simplify
output (need) - Enable signature matrix output with ENSG gene IDs
(want) - Generalize the pipeline for other reference datasets (want) -
Add MuSiC deconvolution algorithm option (want)
