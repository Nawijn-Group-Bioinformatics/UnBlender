# UnBlender: Reliable Cell Type Deconvolution

UnBlender allows respiratory scientists to perform cell type deconvolution with a custom, validated approach. Not all cell type deconvolution analyses yield reliable results, especially when using highly granular (i.e. high resolution, very specific) cell type labels. UnBlender leverages the Human Lung Cell Atlas [(Sikkema et al., 2023)](https://www.nature.com/articles/s41591-023-02327-2) to deconvolute transcriptomics data into cell type subsets tailored to a research question, and validates whether that approach is able to yield accurate results:


<p align="center">
<img src="./UnBlender_workflow.png" alt="The UnBlender workflow" width="550">
</p>

**Why is this important?** 
Deconvolution analyses run a very realistic risk of generating meaningless results. For reliable cell type deconvolution, it is essential to test whether the chosen approach is feasible in the sample type.

**How does UnBlender evaluate accuracy?** 
In short: by deconvoluting pseudo-bulk samples with a known cell type composition. More information can be found in [TODO].

**Which sample types are available?** 
At the moment, UnBlender evaluates deconvolution strategies for nasal brush, bronchial brush, bronchial biopsy and parenchymal resection samples.

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
available for testers, contact us at unblender.info@gmail.com or
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

#### How to run the UnBlender CLI
After downloading the UnBlender CLI code:
1)  Install CIBERSORTx (docker), and activate it.
2)  Install & activate the conda environment supplied with the code.
3)  [Download the .h5ad HLCA file
    (core).](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293)
    into the UnBlender /source directory
4)  Configure your pipeline run using a config.yaml file (see below)
5)  In the command line, run the pipeline as follows from the /source directory: `snakemake -c1
    --configfile my_config_file.yaml`

##### Configuration of parameters
An example config.yaml file is provided with the code, you'll need to specify the following:
-   `config_filename`: configuration .yaml file*
-   `sample_type`: `"parenchyma"`, `"bronchial_brush"`, `"nasal_brush"`,
    or `"bronchial_biopsy"`
-   `output_dir`: the output directory*
-   `DEG_file`: `"no_filter"`, or optionally: a text file* containing genes (HGNC names, one per line) to exclude from the analysis, e.g. because they are differentially expressed in your data. Removing large numbers of genes is *not* recommended.
-   `email`: the email address with which you run CIBERSORTx
-   `token`: your private CIBERSORTx token
-   `cell_types`: the cell types to deconvolute + corresponding HLCA annotation level. (See the "cell_type_names.tsv" file for the available cell types) Format: a nested YAML list, see below.

\* Specify the absolute path

##### Example: cell type selection into three categories
The following example merges the basal and secretory epithelial cells into one shared deconvolution category. This can be helpful to improve deconvolution accuracy for cell types that have (somewhat) similar transcriptomic profiles. Note the extra dash and indentation.

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
unblender.info@gmail.com).

------------------------------------------------------------------------

###### To do:

Future CLI updates may include: - User manual on the UnBlender
workflow + tips & tricks (need) - Clean up the pipeline to simplify
output (need) - Enable signature matrix output with ENSG gene IDs
(want) - Generalize the pipeline for other reference datasets (want) -
Add MuSiC deconvolution algorithm option (want)
