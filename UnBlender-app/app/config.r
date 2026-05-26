# Get the full path of the current file
this_file <- normalizePath("config.r")

# Get the directory containing the file
current_dir <- paste0(dirname(this_file), '/')

outlier_color <- "darkgray"

basedir        <- current_dir
datadir        <- paste0(basedir, 'datafiles/')
pseudobulksdir <- paste0(datadir, 'pseudobulks/')
datadir_music  <- paste0(datadir, 'music_test/')


tissues <- setNames(obj = c("parenchyma",
                            "nasal_brush",
                            "bronchial_brush",
                            "bronchial_biopsy"), 
                    nm  = c("Parenchymal biopsy",
                            "Nasal brush",
                            "Bronchial brush", 
                            "Bronchial biopsy"))

TISSUE_MAPPINGS <- list(
  "parenchyma" = 
    list(arc         = c("parenchyma"),
         sample_type = c("donor_lung", "surgical_resection"),
         pseudo_bulk = "parenchyma.csv"
        ),
  "bronchial_brush" = 
    list(arc         = c("Distal Bronchi"),
         sample_type = c("brush"),
         pseudo_bulk = "bronchial_brush.csv"
         ),
  "nasal_brush" = 
    list(arc         = c('nose','Inferior turbinate'),
         sample_type = c("brush","scraping"),
         pseudo_bulk = "nasal_brush.csv"
          ),
  "bronchial_biopsy" = 
    list(arc         = c('airway','Intermediate Bronchi',"Trachea"),
         sample_type = c("biopsy"),
         pseudo_bulk = "bronchial_biopsy.csv"
        )
)

# SERVER SETTINGS
# Can be set to TRUE if you want to analyze test data. 
DEVELOPMENT <- FALSE
