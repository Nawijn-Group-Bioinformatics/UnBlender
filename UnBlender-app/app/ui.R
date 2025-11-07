#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyjqui)
library(shinyTree)

library(DT)
library(writexl)

source('config.r')
source(paste0(basedir, 'functions.r'))
source(paste0(basedir, 'help_creator_functions.r'))



#### Define the sidebar menu ####
sidebar <- dashboardSidebar(
  sidebarMenu(id="tabs", 
         menuItem("Welcome", tabName = "welcome" ),
          menuItem("Step 1: Create Reference",tabName="",startExpanded = TRUE,
            menuSubItem('Select Tissue', tabName= "select_tissue"),
            menuSubItem("Select cell types", tabName = "create_sets") ,
            menuSubItem("Remove DEGs (optional)", tabName = "remove_degs")
        ),
          
          menuItem("Step 2: Validate", tabName="",startExpanded = TRUE,
             menuSubItem("Run validation",    tabName = "current_sets" ),
             menuSubItem("Inspect output", tabName = "results_evaluation")),
            
          
          menuItem("Step 3 : Deconvolute",     tabName = "",startExpanded = TRUE,
             menuSubItem("Upload data", tabName = "upload_data"),  
             menuSubItem("Results", tabName = "user_results")),
         
         menuItem("Contact",     tabName = "contact")  
                       )
                                        
) #End of sidebar

#### Define the body ####
body <- dashboardBody( 
 # # includeCSS("www/cyberlung.css")
  tags$head(
  tags$script(src="tagger.js"),
  tags$link(rel = "stylesheet", type = "text/css", href = "cyberlung.css")
  ),
 #  ,

# #,
 # includeScript("www/cyberlung.js"),
  shinyjs::useShinyjs(),
  tabItems(
 ###### WELCOME #####
    tabItem(tabName="welcome",
      mybox(title="Welcome", status="primary", width=12, solidHeader=TRUE,
           tabsetPanel(
                     tabPanel("Intro",fluidRow(column(6,uiOutput("intro")))),
                     tabPanel("Details",fluidRow(column(8,includeMarkdown(path = './www/details.md')))
      )
    )
    ) # end of fluidRow
    ), # end of tabItem(tabName="welcome",

 ###### CREATE REFERENCE #####
 ###### Select tissue ######
  tabItem(tabName="select_tissue",
         fluidRow(
           
           column(3,
            mybox(HTML(tissue_select_help()), help='tissue_select_help',title="Select sample type",
            selectizeInput("tissue_type", label="Tissue type", choices = tissues),
            actionButton("confirm_tissue", label="Confirm"),
            uiOutput("tissue_confirm")
                  )),
           column(9,mybox(title="Sample type information"
                          ,includeMarkdown(path = './www/sample_selection.md')))
           )
  ),
 ###### Create sets ######
  tabItem(tabName="create_sets",
      tabsetPanel(id="create_sets",

      tabPanel("Select cell types",
        fluidRow(
          column(3,
           mybox(title="Select cells",
                 HTML(cell_select_help()),
                 help = 'cell_select_help',
                 uiOutput("current_tissue"),
                 uiOutput("no_tissue_error"),
                 hr(),

            # Create the tree for sample selection
            shinyTree("mytree",checkbox =T,theme="proton",themeIcons = F),
            hr(),
            actionButton("add_all_for_deconv",HTML("Add selected cells <br> to deconvolution set"))
            )
           ),

       column(9,
          fluidRow(mybox(title="Current sets for deconvolution",
            column(6,mybox("Delete selected cells from the set",br(),actionButton("delete_collections","Delete selected"))),
            column(6,mybox(
              div(style="display: inline-block;vertical-align:top; width: 350px;",
                  textInput("new_collection_name",placeholder="new name",
                             label="Merge selected rows with new name")),
              div(style="display: inline-block;vertical-align:top; width: 150px;",actionButton("merge_collections","Merge "))
              )),
        fluidRow(column(12, mybox(

          title = "Active sets",
          DT::dataTableOutput("current_collections"))))
        )
       ))
      )
       ),
      tabPanel("Extended Help",
               fluidRow(
                 column(6, mybox(title="Select cell types",includeMarkdown(path = './www/create_sets.md') ))
               )
      ))
      ),

 ###### Remove DEGS ########

  tabItem(tabName="remove_degs",
      tabsetPanel(id="remove_alldegs",   
          tabPanel("Remove DEGs",
       
        fluidRow(
          column(4,mybox(title="Input DEGs",
                         HTML(removedeg_help()), help='removedeg_help',
                   textAreaInput(inputId="inserted_degs", label="Enter your Differentially Expressed Genes", "",
                  rows=20,
                  placeholder=paste0(c("TSPAN6", "TNMD","DPM1","SCYL3","C1orf112"),collapse="\n")
                  ),
              #           verbatimTextOutput("value"),
                         actionButton('remove_degs',"Filter")
          )
          ),
          column(3, mybox(title="Matched DEGs",
                       uiOutput(outputId="matched_degs")
                 )),
          column(3, mybox(title="Unmatched DEGs",
                       uiOutput(outputId="unmatched_degs")
                 )
          )
        )),
        tabPanel("Extended Help", fluidRow(
          column(12, mybox(title="Remove custom genes",
                           includeMarkdown(path = './www/remove_degs.md')
          ))
        )),

   )),


###### VALIDATE ##########
  tabItem(tabName="current_sets",
    fluidRow(
       column(4,
        mybox(
        
           HTML(evaluation_help()),help="evaluation_help",
           uiOutput("pre_validation_summary"),
           uiOutput("pre_validation_error"),
           uiOutput("evaluation_instruction"),
           uiOutput("evaluate_button")
        ,title="Evaluation")),
      column(4,
        mybox(title="Evaluation Progress",
              HTML(progress_help()),help="progress_help",    
        textOutput("evaluation_progress"),
         ))
      )
    ),


###### Evaluate #######
    tabItem(tabName="results_evaluation",uiOutput("no_gt_eval_error"),
      tabsetPanel(
         tabPanel("Overview",
          fluidRow(
            column(6, mybox(title="Proportional Errors",
                            uiOutput("flip_mape"),
                            jqui_resizable(plotOutput("aggregate_mape_plot")))),
            column(6, mybox(title="Correlation",
                            uiOutput("flip_corplot"),
                            jqui_resizable(plotOutput("aggregate_cor_plot"))))
          ))
        ,
        tabPanel("Evaluation data",
           fluidRow(column(4,
             mybox(title="Data", DT::dataTableOutput("gt_stats")
                 )),
             column(6,
                    mybox(title ="Plot",
                         # checkboxInput("show_sample_ids","Show sample names"),
                          uiOutput("show_se"),
                          jqui_resizable(plotOutput("correlation_plot")))
             ))
        )
     )
   ),


########## Deconvolute ##########

######## Results ########

    tabItem(tabName="upload_data",
            fluidRow(column(7,
                            mybox(title='Data upload',
                                  uiOutput("pre_deconvolute_validation"),
                                  uiOutput("upload_error"),
hr(),
                                  HTML("The first column should contain valid gene symbols, the following 
                 columns should contain expression data. See the example below."),
                                  hr(),
             fileInput("user_bulk_file", "Upload your CSV file", accept=".csv"),
hr(),
HTML('You can also choose to use a preloaded example data. These data are obtained from the paper <a href ="https://pubmed.ncbi.nlm.nih.gov/31023846/">Epithelial dysregulation in obese severe asthmatics with gastro-oesophageal reflux.</a>.<br>  
This data set contains 91 samples from bronchial biopsies from moderate and sever asthma patients. Raw data are from <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76225">GEO</a>.'),
            checkboxInput("use_example_file", label="Use example data"),

            uiOutput("deconvolute_button"),
            h2("Sample of your data"),
            DT::dataTableOutput("userdata_snap"),
)),
            column(5,mybox(title="Progress",
                           uiOutput('deconv_error'),
                           textOutput("deconvolution_progress")))

            )

            ),

    tabItem(tabName = "user_results", uiOutput("no_deconvolution_results_error"),
      tabsetPanel(
       
        tabPanel("StackedBar",
                 checkboxInput("flip_stackedbar","Flip chart"),
                 jqui_resizable(plotOutput("cibersort_stackedbar"))
         )
        ,

        tabPanel("Heatmap",
               checkboxInput("flip_heatmap","Flip axis"),
               jqui_resizable(plotOutput("heatmap"))
        ),
        tabPanel("Data",
          fluidRow(column(8, mybox(title="Raw data",
                                   DT::dataTableOutput("music_results")
                                   )
                          ),
                   column(2, mybox(title = "Data handling",
                                   checkboxInput("pivot_music","Pivot data",value=FALSE),
                                   downloadButton("download_music_results", "Download"))),
          )
        ))
      ),

     tabItem(
      tabName="contact",
      fluidRow(
        mybox(title="Contact information",
        column(12,includeMarkdown(path = './www/contact_details.md'))
        )
        )
      )

      # end of tabItem
 )
) # End of body






#### Create the actual page ####
dashboardPage(
  dashboardHeader("title" = "UnBlender v0.4 (beta)"),
  
  sidebar,
  body
  
)
