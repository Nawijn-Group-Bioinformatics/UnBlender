# Functions that create a toggab;e help icon in the headers of boxes. 
# These can be extended in the event that more help is needed.



make_help<- function(id, text){
  paste0(
    '<div id="',id,'"style="display:none;"  onClick=toggle("',
    id,'")><i class="close_cross fa fa-window-close" aria-hidden="true"></i><div class="help_box">',text,'</div></div>'
  )
}




tissue_select_help<-function(){
  make_help('tissue_select_help', "Select the tissue that most closely resembles your input tissue." 
  ) 
}

cell_select_help<-function(){
  make_help('cell_select_help', "Select cells by selecting a number of cells 
            or categories from the tree. After the selection, press the button to confirm the selection"
  ) 
}


removedeg_help<-function(){
  make_help('removedeg_help', "Here you can paste your gene names
  that you want to remove from the analysis")
}

evaluation_help<-function(){
  make_help('evaluation_help', 
            "Pressing the button <b>Evaluate</b> will start the evaluation of your selected tissue and cell types."
  )}

progress_help<-function(){
  make_help('progress_help', "In this window the progress of the analysis is shown. For large data 
  sets this may take up to a few minutes, so please be patient.... "
  )}

