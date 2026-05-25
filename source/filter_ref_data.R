### About this script ###
# Optionally filter out a list of genes from the reference data
# that will be used to generate the signature matrix.

cat('### Filtering reference data\n')

args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = args[2]
DEG_file = args[3]

ref_data <- read.delim(input_file, header = F, row.names = 1)

if (DEG_file != 'no_filter') { # if gene list to filter out is specified
  
  DEGs <- read.delim(DEG_file, header=F)
  DEGs <- DEGs$V1
  
  ref_data_new = ref_data[!(rownames(ref_data) %in% DEGs),]
  write.table(ref_data_new, file = output_file, quote=F, sep = '\t', col.names = F)

} else { # no genes to filter? write reference data as is
  cat('no genes to filter out\n')
  write.table(ref_data, file = output_file, quote=F, sep = '\t', col.names = F)
}
cat('\n\n')