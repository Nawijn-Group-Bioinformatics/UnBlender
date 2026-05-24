# Process results CIBERSORTx on pseudobulk samples
# Compare cibersort pseudobulk results to ground truth.

rm(list=ls())
gc()

library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(ggrepel)
library(cowplot)

# Setup
args = commandArgs(trailingOnly=TRUE)
cell_counts_file = args[1]
deconv_results_file = args[2]
output_dir = args[3]
if (args[4] == 'correct_true') {
  correct = TRUE
  correct_sigmat_only = FALSE
} else if (args[4] == "correct_false") {
  correct = FALSE
} else if (args[4] == "correct_true_sigmatrix") {
  correct=TRUE
  correct_sigmat_only = TRUE
} else {
  print('ERROR: no correct_true or correct_false given')
}
ref_counts_file = args[5]
sigmat_file = args[6]


# 
# # ARGS for testing PARENCHYMA:
# cell_counts_file = "/home/tessa/Documents/20230221_Signature_matrix_validation_webtool/manuscript_analyses/output/SUBANALYSIS_corrected_or_not/parenchyma_uncorrected/cell_counts.csv"
# deconv_results_file = "/home/tessa/Documents/20230221_Signature_matrix_validation_webtool/manuscript_analyses/output/SUBANALYSIS_corrected_or_not/parenchyma_uncorrected/CIBERSORTx_output/CIBERSORTx_Results.txt"
# output_dir = "/home/tessa/Desktop/TEMP/"
# correct=F
# ref_counts_file = "/home/tessa/Documents/20230221_Signature_matrix_validation_webtool/manuscript_analyses/output/SUBANALYSIS_corrected_or_not/parenchyma_uncorrected/nasal_brush_subsampled_matrix_max200cells_HUGO.txt"
# # NEW: ADD TO CONFIG!!! AND SETUP!!
# sigmat_file = '/home/tessa/Documents/20230221_Signature_matrix_validation_webtool/manuscript_analyses/output/SUBANALYSIS_corrected_or_not/parenchyma_uncorrected/nasal_brush/CIBERSORTx_output/CIBERSORTx_nasal_brush_subsampled_matrix_max200cells_HUGO_inferred_phenoclasses.CIBERSORTx_nasal_brush_subsampled_matrix_max200cells_HUGO_inferred_refsample.bm.K999.txt'
# correct_sigmat_only = F
# 
# 
# # ARGS for testing NASAL:
# cell_counts_file = "/home/tessa/Documents/20230221_Signature_matrix_validation_webtool/manuscript_analyses/output/SUBANALYSIS_corrected_or_not/nasal_uncorrected/cell_counts.csv"
# deconv_results_file = "/home/tessa/Documents/20230221_Signature_matrix_validation_webtool/manuscript_analyses/output/SUBANALYSIS_corrected_or_not/nasal_uncorrected/CIBERSORTx_output/CIBERSORTx_Results.txt"
# output_dir = "/home/tessa/Desktop/TEMP/"
# correct=F
# ref_counts_file = "/home/tessa/Documents/20230221_Signature_matrix_validation_webtool/manuscript_analyses/output/SUBANALYSIS_corrected_or_not/nasal_uncorrected/nasal_brush_subsampled_matrix_max200cells_HUGO.txt"
# # NEW: ADD TO CONFIG!!! AND SETUP!!
# sigmat_file = '/home/tessa/Documents/20230221_Signature_matrix_validation_webtool/manuscript_analyses/output/SUBANALYSIS_corrected_or_not/parenchyma_uncorrected/nasal_brush/CIBERSORTx_output/CIBERSORTx_nasal_brush_subsampled_matrix_max200cells_HUGO_inferred_phenoclasses.CIBERSORTx_nasal_brush_subsampled_matrix_max200cells_HUGO_inferred_refsample.bm.K999.txt'
# correct_sigmat_only = F
# 

# Get ground truth proportions
gtruth <- read.table(cell_counts_file, sep=",", header = T)
count = 0
for (s in unique(gtruth$sample)) {
  subset_gt <- gtruth[gtruth$sample == s,]
  subset_gt$sample_proportions <- subset_gt$number / sum(subset_gt$number)
  
  if (count == 0) {
    gtruth_withProportions = subset_gt
  } else {
    gtruth_withProportions = rbind(gtruth_withProportions, subset_gt)
  }
  count = count + 1
} # result: gtruth_withProportions --> data frame with columns: sample, custom_label, number, sample_proportions

# Get proportions deconvolution:
deconv <- read.table(deconv_results_file, sep = "\t", header = T, row.names = 1)
num_columns = length(colnames(deconv)) -3 # columns to include in pivot longer
deconv$sample = rownames(deconv)
deconv_plottable <- pivot_longer(deconv, 1:num_columns, names_to = "custom_labels", values_to = "sample_proportions")
# result: deconv_plottable --> tibble with columns: P.value, Correlation, RMSE, sample, custom_labels, sample_proportions

# Match labels to ground truth formatting
deconv_plottable$custom_labels <- gsub("\\.\\.\\.", " & ", deconv_plottable$custom_labels) 
deconv_plottable$custom_labels <- gsub("\\.", " ", deconv_plottable$custom_labels) 
deconv_plottable$custom_labels <- gsub("Hillock like", "Hillock-like", deconv_plottable$custom_labels)
deconv_plottable$custom_labels <- gsub("CCL3 ", "CCL3+", deconv_plottable$custom_labels)
deconv_plottable$custom_labels <- gsub(" non-nasal ", "(non-nasal)", deconv_plottable$custom_labels)
deconv_plottable$custom_labels <- gsub(" nasal ", "(nasal)", deconv_plottable$custom_labels)
deconv_plottable$custom_labels <- gsub(" bronchial ", "(bronchial)", deconv_plottable$custom_labels)
deconv_plottable$custom_labels <- gsub(" subsegmental ", "(subsegmental)", deconv_plottable$custom_labels)
deconv_plottable$custom_labels <- gsub("MT positive", "MT-positive", deconv_plottable$custom_labels)
deconv_plottable$custom_labels <- gsub("Monocyte derived", "Monocyte-derived", deconv_plottable$custom_labels)
deconv_plottable$custom_labels <- gsub("pre TB", "pre-TB", deconv_plottable$custom_labels)
deconv_plottable$custom_labels <- gsub("2 ", "2_", deconv_plottable$custom_labels)
deconv_plottable$custom_labels <- gsub("4 ", "4_", deconv_plottable$custom_labels)
deconv_plottable$custom_labels <- gsub("3 ", "3_", deconv_plottable$custom_labels)

# Prepare to plot
gtruth_withProportions$which_plot <- "Ground truth"
deconv_plottable$which_plot <- "Deconvolution"
names(gtruth_withProportions)[names(gtruth_withProportions) == "custom_label"] <- "custom_labels"
df_to_plot <- rbind(gtruth_withProportions[, c("sample", "custom_labels", "sample_proportions", "which_plot")], 
                    deconv_plottable[, c("sample", "custom_labels", "sample_proportions", "which_plot")])
# result: df_to_plot --> data frame with columns: sample, custom_labels, sample_proportions, which_plot

# get proportional error
gtruth_withProportions$sample_label <- paste0(gtruth_withProportions$sample, "_", gtruth_withProportions$custom_labels)
deconv_plottable$sample_label <- paste0(deconv_plottable$sample, "_", deconv_plottable$custom_labels)
prop_error_df <- full_join(gtruth_withProportions, deconv_plottable, by = "sample_label")
prop_error_df <- prop_error_df[!(is.na(prop_error_df$P.value)),]
prop_error_df <- prop_error_df[, c("sample.x", "custom_labels.x", "sample_proportions.x", "sample_proportions.y")]
prop_error_df$prop_error <- (prop_error_df$sample_proportions.x - prop_error_df$sample_proportions.y) / prop_error_df$sample_proportions.x
# result: prop_error_df --> coluns: sample.x, custom_labels.x, sample_proportions.x, sample_proportions.y, prop_error
# for x = ground truth, y = deconvolution proportions


# get correlations of deconvoluted cell type proportions with ground truth
df_corr_results = data.frame('cell type' = character(), 'estimate' = numeric(), 'p-value' = numeric(),
                             'avg % ground truth' = numeric(), 'stdev % ground truth' = numeric())

for (celltp in unique(df_to_plot$custom_labels)) {
  truth_temp = df_to_plot[df_to_plot$which_plot == "Ground truth",]
  truth = truth_temp[truth_temp$custom_labels == celltp, c('sample', 'sample_proportions')]
  
  if (correct == TRUE) {
    pred_temp = df_to_plot[df_to_plot$which_plot == "Corrected deconvolution",]
  } else {
    pred_temp = df_to_plot[df_to_plot$which_plot == "Deconvolution",]
  }
  pred = pred_temp[pred_temp$custom_labels == celltp, c('sample', 'sample_proportions')]
  
  corr_df <- full_join(truth, pred, by = "sample")
  corr_df[["sample_proportions.y"]][is.na(corr_df[["sample_proportions.y"]])] <- 0

  # Tukey fence for ground truth sample proportion outliers
  quant <- quantile(corr_df$sample_proportions.x, probs=c(.25, .75))
  H <- 1.5 * IQR(corr_df$sample_proportions.x)
  corr_df = corr_df[corr_df$sample_proportions.x > (quant[1] - H),]
  corr_df = corr_df[corr_df$sample_proportions.x < (quant[2] + H),]

  # test for correlations
  corr_test <- cor.test(corr_df$sample_proportions.x, corr_df$sample_proportions.y)
  
  res <- data.frame(celltp, corr_test$estimate, corr_test$p.value, mean(corr_df$sample_proportions.x), sd(corr_df$sample_proportions.x))
  df_corr_results <- rbind(df_corr_results, res)
}

# MAPE per cell type
MAPE_cell_type = data.frame('cell type' = character(), 'MAPE' = numeric())
for (ct in unique(prop_error_df$custom_labels.x)) {
  temp = prop_error_df[prop_error_df$custom_labels.x == ct,]
  PE = temp$prop_error[is.finite(temp$prop_error)]
  MAPE = mean(abs(PE))
  
  MAPE_cell_type = rbind(MAPE_cell_type, data.frame('cell type' = ct, 'MAPE' = MAPE))
}


############# Generate plots ###############

###### Stacked bar chart ########

# colours to plot
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# plot ground truth and deconvolution proportions (stacked bar chart)
plot_proportions <- function(df, plot_title) {
  ggplot(df, aes(x = factor(sample) ,y = sample_proportions, fill = factor(custom_labels))) + 
    geom_bar(stat = "identity") +
    theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=8), #change legend text font size
          axis.text.x = element_text(size=8, angle=90)) +
    labs(title=plot_title, x ="Sample", y = "Cell type proportions", fill = "annotation") +
    facet_wrap(~which_plot, nrow =2) + scale_fill_manual(values = col_vector)
}

stacked_bars <- plot_proportions(df_to_plot, paste("Pseudobulk composition"))

##### MAPE per cell type #######

MAPE_plot <- ggplot(MAPE_cell_type, aes(cell.type, MAPE)) + 
  ylim(0,NA) + 
  annotate(geom='rect', ymin=-Inf, ymax=1, xmin = -Inf, xmax = Inf, 
           fill='darkolivegreen2', alpha=0.5) +
  annotate(geom='rect', ymin=1, ymax=Inf, xmin = -Inf, xmax = Inf, 
           fill='sandybrown', alpha=0.5) +
  
  geom_point(size=3) + coord_flip()+
  xlab("Cell type") + ylab("Mean absolute proportional error\nof estimated vs. true proportion")


###### Correlation per cell type ######

corr_plot <- ggplot(df_corr_results, aes(celltp, corr_test.estimate)) +
  ylim(0,NA) + 
  annotate(geom='rect', ymin=0.7, ymax=Inf, xmin = -Inf, xmax = Inf, 
           fill='darkolivegreen2', alpha=0.5) +
  annotate(geom='rect', ymin=-Inf, ymax=0.7, xmin = -Inf, xmax = Inf, 
           fill='sandybrown', alpha=0.5) + coord_flip()+
  geom_point(size=3) +
  xlab("Cell type") + ylab("Correlation: estimated vs. true\ncell type proportions")


##### Ground truth vs estimated % plot #####

gt_vs_est_plot <- ggplot(prop_error_df, aes(x=sample_proportions.x, y=sample_proportions.y)) +
  geom_abline(slope=1, intercept = 0, color='chartreuse3') +
  geom_point() + facet_wrap(~custom_labels.x, scales='free') +
  xlab("Ground truth: true sample proportions") + ylab("Deconvolution: estimated proportions") +
  ggtitle("True vs. estimated proportions (green line indicates values where estimated = true)")
gt_vs_est_plot

# Putting things together
title <- ggdraw() + draw_label("Pseudobulk deconvolution accuracy per cell type")
acc_plot <- plot_grid(title, plot_grid(MAPE_plot, corr_plot), rel_heights = c(0.1, 1), ncol=1)

#### write output ####

png(filename=paste0(output_dir, '/pseudobulk_composition.png'))
print(stacked_bars)
dev.off()

png(filename=paste0(output_dir, '/deconvolution_accuracy.png'), width=800)
print(acc_plot)
dev.off()

png(filename=paste0(output_dir, '/estimated_vs_true_proportions.png'), width=700, height=700)
print(gt_vs_est_plot)
dev.off()

output_report = merge(MAPE_cell_type, df_corr_results, by.x='cell.type', by.y='celltp')
colnames(output_report) = c('Cell type', 'MAPE', 'Correlation', 'Corr. p-value', 'Mean pseudobulk proportion', 'Stdev. pseudobulk proportions')
write.table(output_report, sep=',', file = paste0(output_dir, '/deconvolution_accuracy_report.csv'), row.names = F)


# #### For figure S1 biorXiv ####
# 
# 
# 
# # plot stacked bar chart
# plot_proportions2 <- function(df, plot_title) {
#   ggplot(df, aes(x = factor(sample) ,y = sample_proportions, fill = factor(custom_labels))) + 
#     geom_bar(stat = "identity") +
#     theme(strip.text = element_text(size=8),
#           legend.key.size = unit(0.15, 'cm'), #change legend key size
#           title = element_text(size=8), #change legend title font size
#           legend.text = element_text(size=8), #change legend text font size
#           axis.text.x = element_text(size=6, angle=45, vjust = 1, hjust=1)) +
#     labs(title=plot_title, x ="Sample", y = "Cell type proportions", fill = " ") +
#     facet_wrap(~which_plot, nrow =1) + scale_fill_manual(values = col_vector)
# }
# stacked_bars2 <- plot_proportions2(df_to_plot, paste("Pseudobulk composition"))
# stacked_bars2
# 
# # MAPE
# MAPE_plot2 <- ggplot(MAPE_cell_type, aes(cell.type, MAPE)) + 
#   ylim(0,NA) + 
#   annotate(geom='rect', ymin=-Inf, ymax=1, xmin = -Inf, xmax = Inf, 
#            fill='darkolivegreen2', alpha=0.5) +
#   annotate(geom='rect', ymin=1, ymax=Inf, xmin = -Inf, xmax = Inf, 
#            fill='sandybrown', alpha=0.5) +
#   geom_point(size=1) + coord_flip()+
#   xlab("Cell type") + ylab("Mean absolute proportional error\nof estimated vs. true proportion") +
#   theme(axis.title.x = element_text(size=8),
#         axis.title.y = element_blank())
# 
# 
# # Correlation per cell type
# corr_plot2 <- ggplot(df_corr_results, aes(celltp, corr_test.estimate)) +
#   ylim(0,NA) + 
#   annotate(geom='rect', ymin=0.7, ymax=Inf, xmin = -Inf, xmax = Inf, 
#            fill='darkolivegreen2', alpha=0.5) +
#   annotate(geom='rect', ymin=-Inf, ymax=0.7, xmin = -Inf, xmax = Inf, 
#            fill='sandybrown', alpha=0.5) + coord_flip()+
#   geom_point(size=1) +
#   xlab("Cell type") + ylab("Correlation of estimated vs.\ntrue cell type proportions") +
#   theme(axis.title.x = element_text(size=8),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank()
#         )
# 
# # Ground truth vs estimated % plot
# gt_vs_est_plot2 <- ggplot(prop_error_df, aes(x=sample_proportions.x, y=sample_proportions.y)) +
#   geom_abline(slope=1, intercept = 0, color='chartreuse3') +
#   geom_point(size=0.8) + facet_wrap(~custom_labels.x, scales='free') +
#   xlab("Ground truth: true sample proportions") + ylab("Deconvolution: estimated proportions") +
#   ggtitle("True vs. estimated proportions\n(Green line indicates values where estimated equals true)") +
#   theme(strip.text = element_text(size=8),
#         title = element_text(size=8),
#         axis.text.x = element_text(size=7),
#         axis.text.y = element_text(size=7))
# gt_vs_est_plot2
# 
# # Putting accuracy plot together
# title <- ggdraw() + draw_label("Pseudobulk deconvolution accuracy per cell type", size = 10)
# acc_plot2 <- plot_grid(title, 
#                        plot_grid(NULL, MAPE_plot2, corr_plot2, NULL, rel_widths = c(0.2,2,1.1, 0.2), nrow=1), 
#                        rel_heights = c(0.1, 1), ncol=1)
# 
# 
# plot_grid(plot_grid(NULL, stacked_bars2, ncol=2, rel_widths = c(0.1,1)), acc_plot2, gt_vs_est_plot2, ncol=1, rel_heights = c(1.2,0.8,2),
#           labels = LETTERS[1:3])
