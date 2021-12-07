args <- commandArgs(trailingOnly = T);

refdir=args[1]
datadir <- args[2]
fout <- args[3]

#source('/home/meisl/bin/FunctionLib/Lib/pagodaLib.r')

library(ggplot2)


dgeToSeurat2=function (refdir) {
  #refdir='/home/meisl/tools/slide-seq/test/PCA6T/'
  dge_file = file.path(refdir, "dge.rds")
  raw.data= readRDS(dge_file)
  
  
  meta_data = read.csv(file.path(refdir, "meta_data.csv"))
  rownames(meta_data) = meta_data$barcode
  meta_data$barcode = NULL
  common_barcodes = intersect(colnames(raw.data), rownames(meta_data))
  raw.data = raw.data[, common_barcodes]
  meta_data = meta_data[common_barcodes, ]
  cell_dict_file <- paste(refdir, "cell_type_dict.csv", sep = "/")
  true_type_names <- RCTD:::remap_celltypes(cell_dict_file, meta_data$cluster)
  meta_data$liger_ident_coarse = true_type_names
  reference = Seurat::CreateSeuratObject(raw.data, meta.data = meta_data)
  saveRDS(reference, paste(refdir, "SCRef.RDS", sep = "/"))
  dref <- RCTD:::create_downsampled_data(reference, refdir)
  return(dref)
}


read.SpatialRNA2=function (datadir, count_file = "MappedDGEForR.rds") 
{
  coords <- readr::read_csv(file = paste(datadir, "BeadLocationsForR.csv", 
                                         sep = "/"))
  counts <- readRDS(file = paste(datadir, count_file, 
                                 sep = "/"))
  colnames(coords)[2] = "x"
  colnames(coords)[3] = "y"
  #counts = tibble::column_to_rownames(counts, var = colnames(counts)[1])
  coords = tibble::column_to_rownames(coords, var = "barcodes")
  coords$barcodes <- NULL
  #counts = counts[, 2:dim(counts)[2]]
  puck = RCTD:::SpatialRNA(coords, counts)
  restrict_puck(puck, colnames(puck@counts))
}


library(RCTD)
library(Matrix)


# refdir <- '/home/meisl/tools/RCTD/Example/scRNA.refernece/'
reference <- dgeToSeurat2(refdir)

#datadir <-'/home/meisl/tools/RCTD/Example/input.slide.seq/'

puck <- read.SpatialRNA2(datadir) #


barcodes <- colnames(puck@counts) #pixels to be used (a list of barcode names). 
# This list can be restricted if you want to crop the puck e.g. 
# puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
# on the plot:
p1=plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI') 


ggsave(paste(datadir,'/',fout,'.nUMI.png',sep=''),p1)




myRCTD <- create.RCTD(puck, reference, max_cores = 10)
myRCTD <- run.RCTD(myRCTD, doublet_mode = TRUE)





results <- myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = sweep(results$weights, 1, Matrix::rowSums(results$weights), '/') 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA


resultsdir=fout

dir.create(resultsdir)

saveRDS(myRCTD,paste('./',resultsdir,'/',fout,'.res.rds',sep=''))


# make the plots
# Plots the confident weights for each cell type as in full_mode (saved as 
# 'results/cell_type_weights_unthreshold.pdf')
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
# Plots all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
# Plots the weights for each cell type as in doublet_mode. (saved as 
# 'results/cell_type_weights_doublets.pdf')
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                     results$results_df) 
# Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)


p2=plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir) 

ggsave(paste(datadir,'/',fout,'.ano.png',sep=''),p2)





# doublets
#obtain a dataframe of only doublets
doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",] 
# Plots all doublets in space (saved as 
# 'results/all_doublets.pdf')
plot_doublets(spatialRNA, doublets, resultsdir, cell_type_names) 



# Plots all doublets in space for each cell type (saved as 
# 'results/all_doublets_type.pdf')
plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names) 
# a table of frequency of doublet pairs 
doub_occur <- table(doublets$second_type, doublets$first_type) 
# Plots a stacked bar plot of doublet ocurrences (saved as 
# 'results/doublet_stacked_bar.pdf')
plot_doub_occur_stack(doub_occur, resultsdir, cell_type_names) 

