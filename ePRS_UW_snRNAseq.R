library(readxl)
library(Matrix)
library(rqdatatable)
library(synapser)

#to login to synapse manually: synLogin("username", "password")
synLogin()

system( paste0('tar -xvf ', synapser::synGet('syn21614190')$path) )
system( 'gunzip home/dnanexus/6672/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6672/outs/filtered_feature_bc_matrix/features.tsv.gz')
system( 'gunzip home/dnanexus/6672/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6672 <- readMM('home_UWeprs/dnanexus/6672/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6672  <- read.delim('home_UWeprs/dnanexus/6672/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)
#the gene names file is identical for all patients
gene_metadata <- read.delim('home_UWeprs/dnanexus/6672/outs/filtered_feature_bc_matrix/features.tsv', header=FALSE)


system( paste0('tar -xvf ', synapser::synGet('syn21614188')$path) )
system( 'gunzip home/dnanexus/6687/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6687/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6687 <- readMM('home_UWeprs/dnanexus/6687/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6687  <- read.delim('home_UWeprs/dnanexus/6687/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

system( paste0('tar -xvf ', synapser::synGet('syn21614184')$path) )
system( 'gunzip home/dnanexus/6726/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6726/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6726 <- readMM('home_UWeprs/dnanexus/6726/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6726  <- read.delim('home_UWeprs/dnanexus/6726/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

system( paste0('tar -xvf ', synapser::synGet('syn21614182')$path) )
system( 'gunzip home/dnanexus/6774/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6774/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6774 <- readMM('home_UWeprs/dnanexus/6774/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6774  <- read.delim('home_UWeprs/dnanexus/6774/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

system( paste0('tar -xvf ', synapser::synGet('syn21614179')$path) )
system( 'gunzip home/dnanexus/6802/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6802/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6802 <- readMM('home_UWeprs/dnanexus/6802/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6802  <- read.delim('home_UWeprs/dnanexus/6802/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

system( paste0('tar -xvf ', synapser::synGet('syn21614178')$path) )
system( 'gunzip home/dnanexus/6829/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6829/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6829 <- readMM('home_UWeprs/dnanexus/6829/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6829  <- read.delim('home_UWeprs/dnanexus/6829/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

system( paste0('tar -xvf ', synapser::synGet('syn21614192')$path) )
system( 'gunzip home/dnanexus/6845/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6845/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6845 <- readMM('home_UWeprs/dnanexus/6845/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6845  <- read.delim('home_UWeprs/dnanexus/6845/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

system( paste0('tar -xvf ', synapser::synGet('syn21614191')$path) )
system( 'gunzip home/dnanexus/6874/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
dat_6874 <- readMM('home_UWeprs/dnanexus/6874/outs/filtered_feature_bc_matrix/matrix.mtx')
system( 'gunzip home/dnanexus/6874/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
tags_6874  <- read.delim('home_UWeprs/dnanexus/6874/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

#need to add column names and row names to each matrix; row names are the gene names, cols are barcodes + patientID
#only need column 2 of gene_meta1 (gene names)
#delete extra column in gene metadata:
gene_metadata <- subset(gene_metadata, select=c(1,2))
colnames(gene_metadata)[2] <- "gene_short_name"
#need to make duplicated gene short names unique
gene_metadata$gene_short_name <- make.names(gene_metadata$gene_short_name, unique=TRUE)

#each rowname in the matrix is a gene, from the features file:
rownames(dat_6672) <- gene_metadata[,2]
rownames(dat_6687) <- gene_metadata[,2]
rownames(dat_6726) <- gene_metadata[,2]
rownames(dat_6774) <- gene_metadata[,2]
rownames(dat_6802) <- gene_metadata[,2]
rownames(dat_6829) <- gene_metadata[,2]
rownames(dat_6845) <- gene_metadata[,2]
rownames(dat_6874) <- gene_metadata[,2]

#each column in the matrix is one of the barcodes, or cell tags:
#colnames(barcodes1)[1] <- "tags"
colnames(dat_6672) <- tags_6672[,1]
colnames(dat_6687) <- tags_6687[,1]
colnames(dat_6726) <- tags_6726[,1]
colnames(dat_6774) <- tags_6774[,1]
colnames(dat_6802) <- tags_6802[,1]
colnames(dat_6829) <- tags_6829[,1]
colnames(dat_6845) <- tags_6845[,1]
colnames(dat_6874) <- tags_6874[,1]

#append the patient id to each column name:
colnames(dat_6672) <- paste(colnames(dat_6672), "6672", sep="_")
colnames(dat_6687) <- paste(colnames(dat_6687), "6687", sep="_")
colnames(dat_6726) <- paste(colnames(dat_6726), "6726", sep="_")
colnames(dat_6774) <- paste(colnames(dat_6774), "6774", sep="_")
colnames(dat_6802) <- paste(colnames(dat_6802), "6802", sep="_")
colnames(dat_6829) <- paste(colnames(dat_6829), "6829", sep="_")
colnames(dat_6845) <- paste(colnames(dat_6845), "6845", sep="_")
colnames(dat_6874) <- paste(colnames(dat_6874), "6874", sep="_")

#join the individual patient matrices into one (they all have the same rows)
dat <- cbind2(dat_6672,dat_6687)
dat <- cbind2(dat,dat_6726)
dat <- cbind2(dat,dat_6774)
dat <- cbind2(dat,dat_6802)
dat <- cbind2(dat,dat_6829)
dat <- cbind2(dat,dat_6845)
dat <- cbind2(dat,dat_6874)
dim(dat)
head(colnames(dat))
head(rownames(dat))

#read in metadata
syn21598855 <- synapser::synGet(entity='syn21598855')
metadata <- read_excel(syn21598855$path)
head(metadata)
metadata <- as.data.frame(metadata)
names(metadata)[names(metadata) == "Sample ID"] <- "ids"
names(metadata)[names(metadata) == "Clinical DX"] <- "Clinical.DX"
names(metadata)[names(metadata) == "path DX"] <- "path.DX"
names(metadata)[names(metadata) == "Apo E"] <- "ApoE"
names(metadata)[names(metadata) == "SEX"] <- "Sex"
metadata$ids <- as.character(metadata$ids)

head(metadata)


#this count threshold can be changed
counts <- dat[rowSums(dat != 0) >= 500,]
dim(counts)

#Make the cell_data_set (CDS) object for monocle:
#cds <- new_cell_data_set(expression_matrix (counts),
#                        cell_metadata = cell_metadata,
#                       gene_metadata = gene_annotation)
#The expression value matrix must:
#have the same number of columns as the cell_metadata has rows.
#have the same number of rows as the gene_metadata has rows.
#Additionally:
#row names of the cell_metadata object should match the column names of the expression matrix.
#row names of the gene_metadata object should match row names of the expression matrix.
#one of the columns of the gene_metadata should be named "gene_short_name", which represents the gene symbol or simple name (generally used for plotting) for each gene.

counts<-counts2

Labels = data.frame(colnames(counts))
Labels$ids <- ifelse(grepl("6672", Labels$colnames.counts.)==T, "6672",
                     ifelse(grepl("6687", Labels$colnames.counts.)==T, "6687",
                            ifelse(grepl("6726", Labels$colnames.counts.)==T, "6726",
                                   ifelse(grepl("6774", Labels$colnames.counts.)==T, "6774",
                                          ifelse(grepl("6802", Labels$colnames.counts.)==T, "6802",
                                                 ifelse(grepl("6829", Labels$colnames.counts.)==T, "6829",
                                                        ifelse(grepl("6845", Labels$colnames.counts.)==T, "6845", "6874")))))))


Labels <- natural_join(Labels, metadata, 
                       by = "ids",
                       jointype = "FULL")
rownames(Labels) = colnames(counts)


head(Labels)
dim(Labels)



detach("package:synapser", unload=TRUE)

unloadNamespace("PythonEmbedInR") # must unload synapser because causes multiple definitions of S4Vectors
dat2 <- CreateSeuratObject(counts, project = "PU1", min.cells = 1, min.features = 3)
dat2
dim(dat2)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#find the percent mitochondrial genes for each cell
dat2[["percent.mt"]] <- PercentageFeatureSet(dat2, pattern = "^MT.")

# head(dat2@meta.data, 5)
# hist(dat2@meta.data$nCount_RNA)
# hist(dat2@meta.data$nFeature_RNA)
# median(dat2@meta.data$nCount_RNA)
# median(dat2@meta.data$nFeature_RNA)
# Visualize QC metrics as a violin plot
VlnPlot(dat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(dat2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(dat2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dat2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

dat3 <- subset(dat2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 10)
dim(dat3)


#use Rebecca's preprocessed cds for the UW data, which she normalized using scran and seurat:
cds_uw<-0
sce<-0
dat<-0
#p <- synapser::synGet('syn22107806')
#sce = readRDS(p$path)
rebecca_cds <- readRDS(file="~/scAD_analysis/UW_ePRS/UW_monocle_preprocessed_cds.rds")
dim(rebecca_cds)
sce <- rebecca_cds
sce= sce[,sce$ids!='6672'] #remove the single control patient for grant figures
counts = counts(sce)
counts <- counts[rowSums(counts != 0) >= 500,]
dim(counts)
gene_short_name <- data.frame(rownames(counts))
rownames(gene_short_name) <- rownames(counts)
gene_short_name$gene_short_name <- rownames(counts)
Labels <- subset(Labels, Labels$ids!='6672')
rownames(Labels)<-colnames(counts)
head(gene_short_name)
length(gene_short_name)

saveRDS(Labels, file="~/scAD_analysis/UW_ePRS/UW_Labels.rds")
saveRDS(gene_short_name, file="~/scAD_analysis/UW_ePRS/gene_short_name.rds")
Labels <- readRDS(file="~/scAD_analysis/UW_ePRS/UW_Labels.rds")
gene_short_name <- readRDS(file="~/scAD_analysis/UW_ePRS/gene_short_name.rds")


detach("package:synapser", unload=TRUE)
unloadNamespace("PythonEmbedInR") # must unload synapser because causes multiple definitions of S4Vectors

cds_uw <- new_cell_data_set(counts,
                            cell_metadata = Labels,
                            gene_metadata = gene_short_name)


cds_subset<-0
counts2<-0
cds<-0
counts<-0
sce<-0
#cds_uw<-cds_uw[genes,]
cds_uw = preprocess_cds(cds_uw, num_dim = 30,residual_model_formula_str="~PMI", norm_method="log")
plot_pc_variance_explained(cds_uw)
cds_uw = reduce_dimension(cds_uw)
cds_uw = cluster_cells(cds_uw, method="louvain")


cds_uw$Diagnosis = cds_uw$Clinical.DX
cds_uw$Sex = cds_uw$Sex
cds_uw$Samples = cds_uw$ids

cds_uw$Sex = Labels$Sex
head(cds_uw$Sex)



p1<-plot_cells(cds_uw, color_cells_by="partition",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10))+theme(legend.position = "none")
p2<-plot_cells(cds_uw, color_cells_by="Sex",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10))
p3<-plot_cells(cds_uw, color_cells_by="Diagnosis",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8))
p4<-plot_cells(cds_uw, color_cells_by="Samples",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)
p6<-plot_cells(cds_uw, color_cells_by="cluster",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm") )+theme(legend.position = "none")
#pdf(paste0("~/UW_sn_summary.pdf"))
grid.arrange(arrangeGrob(p1,p6, ncol=2),arrangeGrob(p3,p2,ncol=2),p4, heights=c(2,2,4), ncol=1)
#dev.off()

#label marker genes, defined by mathys et al
# p <- synapser::synGet('syn21618703')
# mathy_marker_genes <- read.csv(p$path,stringsAsFactors = F)
# mathy_marker_genes <- read.csv(file="~/scAD_analysis/mathys_marker_genes.csv")
# head(mathy_marker_genes)
# genes<-c()
# for (gene in unique(c(as.vector(mathy_marker_genes$gene.name),c("SYT1","SNAP25","GRIN1","GAD1","GAD2","SLC17A7","CAMK2A","NRGN","AQP4",
#                                                                 "GFAP","MBP","MOBP","PLP1","PDGFRA","VCAN","CD74","CSF1R","C3","FLT1","CLDN5")))){
#   if (gene %in% rownames(cds_uw)){
#     genes <- c(genes,which(rownames(cds_uw)==gene))
#   }
# }
# length(genes)
# 
# cds_subset = cds_uw[genes,]
# cds_subset = preprocess_cds(cds_subset, num_dim = 30,method="PCA",residual_model_formula_str="~PMI",norm_method="size_only")
# cds_subset = reduce_dimension(cds_subset)
# cds_subset = cluster_cells(cds_subset)

p1<-plot_cells(cds_subset, color_cells_by="partition",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10))+theme(legend.position = "none")
p2<-plot_cells(cds_subset, color_cells_by="Sex",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10))
p3<-plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8))
p4<-plot_cells(cds_subset, color_cells_by="Samples",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)
p6<-plot_cells(cds_subset, color_cells_by="cluster",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm") )+theme(legend.position = "none")
#pdf(paste0("/Users/relyanow/Documents/Sage/ROSMAP_sn_writeup/figures/UW_sn_summary_marker_genes.pdf"))
grid.arrange(arrangeGrob(p1,p6, ncol=2),arrangeGrob(p3,p2,ncol=2),p4, heights=c(2,2,4), ncol=1)
#dev.off()

### save progress here to avoid rebuilding monocle object if R is terminated
saveRDS(cds_subset, file = "~/scAD_analysis/UW_ePRS/UW_monocle_preprocessed_cds.rds")
cds_uw = readRDS(file = "~/scAD_analysis/UW_ePRS/UW_monocle_preprocessed_cds.rds")
cds_subset <- cds_uw


png(paste0("~/scAD_analysis/UW_ePRS/UW_sn_diffexp_markers.png"),width = 150, height = 150, units='mm', res = 300)
plot_cells(cds_subset, genes=c("SYT1","SNAP25","GRIN1","GAD1","GAD2","SLC17A7","CAMK2A","NRGN","AQP4",
                               "GFAP","MBP","MOBP","PLP1","PDGFRA","VCAN","CD74","CSF1R","C3","FLT1","CLDN5"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=.1)
dev.off()

plot_cells(cds_subset, genes=c("SYT1","SNAP25","GRIN1","GAD1","GAD2","SLC17A7","CAMK2A","NRGN","AQP4",
                               "GFAP","MBP","MOBP","PLP1","PDGFRA","VCAN","CD74","CSF1R","C3","FLT1","CLDN5"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=.1)
plot_cells(cds_subset, color_cells_by="partition",cell_size=.001,label_cell_groups=1)

#define braod cell types based on mathys marker genes:
cds_subset$broad.cell.type = cds_subset$Sex
cds_subset$broad.cell.type<-as.character(cds_subset$broad.cell.type)
for (partition in unique(partitions(cds_subset))){
  interneuron <- mean(counts(cds_subset)[rownames(cds_subset)=="GAD1",partitions(cds_subset)==partition])
  interneuron <- mean(c(interneuron,mean(counts(cds_subset)[rownames(cds_subset)=="GAD2",partitions(cds_subset)==partition])))
  Ex <- mean(counts(cds_subset)[rownames(cds_subset)=="CAMK2A",partitions(cds_subset)==partition])
  Ex <- mean(c(Ex,mean(counts(cds_subset)[rownames(cds_subset)=="NRGN",partitions(cds_subset)==partition])))
  Ast <- mean(counts(cds_subset)[rownames(cds_subset)=="AQP4",partitions(cds_subset)==partition])
  Ast <- mean(c(Ast,mean(counts(cds_subset)[rownames(cds_subset)=="GFAP",partitions(cds_subset)==partition])))/2
  Oli <- mean(counts(cds_subset)[rownames(cds_subset)=="MBP",partitions(cds_subset)==partition])
  Oli <- mean(c(Oli,mean(counts(cds_subset)[rownames(cds_subset)=="PLP1",partitions(cds_subset)==partition])))/3
  OPC <- mean(counts(cds_subset)[rownames(cds_subset)=="PDGFRA",partitions(cds_subset)==partition])
  OPC <- mean(c(OPC,mean(counts(cds_subset)[rownames(cds_subset)=="VCAN",partitions(cds_subset)==partition])))
  Mic <- mean(counts(cds_subset)[rownames(cds_subset)=="CD74",partitions(cds_subset)==partition])
  Mic <- mean(c(Mic,mean(counts(cds_subset)[rownames(cds_subset)=="CSF1R",partitions(cds_subset)==partition])))
  End <- mean(counts(cds_subset)[rownames(cds_subset)=="FLT1",partitions(cds_subset)==partition])
  End <- mean(c(End,mean(counts(cds_subset)[rownames(cds_subset)=="CLDN5",partitions(cds_subset)==partition])))
  names <- c('In','Ex','Ast','Oli','Opc','Mic','End')
  means <- c(interneuron, Ex, Ast, Oli, OPC, Mic, End)
  best_name <- names[which(means == max(means))]
  cds_subset$broad.cell.type[partitions(cds_subset)==partition] = best_name
}
plot_cells(cds_subset, color_cells_by="broad.cell.type",cell_size=.001,label_cell_groups=TRUE,show_trajectory_graph=FALSE)
for (celltype in unique(cds_subset$broad.cell.type)){
  print(celltype)
  print(length(cds_subset$broad.cell.type[cds_subset$broad.cell.type==celltype])/length(colnames(cds_subset)))
}


for (partition in unique(partitions(cds_subset))){
  interneuron <- mean(counts(cds_subset)[rownames(cds_subset)=="GAD1",partitions(cds_subset)==partition])
  interneuron <- mean(c(interneuron,mean(counts(cds_subset)[rownames(cds_subset)=="GAD2",partitions(cds_subset)==partition])))
  Ex <- mean(counts(cds_subset)[rownames(cds_subset)=="CAMK2A",partitions(cds_subset)==partition])
  Ex <- mean(c(Ex,mean(counts(cds_subset)[rownames(cds_subset)=="NRGN",partitions(cds_subset)==partition])))
  Ast <- mean(counts(cds_subset)[rownames(cds_subset)=="AQP4",partitions(cds_subset)==partition])
  Ast <- mean(c(Ast,mean(counts(cds_subset)[rownames(cds_subset)=="GFAP",partitions(cds_subset)==partition])))/2
  Oli <- mean(counts(cds_subset)[rownames(cds_subset)=="MBP",partitions(cds_subset)==partition])
  Oli <- mean(c(Oli,mean(counts(cds_subset)[rownames(cds_subset)=="PLP1",partitions(cds_subset)==partition])))/3
  OPC <- mean(counts(cds_subset)[rownames(cds_subset)=="PDGFRA",partitions(cds_subset)==partition])
  OPC <- mean(c(OPC,mean(counts(cds_subset)[rownames(cds_subset)=="VCAN",partitions(cds_subset)==partition])))
  Mic <- mean(counts(cds_subset)[rownames(cds_subset)=="CD74",partitions(cds_subset)==partition])
  Mic <- mean(c(Mic,mean(counts(cds_subset)[rownames(cds_subset)=="CSF1R",partitions(cds_subset)==partition])))
  names <- c('In','Ex','Ast','Oli','Opc','Mic')
  means <- c(interneuron, Ex, Ast, Oli, OPC, Mic)
  best_name <- names[which(means == max(means))]
  cds_subset$broad.cell.type[partitions(cds_subset)==partition] = best_name
}
for (celltype in unique(cds_subset$broad.cell.type)){
  print(celltype)
  print(length(cds_subset$broad.cell.type[cds_subset$broad.cell.type==celltype])/length(colnames(cds_subset)))
}

plot_cells(cds_subset, color_cells_by="broad.cell.type",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm") )#+theme(legend.position = "none")


# cds$broad.cell.type[cds$broad.cell.type=='Per'] ='End'
# p1<-plot_cells(cds, color_cells_by="broad.cell.type",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+ggtitle('ROSMAP')+theme(axis.title=element_text(size=10))
# p2<-plot_cells(cds_subset, color_cells_by="broad.cell.type",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+ggtitle('UW')+theme(axis.title=element_text(size=10))


# celltype <- c(rep("Ex" , 2) , rep("Oli" , 2) ,rep("In",2),rep("Ast",2),rep("Opc",2),rep("Mic",2),rep("End",2) )
# dataset <- rep(c("ROSMAP" , "UW") , 7)
# Percent <- c(.50,.39,.26,.33,.13,.13,.05,.06,.04,.04,.03,.04,.01,.01)
# data <- data.frame(celltype,dataset,Percent)
# #Turn your 'treatment' column into a character vector
# data$celltype <- as.character(data$celltype)
# #Then turn it back into a factor with the levels in the correct order
# data$celltype <- factor(data$celltype, levels=unique(data$celltype))
# # Grouped
# p3<-ggplot(data, aes(fill=dataset, y=Percent, x=celltype)) + 
#   geom_bar(position="dodge", stat="identity")+ theme_classic()+theme(text=element_text(size=21))
# pdf("/Users/relyanow/Documents/Sage/ROSMAP_sn_writeup/figures/celltypes_by_dataset_percentage.pdf") 
# grid.arrange(arrangeGrob(p1,p2,ncol=2),p3,ncol=1,heights=c(.6,1))
# dev.off()
# 
# 

#recode diagnosis to control vs ad (pt 6672 was the only control)
cds_subset$dementia = cds_subset$Diagnosis
cds_subset$Diagnosis = cds_subset$dementia
cds_subset$Diagnosis[cds_subset$ids==6672]='Control'
cds_subset$Diagnosis[cds_subset$ids!=6672]='AD'

#keep a non-manipulated cds with broad cell types labeled to store, and to retrieve:
cds_preserve <- cds_subset
#cds_uw <- cds_preserve

## Identify Mic1 subcluster in UW data
mathy_marker_genes <-read.csv(file="~/scAD_analysis/mathys_marker_genes.csv")
#l = c(l,as.vector(mathy_marker_genes$gene.name[mathy_marker_genes$subpopulation=='Mic0']))
l = as.vector(mathy_marker_genes$gene.name[mathy_marker_genes$subpopulation=='Mic1'])
#l = c(l,as.vector(mathy_marker_genes$gene.name[mathy_marker_genes$subpopulation=='Mic2']))
l = unique(l)
length(l)
inds = c()
for (gene in l){
  if ((gene %in% rownames(cds_uw))){
    inds = c(inds,which(rownames(cds_uw)==gene))
  }
}
cds_subset <- cds_uw[inds,cds_uw$broad.cell.type=='Mic']
cds_subset <- preprocess_cds(cds_subset, num_dim = 30,residual_model_formula_str="~PMI")
cds_subset = reduce_dimension(cds_subset)
cds_subset = cluster_cells(cds_subset, method="louvain")
plot_cells(cds_subset, color_cells_by="ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)


# get all Mic cells
cds_subset <- cds_preserve
cs = cds_subset[,cds_subset$broad.cell.type=='Mic']
dim(cs)
cs$Mic1 = cs$Sex
# cs$Mic1[clusters(cds_subset)==7]='Mic1'
# cs$Mic1[clusters(cds_subset)==11]='Mic1'
# cs$Mic1[clusters(cds_subset)==5]='Mic1'
cs$Mic1[log(as.vector(counts(cs)[rownames(cs)=='FTL',]))>1.5]='Mic1'
cs$Mic1[cs$Mic1!='Mic1']='Mic0'
cds_subset <- cs

# plot_cells(cds_subset, genes=c("FTL","APOE","TPT1","RPL13","SLC5A11","NLGN1"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
# 
# 
# 
# 
cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
dim(cds_subset)
cds_subset <- preprocess_cds(cds_subset, num_dim = 30,residual_model_formula_str="~PMI")
cds_subset = reduce_dimension(cds_subset, reduction_method="PCA")
cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
# 
# plot_cells(cds_subset, color_cells_by="ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
# plot_cells(cds_subset, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
# plot_cells(cds_subset, genes="FTL",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
# plot_cells(cds_subset, color_cells_by="ePRS",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
# plot_cells(cds_subset, color_cells_by='Clinical.DX',cell_size=1,label_cell_groups=1,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
# plot_cells(cds_subset, color_cells_by='Diagnosis',cell_size=1,label_cell_groups=1,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
# 
# p1<-plot_cells(cds_subset, color_cells_by="ePRS",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
# plot_cells(cds_subset, genes=c("FTL","APOE","TPT1","RPL13","SLC5A11","NLGN1"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
# p2<-plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
# p3<-plot_cells(cds_subset, color_cells_by='Clinical.DX',cell_size=1,label_cell_groups=1,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
# p4<-plot_cells(cds_subset, color_cells_by="ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
# pdf("~/scAD_analysis/UW_ePRS/UW_eprs_Mic1.pdf") 
# grid.arrange(p1,p2,p3,p4,ncol=2)
# dev.off()
# 
# #there are two partitions forming the microglia cell clusters:
# cs = cds_subset[,cds_subset$broad.cell.type=='Mic']
# cs = cs[rowSums(exprs(cs) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
# dim(cs)
# cs = reduce_dimension(cs, reduction_method="PCA")
# cs = cluster_cells(cst, cluster_method="louvain")
# 
plot_cells(cs, color_cells_by="ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cs, color_cells_by="ePRS",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)


cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset, color_cells_by="ePRS",cell_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"))

cds_subset<-order_cells(cds_subset)
plot_cells(cds_subset, color_cells_by="ePRS",cell_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

cds_subset$pseudotime = pseudotime(cds_subset)
plot_cells(cds_subset, color_cells_by="pseudotime",cell_size=2)

# cds_subset$HighLow <- ifelse(cds_subset$ePRS>0, 1, 0)
# 
# 
# pval_data <- data.frame(
#   pmi = colData(cds_subset)$PMI,
#   Diagnosis = cds_subset$Diagnosis,
#   ids = cds_subset$ids,
#   ePRS = cds_subset$ePRS)
# pval_data$pseudotime = cds_subset$pseudotime
# pval_data$HighLow <- ifelse(pval_data$ePRS>0, 1, 0)
# 
# fit <- lm(ePRS~pseudotime+pmi, data=pval_data)
# summary(fit)
# 
# fitHL <- lm(HighLow~pseudotime+pmi, data=pval_data)
# summary(fitHL)
# 
# dm <- model.matrix(~pseudotime+pmi, data=pval_data)
# fit1 <- lmFit(pval_data$ePRS,dm)
# fit2 <- eBayes(fit1)
# pval=topTable(fit2,coef=2)$adj.P.Val   
# logfc = topTable(fit2,coef=2)$logFC
# 
# stars=''
# 
# 
# pval_data$HighLow = as.factor(pval_data$HighLow)
# dp <- ggplot(pval_data, aes(x=HighLow, y=pseudotime, fill=HighLow)) + 
#   geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
# p1<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#   theme(axis.title=element_text(size=14))+
#   ggtitle(paste0('p=',formatC(pval, format = "e", digits = 2)))+
#   theme(text = element_text(size = 15))
# p2<-plot_cells(cds_subset, color_cells_by="HighLow",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)+
#   theme(legend.position = "none")+
#   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#   theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
#   ggtitle(paste0('logFC=',formatC(logfc, format = "f", digits = 3),''))+
#   theme(text = element_text(size = 15))
# pdf(paste0("~/scAD_analysis/UW_ePRS/UW_mic1_pseudotime_eprs.pdf"))
# p3<-grid.arrange(arrangeGrob(p1,p2,ncol=2,top = textGrob(paste0('Mic1',stars),gp=gpar(fontsize=18,font=2))))
# dev.off()



##################### DEG analysis: pseudotime associations  #################

gene_fits = fit_models(cds_subset, model_formula_str = "~pseudotime+PMI")
fit_coefs = coefficient_table(gene_fits)
fit_coefs  <- subset(fit_coefs, term != "(Intercept)")
fit_coefs2 <- subset(fit_coefs, term == "pseudotime")
fit_coefs2 <- subset(fit_coefs2, status == 'OK')
fit_coefs2 <- fit_coefs2[order(fit_coefs2$q_value),]
fit_coefs2 <- subset(fit_coefs2, q_value < .05)
fit_coefs2
#fit_coefs[which(fit_coefs$gene_short_name=='CST3'),]
keeps = c("gene_short_name","test_val","p_value","normalized_effect","q_value")
fit_coefs2 = fit_coefs2[keeps]
write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='~/scAD_analysis/UW_ePRS/UW_allmic_degs_pseudotime.csv')

pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_top4DEGS.pdf")
p <- plot_cells(cds_subset, genes=c("FRMD4A","HS3ST4","ARHGAP18",'SLC26A3'),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
dev.off()


#### finding genes that are differentially expressed on different paths through the trajectory ###
mic1_graphtest <- graph_test(cds_subset, neighbor_graph="principal_graph")
mic1_deg_ids <- subset(mic1_graphtest, q_value < .05)
write.csv(mic1_deg_ids,file='~/scAD_analysis/UW_ePRS/UW_allmic_graphtestdegs.csv')

#cds_top4 <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_genes,]
pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_top4_graphtestDEGS.pdf")
p <- plot_cells(cds_subset, genes=c("MAP7","CLDN11", "EDIL3", "PEX5L"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
dev.off()

top_genes1 <- c("FRMD4A", "HS3ST4", "ARHGAP18", "SLC26A3", "MAP7", "CLDN11", "EDIL3", "PEX5L")
Mic1_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_genes1,]
#pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_genesXpseudotime.pdf")
p <- plot_genes_in_pseudotime(Mic1_lineage_cds, color_cells_by="ids",min_expr=0.5, vertical_jitter=TRUE)
p
dev.off()

top_genes1 <- c("MAP7", "CLDN11")
Mic1_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_genes1,]
#pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_genesXpseudotime.pdf")
p <- plot_genes_in_pseudotime(Mic1_lineage_cds, color_cells_by="ids",min_expr=0.5)
p

#by patient:
top_genes1 <- c("FRMD4A", "HS3ST4")
cds6874 = cds_subset[,cds_subset$ids=='6874']
cds6687 = cds_subset[,cds_subset$ids=='6687']
cds6726 = cds_subset[,cds_subset$ids=='6726']
cds6774 = cds_subset[,cds_subset$ids=='6774']
cds6802 = cds_subset[,cds_subset$ids=='6802']
cds6829 = cds_subset[,cds_subset$ids=='6829']
cds6845 = cds_subset[,cds_subset$ids=='6845']
Mic1_lineage_cds6874 <- cds6874[rowData(cds6874)$gene_short_name %in% top_genes1,]
Mic1_lineage_cds6687 <- cds6687[rowData(cds6687)$gene_short_name %in% top_genes1,]
Mic1_lineage_cds6726 <- cds6726[rowData(cds6726)$gene_short_name %in% top_genes1,]
Mic1_lineage_cds6774 <- cds6774[rowData(cds6774)$gene_short_name %in% top_genes1,]
Mic1_lineage_cds6802 <- cds6802[rowData(cds6802)$gene_short_name %in% top_genes1,]
Mic1_lineage_cds6829 <- cds6829[rowData(cds6829)$gene_short_name %in% top_genes1,]
Mic1_lineage_cds6845 <- cds6845[rowData(cds6845)$gene_short_name %in% top_genes1,]

#pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_genesXpseudotime.pdf")
plot_genes_in_pseudotime(Mic1_lineage_cds6874, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(Mic1_lineage_cds6687, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(Mic1_lineage_cds6726, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(Mic1_lineage_cds6774, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(Mic1_lineage_cds6802, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(Mic1_lineage_cds6829, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(Mic1_lineage_cds6845, color_cells_by="ids",min_expr=0.5, cell_size=2)

FRMD4A_cds6874 <- cds6874[rowData(cds6874)$gene_short_name=="FRMD4A"]
FRMD4A_cds6687 <- cds6687[rowData(cds6687)$gene_short_name=="FRMD4A"]
FRMD4A_cds6726 <- cds6726[rowData(cds6726)$gene_short_name=="FRMD4A"]
FRMD4A_cds6774 <- cds6774[rowData(cds6774)$gene_short_name=="FRMD4A"]
FRMD4A_cds6802 <- cds6802[rowData(cds6802)$gene_short_name=="FRMD4A"]
FRMD4A_cds6829 <- cds6829[rowData(cds6829)$gene_short_name=="FRMD4A"]
FRMD4A_cds6845 <- cds6845[rowData(cds6845)$gene_short_name=="FRMD4A"]

plot_genes_in_pseudotime(FRMD4A_cds6874, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(FRMD4A_cds6687, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(FRMD4A_cds6726, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(FRMD4A_cds6774, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(FRMD4A_cds6802, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(FRMD4A_cds6829, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(FRMD4A_cds6845, color_cells_by="ids",min_expr=0.5, cell_size=2)

HS3ST4_cds6874 <- cds6874[rowData(cds6874)$gene_short_name=="HS3ST4"]
HS3ST4_cds6687 <- cds6687[rowData(cds6687)$gene_short_name=="HS3ST4"]
HS3ST4_cds6726 <- cds6726[rowData(cds6726)$gene_short_name=="HS3ST4"]
HS3ST4_cds6774 <- cds6774[rowData(cds6774)$gene_short_name=="HS3ST4"]
HS3ST4_cds6802 <- cds6802[rowData(cds6802)$gene_short_name=="HS3ST4"]
HS3ST4_cds6829 <- cds6829[rowData(cds6829)$gene_short_name=="HS3ST4"]
HS3ST4_cds6845 <- cds6845[rowData(cds6845)$gene_short_name=="HS3ST4"]

plot_genes_in_pseudotime(HS3ST4_cds6874, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(HS3ST4_cds6687, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(HS3ST4_cds6726, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(HS3ST4_cds6774, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(HS3ST4_cds6802, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(HS3ST4_cds6829, color_cells_by="ids",min_expr=0.5, cell_size=2)
plot_genes_in_pseudotime(HS3ST4_cds6845, color_cells_by="ids",min_expr=0.5, cell_size=2)






pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_ePRS.pdf")
p <- plot_cells(cds_subset, color_cells_by="ids",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
p
dev.off()
plot_cells(cds_subset, genes=c("FRMD4A", "HS3ST4"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)







##### oligodendrocytes
cds_subset <- cds_preserve
cs = cds_subset[,cds_subset$broad.cell.type=='Oli']
cs = cs[rowSums(exprs(cs) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
dim(cs)
cs <- preprocess_cds(cs, num_dim = 30,residual_model_formula_str="~PMI")
cs = reduce_dimension(cs, reduction_method="PCA")
cs = cluster_cells(cs, cluster_method="louvain")

plot_cells(cs, color_cells_by="ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cs, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cs, color_cells_by="cluster")#,cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)

#need to remove weird tiny clusters:
dim(cs)
cs$cluster <- clusters(cs)
cs = cs[,cs$cluster!=20]
cs = cs[,cs$cluster!=2]
cs = cs[,cs$cluster!=34]
cs = cs[,cs$cluster!=16]

dim(cs)

cs = reduce_dimension(cs, reduction_method="PCA")
cs = cluster_cells(cs, cluster_method="louvain")

plot_cells(cs, color_cells_by="ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cs, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cs, color_cells_by="ePRS",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)

cds_subset <- cs 
cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset, color_cells_by="ePRS",cell_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"))

cds_subset<-order_cells(cds_subset)
plot_cells(cds_subset, color_cells_by="ePRS",cell_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

cds_subset$pseudotime = pseudotime(cds_subset)
plot_cells(cds_subset, color_cells_by="pseudotime",cell_size=.5)
##################### DEG analysis: pseudotime associations  #################

gene_fits = fit_models(cds_subset, model_formula_str = "~pseudotime+PMI")
fit_coefs = coefficient_table(gene_fits)
fit_coefs  <- subset(fit_coefs, term != "(Intercept)")
fit_coefs2 <- subset(fit_coefs, term == "pseudotime")
fit_coefs2 <- subset(fit_coefs2, status == 'OK')
fit_coefs2 <- fit_coefs2[order(fit_coefs2$q_value),]
fit_coefs2 <- subset(fit_coefs2, q_value < .05)
fit_coefs2
#fit_coefs[which(fit_coefs$gene_short_name=='CST3'),]
keeps = c("gene_short_name","test_val","p_value","normalized_effect","q_value")
fit_coefs2 = fit_coefs2[keeps]
write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='~/scAD_analysis/UW_ePRS/UW_allOli_degs_pseudotime.csv')

pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_top4DEGS.pdf")
p <- plot_cells(cds_subset, genes=c("CRYAB","QDPR","FTH1",'GPM6B'),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
dev.off()


#### finding genes that are differentially expressed on different paths through the trajectory ###
oli_graphtest <- graph_test(cds_subset, neighbor_graph="principal_graph")
oli_deg_ids <- subset(oli_graphtest, q_value < .05)
write.csv(oli_deg_ids,file='~/scAD_analysis/UW_ePRS/UW_allOli_graphtestdegs.csv')

#cds_top4 <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_genes,]
pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_top4_graphtestDEGS.pdf")
p <- plot_cells(cds_subset, genes=c("FTH1","PLP1", "SLC26A3", "FTL"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
dev.off()

top_genes1 <- c("CRYAB", "QDPR", "FTH1", "GPM6B", "FTH1", "PLP1", "SLC26A3", "FTL")
Oli_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_genes1,]
#pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_genesXpseudotime.pdf")
p <- plot_genes_in_pseudotime(Oli_lineage_cds, color_cells_by="ids",min_expr=0.5)
p
dev.off()

top_genes1 <- c("CRYAB", "QDPR", "FTH1", "FTL")
Oli_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_genes1,]
#pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_genesXpseudotime.pdf")
p <- plot_genes_in_pseudotime(Oli_lineage_cds, color_cells_by="ids",min_expr=0.5, vertical_jitter=TRUE)
p
#dev.off()
top_genes1 <- c("CRYAB", "QDPR", "FTH1", "FTL")
cs = cds_subset[,cds_subset$broad.cell.type=='Oli']
cds_subset2 = cds_subset[,cds_subset$ids!='6672']
Oli_lineage_cds <- cds_subset2[rowData(cds_subset2)$gene_short_name %in% top_genes1,]
#pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_genesXpseudotime.pdf")
p <- plot_genes_in_pseudotime(Oli_lineage_cds, color_cells_by="ids",min_expr=0.5, vertical_jitter=TRUE)
p
#dev.off()

pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_ePRS.pdf")
p <- plot_cells(cds_subset, color_cells_by="ePRS",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
p
dev.off()

# Remove the control patient and repeat clustering:
cds_subset <- cds_preserve
cs = cds_subset[,cds_subset$broad.cell.type=='Oli']
cs = cs[,cs$ids!='6672']
dim(cs)
cs = cs[rowSums(exprs(cs) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
dim(cs)
cs <- preprocess_cds(cs, num_dim = 30,residual_model_formula_str="~PMI", norm_method="log")
cs = reduce_dimension(cs, reduction_method="PCA")
cs = cluster_cells(cs, cluster_method="louvain")

plot_cells(cs, color_cells_by="ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cs, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cs, color_cells_by="cluster")#,cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)

#need to remove two weird tiny clusters:
dim(cs)
cs$cluster <- clusters(cs)
cs = cs[,cs$cluster!=5]
cs = cs[,cs$cluster!=9]
dim(cs)
cs = reduce_dimension(cs, reduction_method="PCA")
cs = cluster_cells(cs, cluster_method="louvain")

plot_cells(cs, color_cells_by="ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cs, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cs, color_cells_by="ePRS",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)

cds_subset <- cs 
cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset, color_cells_by="ePRS",cell_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"))

cds_subset<-order_cells(cds_subset)
plot_cells(cds_subset, color_cells_by="ePRS",cell_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

cds_subset$pseudotime = pseudotime(cds_subset)
plot_cells(cds_subset, color_cells_by="pseudotime",cell_size=.5)
##################### DEG analysis: pseudotime associations  #################

gene_fits = fit_models(cds_subset, model_formula_str = "~pseudotime+PMI")
fit_coefs = coefficient_table(gene_fits)
fit_coefs  <- subset(fit_coefs, term != "(Intercept)")
fit_coefs2 <- subset(fit_coefs, term == "pseudotime")
fit_coefs2 <- subset(fit_coefs2, status == 'OK')
fit_coefs2 <- fit_coefs2[order(fit_coefs2$q_value),]
fit_coefs2 <- subset(fit_coefs2, q_value < .05)
fit_coefs2
#fit_coefs[which(fit_coefs$gene_short_name=='CST3'),]
keeps = c("gene_short_name","test_val","p_value","normalized_effect","q_value")
fit_coefs2 = fit_coefs2[keeps]
write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='~/scAD_analysis/UW_ePRS/UW_allOli_degs_pseudotime.csv')

pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_top4DEGS.pdf")
p <- plot_cells(cds_subset, genes=c("CRYAB","QDPR","FTH1",'GPM6B'),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
dev.off()


#### finding genes that are differentially expressed on different paths through the trajectory ###
oli_graphtest <- graph_test(cds_subset, neighbor_graph="principal_graph")
oli_deg_ids <- subset(oli_graphtest, q_value < .05)
write.csv(oli_deg_ids,file='~/scAD_analysis/UW_ePRS/UW_allOli_graphtestdegs.csv')

#cds_top4 <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_genes,]
pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_top4_graphtestDEGS.pdf")
p <- plot_cells(cds_subset, genes=c("FTH1","PLP1", "SLC26A3", "FTL"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
dev.off()

top_genes1 <- c("CRYAB", "QDPR", "FTH1", "GPM6B", "FTH1", "PLP1", "SLC26A3", "FTL")
Oli_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_genes1,]
#pdf("~/scAD_analysis/UW_ePRS/UW_Mic1_genesXpseudotime.pdf")
p <- plot_genes_in_pseudotime(Oli_lineage_cds, color_cells_by="ids",min_expr=0.5)
p
dev.off()





#Rebecca's DEG analysis for Mic cluster
cs = cds_preserve[,cds_preserve$broad.cell.type=='Mic']
cs$Mic1 = cs$Sex
# cs$Mic1[clusters(cds_subset)==7]='Mic1'
# cs$Mic1[clusters(cds_subset)==11]='Mic1'
# cs$Mic1[clusters(cds_subset)==5]='Mic1'
cs$Mic1[log(as.vector(counts(cs)[rownames(cs)=='FTL',]))>1.5]='Mic1'
cs$Mic1[cs$Mic1!='Mic1']='Mic0'
cds_subset$Mic1 = cs$Mic1
cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset, genes=c("FTL","APOE","TPT1","RPL13","SLC5A11","NLGN1"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p1<-plot_cells(cds_subset, color_cells_by="ePRS",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
p2<-plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
p3<-plot_cells(cds_subset, color_cells_by='Mic1',cell_size=1,label_cell_groups=1,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
p4<-plot_cells(cds_subset, color_cells_by="ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
grid.arrange(p1,p2,p3,p4,ncol=2)

cs2 <- cs[,cs$Sex=='M']
gene_fits = fit_models(cs2, model_formula_str = "~Mic1+PMI")
fit_coefs = coefficient_table(gene_fits)
fit_coefs  <- subset(fit_coefs, term != "(Intercept)")
fit_coefs2 <- subset(fit_coefs, term == "Mic1Mic1")
fit_coefs2 <- subset(fit_coefs2, status == 'OK')
fit_coefs2 <- fit_coefs2[order(fit_coefs2$q_value),]
#fit_coefs2 <- subset(fit_coefs2, q_value < .05)
fit_coefs2
fit_coefs[which(fit_coefs$gene_short_name=='CST3'),]
keeps = c("gene_short_name","test_val","p_value","normalized_effect","q_value")
fit_coefs2 = fit_coefs2[keeps]
write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='~/scAD_analysis/UW_ePRS/UW_mic1.csv')







###get pseudotime for one subcluster at a time, starting with Mic1:
cds_subset <- cds_preserve
#cds_subset = cds_subset[,cds_subset$sex=='female']
cds_subset = cds_subset[,cds_subset$sex=='male']
cds_subset = cds_subset[,cds_subset$Subcluster=='Mic1'] 

cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells


dim(cds_subset)
cds_subset = reduce_dimension(cds_subset)
cds_subset = cluster_cells(cds_subset, cluster_method="louvain")

#for males, tweak clustering algorithm for continuous pseudotime tree
#cds_subset = cluster_cells(cds_subset, cluster_method="louvain", k=30)

plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

cds_subset<-order_cells(cds_subset)
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

cds_subset$pseudotime = pseudotime(cds_subset)
plot_cells(cds_subset, color_cells_by="pseudotime",cell_size=2)

cds_subset$Diagnosis[cds_subset$Diagnosis=='Control']='Cont'
cds_subset$Diagnosis[cds_subset$Diagnosis=='Late-path AD']='Late'
cds_subset$Diagnosis[cds_subset$Diagnosis=='Early-path AD']='Early'
pval_data <- data.frame(
  batch = colData(cds_subset)$batch,
  pmi = colData(cds_subset)$pmi,
  educ = colData(cds_subset)$educ,
  Diagnosis = cds_subset$Diagnosis,
  ros_ids = cds_subset$ros_ids)
pval_data$diagnosis2 <- ifelse(pval_data$Diagnosis=='Cont',0,1)
pval_data$pseudotime = cds_subset$pseudotime

fit <- lm(diagnosis2~pseudotime+pmi+batch+educ, data=pval_data)
summary(fit)


dm <- model.matrix(~pseudotime+pmi+batch+educ, data=pval_data)
fit1 <- lmFit(pval_data$diagnosis2,dm)
fit2 <- eBayes(fit1)
pval=topTable(fit2,coef=2)$adj.P.Val   
logfc = topTable(fit2,coef=2)$logFC

stars=''
# pval_data$Diagnosis = as.factor(cds_subset$Diagnosis)
# dp <- ggplot(pval_data, aes(x=Diagnosis, y=pseudotime, fill=Diagnosis)) + 
#   geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
# p1<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#   theme(axis.title=element_text(size=14))+
#   ggtitle(paste0('p=',formatC(pval, format = "e", digits = 2)))+
#   theme(plot.title = element_text(size = 14))
# p2<-plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)+
#   theme(legend.position = "none")+
#   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#   theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
#   ggtitle(paste0('logFC=',formatC(logfc, format = "f", digits = 3),''))+
#   theme(plot.title = element_text(size = 14))
# pdf(paste0("~/scAD_analysis/figures2/mathys_Mic1_M_dx_pseudotime.pdf"))
# p3<-grid.arrange(arrangeGrob(p1,p2,ncol=2,top = textGrob(paste0('Mic1',stars),gp=gpar(fontsize=15,font=2))))
# dev.off()
# 
pval_data$Diagnosis = as.factor(cds_subset$Diagnosis)
dp <- ggplot(pval_data, aes(x=Diagnosis, y=pseudotime, fill=Diagnosis)) + 
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
p1<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(axis.title=element_text(size=14))+
  ggtitle(paste0('p=',formatC(pval, format = "e", digits = 2)))+
  theme(text = element_text(size = 15))
p2<-plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  ggtitle(paste0('logFC=',formatC(logfc, format = "f", digits = 3),''))+
  theme(text = element_text(size = 15))
pdf(paste0("~/scAD_analysis/figures2/mathys_Mic0_M_dx_pseudotime.pdf"))
p3<-grid.arrange(arrangeGrob(p1,p2,ncol=2,top = textGrob(paste0('Mic0 Males',stars),gp=gpar(fontsize=18,font=2))))
dev.off()

#DEG analysis for Mic cluster
cs = cds_preserve[,cds_preserve$broad.cell.type=='Mic']
cs$Mic1 = cs$Sex
# cs$Mic1[clusters(cds_subset)==7]='Mic1'
# cs$Mic1[clusters(cds_subset)==11]='Mic1'
# cs$Mic1[clusters(cds_subset)==5]='Mic1'
cs$Mic1[log(as.vector(counts(cs)[rownames(cs)=='FTL',]))>1.5]='Mic1'
cs$Mic1[cs$Mic1!='Mic1']='Mic0'
cds_subset$Mic1 = cs$Mic1
cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset, genes=c("FTL","APOE","TPT1","RPL13","SLC5A11","NLGN1"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p1<-plot_cells(cds_subset, color_cells_by="ePRS",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
p2<-plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
p3<-plot_cells(cds_subset, color_cells_by='Mic1',cell_size=1,label_cell_groups=1,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
p4<-plot_cells(cds_subset, color_cells_by="ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
grid.arrange(p1,p2,p3,p4,ncol=2)

cs2 <- cs[,cs$Sex=='M']
gene_fits = fit_models(cs2, model_formula_str = "~Mic1+PMI")
fit_coefs = coefficient_table(gene_fits)
fit_coefs  <- subset(fit_coefs, term != "(Intercept)")
fit_coefs2 <- subset(fit_coefs, term == "Mic1Mic1")
fit_coefs2 <- subset(fit_coefs2, status == 'OK')
fit_coefs2 <- fit_coefs2[order(fit_coefs2$q_value),]
#fit_coefs2 <- subset(fit_coefs2, q_value < .05)
fit_coefs2
fit_coefs[which(fit_coefs$gene_short_name=='CST3'),]
keeps = c("gene_short_name","test_val","p_value","normalized_effect","q_value")
fit_coefs2 = fit_coefs2[keeps]
write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='~/scAD_analysis/UW_ePRS/UW_mic1.csv')
                   