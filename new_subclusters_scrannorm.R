library(scran)
library(scater)
library(Matrix)
synapser::synLogin()

#counts <- readRDS(file="~/scAD_analysis/mathys_counts_matrix.rds")
Labels <- readRDS(file="~/scAD_analysis/mathys_Labels.rds")
gene_short_name <- readRDS(file="~/scAD_analysis/mathys_gene_short_name.rds")


# 
# sce <- SingleCellExperiment(list(counts=counts))
# dim(sce)
# clusters <- quickCluster(sce, min.size=100)
# sce <- computeSumFactors(sce, cluster=clusters)
# #sce <- computeSumFactors(sce, cluster=clusters, BPPARAM=MulticoreParam(8))
# summary(sizeFactors(sce))
# 
# sce <- logNormCounts(sce)
# head(logcounts(sce[0:20,0:20]))
# head(rownames(logcounts(sce)))
# 
# dim(logcounts(sce))
# counts2 = logcounts(sce)

#scran normalized count matrix:
#counts2 <- readRDS(file="mathys_scran_normalized_counts2.rds")
p <- synapser::synGet('syn22363128')
counts <- readRDS(p$path)
dim(counts)
head(counts)
#extract count matrix from the cds
#counts <- exprs(cds)
dim(counts)

cds <- new_cell_data_set(counts,
                         cell_metadata = Labels,
                         gene_metadata = gene_short_name)

saveRDS(cds, "~/scAD_analysis/Mathys_scrannorm_cds.rds")
cds <- readRDS(file="~/scAD_analysis/Mathys_scrannorm_cds.rds")

cds$educ = as.numeric(cds$educ)
cds$Education = cds$educ
cds$Sex = cds$sex
cds$CERAD = cds$ceradsc
cds$Braak = cds$braaksc
cds$APOE_genotype = cds$apoe_genotype
cds$batch = as.factor(cds$batch)


#Preprocessing the data using PCA. Data is already size-normalized and log transformed.
cds <- clear_cds_slots(cds)
cds = preprocess_cds(cds, num_dim = 50,method="PCA", norm_method="none")
plot_pc_variance_explained(cds)

#adjust for batch effects
cds <- align_cds(cds, 
                 preprocess_method="PCA",
                 alignment_group="batch",
                 residual_model_formula_str="~educ+pmi")
cds = reduce_dimension(cds)
cds = cluster_cells(cds, cluster_method="louvain")
plot_cells(cds, color_cells_by="partition", cell_size=.1,label_cell_groups=1,show_trajectory_graph=FALSE)#+theme(
#legend.title = element_text(size = 10),
#legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
#legend.key.width = unit(0.5,"cm"))
plot_cells(cds, color_cells_by="broad.cell.type",cell_size=.1,label_cell_groups=1,show_trajectory_graph=FALSE)
plot_cells(cds, color_cells_by="Subcluster",cell_size=.1,label_cell_groups=1,show_trajectory_graph=FALSE)
plot_cells(cds, color_cells_by="cluster",cell_size=.1,label_cell_groups=1,show_trajectory_graph=FALSE)

markergenes <- c("GAD1", "GAD2", "CAMK2A", "NRGN", "AQP4", "GFAP", "MBP", "PLP1", "PDGFRA", "VCAN", "CD74", "CSF1R", "FLT1", "CLDN5", "PDGFRB", "ZIC1", "C3", "CSPG4", "MOBP", "SLC17A7")
plot_genes_by_group(cds,
                    markergenes,
                    group_cells_by="partition",
                    ordering_type="cluster_row_col",
                    max.size=8)

#Examine expression of marker genes for each partition and assign broad cell types as follows:
# braod cell type Ast: marker genes = AQP4, GFAP
# Ex: NRGN, CAMK2A, SLC17A7
# In: GAD2, GAD1
# Mic: CD74, CSF1R, C3
# Oli: MBP, MOBP, PLP1
# Opc: VCAN, PDGFRA, CSPG4
# End: FLT1, CLDN5
# Per: ZIC1, PDGFRB, AMBP
colData(cds)$broad.cell.type2 <- as.character(partitions(cds))
colData(cds)$broad.cell.type2 = dplyr::recode(colData(cds)$broad.cell.type2,
                                              "1"="Ex",
                                              "2"="Ex",
                                              "3"="Ex",
                                              "4"="Ex",
                                              "5"="Ex",
                                              "6"="Ex",
                                              "7"="Oli",
                                              "8"="In",
                                              "9"="In",
                                              "10"="In",
                                              "11"="unknown",
                                              "12"="Mic",
                                              "13"="Opc",
                                              "14"="Ast",
                                              "15"="Per",
                                              "16"="Per",
                                              "17"="End",)
plot_cells(cds, color_cells_by="broad.cell.type2",cell_size=.1,label_cell_groups=1,show_trajectory_graph=FALSE)
plot_cells(cds, color_cells_by="broad.cell.type",cell_size=.1,label_cell_groups=1,show_trajectory_graph=FALSE)
plot_cells(cds, color_cells_by="partition",cell_size=.1,label_cell_groups=1,show_trajectory_graph=FALSE)

plot_cells(cds, color_cells_by="broad.cell.type2",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

celltypes <- data.frame(
  broad.cell.type = colData(cds)$broad.cell.type,
  broad.cell.type2 = colData(cds)$broad.cell.type2,
  Subcluster = colData(cds)$Subcluster,
  ros_ids = colData(cds)$ros_ids
)
celltypes$broad.cell.type <- as.character(celltypes$broad.cell.type)
mismatch <- celltypes[celltypes$broad.cell.type!=celltypes$broad.cell.type2,]
colData(cds)$concordance <- ifelse(colData(cds)$broad.cell.type==colData(cds)$broad.cell.type2, "match", "mismatch")
table(cds$concordance)

plot_cells(cds, color_cells_by="concordance",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))


################ find top genes within partitions in each broad cell type to classify subclusters #####################
# Start with microglia

cds_mic = cds[,cds$broad.cell.type2=='Mic']
plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#cds_mic <- clear_cds_slots(cds_mic)
#cds_mic = preprocess_cds(cds_mic, num_dim = 30,method="PCA", norm_method="none")
#adjust for batch effects
cds_mic$batch <- as.character(cds_mic$batch)
cds_mic$batch <- as.factor(cds_mic$batch)
cds_mic <- align_cds(cds_mic, 
                     preprocess_method="PCA",
                     alignment_group="batch",
                     residual_model_formula_str="~educ+pmi")
cds_mic = reduce_dimension(cds_mic)
cds_mic = cluster_cells(cds_mic, cluster_method="louvain", k=500)
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=1,show_trajectory_graph=FALSE)#+theme(
plot_cells(cds_mic, color_cells_by="ros_ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#Choose k=500, which yields 4 microglia clusters
cds_mic$cluster <- clusters(cds_mic)
table(cds_mic$cluster)
cds_mic$Subcluster2[cds_mic$cluster==1] <- 'Mic1'
cds_mic$Subcluster2[cds_mic$cluster==2] <- 'Mic0'
cds_mic$Subcluster2[cds_mic$cluster==3] <- 'Mic3'
cds_mic$Subcluster2[cds_mic$cluster==4] <- 'Mic2'
plot_cells(cds_mic, color_cells_by="Subcluster2",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#extract the new subcluster labels
Mic_subs <- data.frame(
  TAG = colData(cds_mic)$TAG,
  Subcluster2 = colData(cds_mic)$Subcluster2
)
write.csv(Mic_subs, file="~/scAD_analysis/Mic_subs.csv")

#find top marker genes for each subcluster
# marker_test_res <- top_markers(cds_mic, group_cells_by="Subcluster2", cores=8)
# cds_mic2 <- cds_mic[,cds_mic$Subcluster!='Opc1']
# cds_mic2 <- cds_mic2[,cds_mic2$Subcluster!='Mic2']
# marker_test_res <- top_markers(cds_mic2, group_cells_by="Subcluster", cores=8)
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.25) %>%
#   group_by(cell_group) %>%
#   top_n(5, marker_test_q_value)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# plot_genes_by_group(cds_mic,
#                     top_specific_marker_ids,
#                     group_cells_by="Subcluster2",
#                     ordering_type="cluster_row_col",
#                     max.size=5)

cds_subset <- cds_mic[,cds_mic$sex=='female'] 
cds_subset <- cds_subset[,cds_subset$Subcluster2=='Mic0']
#cds_subset <- clear_cds_slots(cds_subset)
cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,]
cds_subset$batch <- as.character(cds_subset$batch)
cds_subset$batch <- as.factor(cds_subset$batch)
cds_subset <- align_cds(cds_subset, 
                        preprocess_method="PCA",
                        alignment_group="batch",
                        residual_model_formula_str="~educ+pmi")
cds_subset = reduce_dimension(cds_subset)
plot_cells(cds_subset, color_cells_by="Subcluster", cell_size=1)
cds_subset = cluster_cells(cds_subset, cluster_method="louvain")

# if (length(unique(cds_subset@clusters$UMAP$partitions))>1){
#   m<-0
#   p<-"1"
#   for (partition in unique(cds_subset@clusters$UMAP$partitions)){
#     if (length(cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions==partition])>m){
#       m<-length(cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions==partition])
#       p<-partition
#     }
#   }
#   cds_subset = cds_subset[,cds_subset@clusters$UMAP$partitions==p]
#   cds_subset = reduce_dimension(cds_subset)
#   cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
# }

cds_subset <- learn_graph(cds_subset)

plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2,label_cell_groups=0)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))


#The function below will pick the root at the "beginning" of the trajectory by first grouping the cells according to which trajectory graph node they are nearest to. It calcuates what fraction of the cells at each node come from Control samples. Then it picks the node that is most heavily occupied by Control cells and returns that as the root.

#After setting the root, call the order_cells function to calculate pseudotime based on order, and plot by pseudotime:


get_earliest_principal_node <- function(cds_subset, diagnosis="Control"){
  cell_ids <- which(colData(cds_subset)[, "Diagnosis"]==diagnosis)
  
  closest_vertex <-
    cds_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds_subset), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds_subset)[["UMAP"]])$name[as.numeric(names
                                                                     (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
  
}
cds_subset <- order_cells(cds_subset, root_pr_nodes=get_earliest_principal_node(cds_subset))
cds_subset$pseudotime = pseudotime(cds_subset)

plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2,label_cell_groups=0)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds_subset, color_cells_by="ros_ids",cell_size=2,label_cell_groups=0)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#Create a dataframe that will contain variables necessary to regress: diagnosis ~ pseudotime + educ + pmi
#Note: we are not adjusting for batch because the align_cds function already effectively subtracted out batch effects and we do not want to overcorrect. the align_cds residual_model_formula_str function only miniimally affects the clustering of the cells, so these variables should be included.


pval_data <- data.frame(
  pmi = colData(cds_subset)$pmi,
  educ = colData(cds_subset)$educ,
  Diagnosis = cds_subset$Diagnosis,
  ros_ids = cds_subset$ros_ids)
pval_data$diagnosis2 <- ifelse(pval_data$Diagnosis=='Control',0,1)
pval_data$pseudotime = cds_subset$pseudotime
fit <- lm(diagnosis2~pseudotime+pmi+educ, data=pval_data)
summary(fit)

dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
fit1 <- limma::lmFit(pval_data$diagnosis2,dm)
fit2 <- limma::eBayes(fit1)
pval=limma::topTable(fit2,coef=2)$adj.P.Val   
logfc = limma::topTable(fit2,coef=2)$logFC
stars=''
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
p4<-plot_cells(cds_subset, color_cells_by="ros_ids",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE)+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())
#pdf(paste0("~/scAD_analysis/Figures_LH/F_Ast0_dx_pseudotime.pdf"))
p <- grid.arrange(arrangeGrob(p1,p2,p4,ncol=3,top = textGrob(paste0('Ast0 Female',stars),gp=gpar(fontsize=18,font=2))))
p
#dev.off()


#We want to see which genes change expression across pseudotime by performing regression using the Monocle3 function fit_models, which fits the following regression for each gene: expr~pseudotime+pmi+educ
#We will output a csv file with the pseudotime coefficients and p-values for all genes with qvalue<0.05.
gene_fits = fit_models(cds_subset, model_formula_str = "~pseudotime+pmi+educ")
fit_coefs = coefficient_table(gene_fits)
fit_coefs  <- subset(fit_coefs, term != "(Intercept)")
fit_coefs2 <- subset(fit_coefs, term == "pseudotime")
fit_coefs2 <- subset(fit_coefs2, status == 'OK')
fit_coefs2 <- fit_coefs2[order(fit_coefs2$q_value),]
fit_coefs2 <- subset(fit_coefs2, q_value < .05)
keeps = c("gene_short_name","test_val","p_value","normalized_effect","q_value")
fit_coefs2 = fit_coefs2[keeps]
#write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='~/scRNAseq-workshop2020/Oli1_F_genes.csv')
head(fit_coefs2)

#Pull the top 6 most significant genes that are positively associated with pseudotime (i.e. increasing expression through pseudotime) 


up_oli1 <- subset(fit_coefs2, fit_coefs2$test_val>0)
top_up_oli1 <- up_oli1 %>%
  top_n(-6, q_value)
genelist_up_oli <- unique(top_up_oli1 %>% pull(gene_short_name))
plot_cells(cds_subset, genes=genelist_up_oli,cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)

#Repeat for the top six most significant genes that are negatively associated with pseudotime (i.e. decreasing expression through pseudotime)


down_oli1 <- subset(fit_coefs2, fit_coefs2$test_val<0)
top_down_oli1 <- down_oli1 %>%
  top_n(-6, q_value)
genelist_down_oli <- unique(top_down_oli1 %>% pull(gene_short_name))
plot_cells(cds_subset, genes=genelist_down_oli,cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)

oli1_up_plotgenes <- cds_subset[rowData(cds_subset)$gene_short_name %in% genelist_up_oli,]
plot_genes_in_pseudotime(oli1_up_plotgenes, color_cells_by="Diagnosis",min_expr=0.5)

#Plot pseudotime across the x axis and expression on the y axis for each of the top 6 genes that are negatively associated with pseudotime, using the Monocle3 function plot_genes_in_pseudotime 

oli1_down_plotgenes <- cds_subset[rowData(cds_subset)$gene_short_name %in% genelist_down_oli,]
plot_genes_in_pseudotime(oli1_down_plotgenes, color_cells_by="Diagnosis",min_expr=0.5)




################ find top genes within partitions in each broad cell type to classify subclusters #####################
#Excitatory neurons
cds <- cds_all
cds_mic = cds[,cds$broad.cell.type2=='Ex']
plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#cds_mic <- clear_cds_slots(cds_mic)
#cds_mic = preprocess_cds(cds_mic, num_dim = 30,method="PCA", norm_method="none")
#adjust for batch effects
cds_mic$batch <- as.character(cds_mic$batch)
cds_mic$batch <- as.factor(cds_mic$batch)
cds_mic <- align_cds(cds_mic, 
                     preprocess_method="PCA",
                     alignment_group="batch",
                     residual_model_formula_str="~educ+pmi")
cds_mic = reduce_dimension(cds_mic)
cds_mic = cluster_cells(cds_mic, cluster_method="louvain", k=500)
plot_cells(cds_mic, color_cells_by="cluster",cell_size=0.5,label_cell_groups=1,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=1,show_trajectory_graph=FALSE)#+theme(
plot_cells(cds_mic, color_cells_by="ros_ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#Choose k=500, which yields 17 excitatory neuorn clusters
cds_mic$cluster <- clusters(cds_mic)
table(cds_mic$Subcluster2)
cds_mic$Subcluster2[cds_mic$cluster==1] <- 'Ex11'
cds_mic$Subcluster2[cds_mic$cluster==2] <- 'Ex3'
cds_mic$Subcluster2[cds_mic$cluster==3] <- 'Ex7'
cds_mic$Subcluster2[cds_mic$cluster==4] <- 'Ex9'
cds_mic$Subcluster2[cds_mic$cluster==5] <- 'Ex2'
cds_mic$Subcluster2[cds_mic$cluster==6] <- 'Ex8'
cds_mic$Subcluster2[cds_mic$cluster==7] <- 'Ex1'
cds_mic$Subcluster2[cds_mic$cluster==8] <- 'Ex4'
cds_mic$Subcluster2[cds_mic$cluster==9] <- 'Ex10'
cds_mic$Subcluster2[cds_mic$cluster==10] <- 'Ex12'
cds_mic$Subcluster2[cds_mic$cluster==11] <- 'Ex5'
cds_mic$Subcluster2[cds_mic$cluster==12] <- 'Ex13'
cds_mic$Subcluster2[cds_mic$cluster==13] <- 'Ex14'
cds_mic$Subcluster2[cds_mic$cluster==14] <- 'Ex15'
cds_mic$Subcluster2[cds_mic$cluster==15] <- 'Ex0'
cds_mic$Subcluster2[cds_mic$cluster==16] <- 'Ex6'
cds_mic$Subcluster2[cds_mic$cluster==17] <- 'Ex16'

#extract the new subcluster labels
Ex_subs <- data.frame(
  TAG = colData(cds_mic)$TAG,
  Subcluster2 = colData(cds_mic)$Subcluster2
)
write.csv(Ex_subs, file="~/scAD_analysis/Ex_subs.csv")

plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=0.4,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))


################ find top genes within partitions in each broad cell type to classify subclusters #####################
#Inhibitory neurons

cds_mic = cds[,cds$broad.cell.type2=='In']
plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="sex",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#cds_mic <- clear_cds_slots(cds_mic)
#cds_mic = preprocess_cds(cds_mic, num_dim = 30,method="PCA", norm_method="none")
#adjust for batch effects
cds_mic$batch <- as.character(cds_mic$batch)
cds_mic$batch <- as.factor(cds_mic$batch)
cds_mic <- align_cds(cds_mic, 
                     preprocess_method="PCA",
                     alignment_group="batch",
                     residual_model_formula_str="~educ+pmi")
cds_mic = reduce_dimension(cds_mic)
cds_mic = cluster_cells(cds_mic, cluster_method="louvain", k=300)
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster2",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=1,show_trajectory_graph=FALSE)#+theme(
plot_cells(cds_mic, color_cells_by="ros_ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#Choose k=300, which yields 11 inhibitory neuron clusters
cds_mic$cluster <- clusters(cds_mic)
table(cds_mic$cluster)
cds_mic$Subcluster2[cds_mic$cluster==1] <- 'In3'
cds_mic$Subcluster2[cds_mic$cluster==2] <- 'In5'
cds_mic$Subcluster2[cds_mic$cluster==3] <- 'In0'
cds_mic$Subcluster2[cds_mic$cluster==4] <- 'In1'
cds_mic$Subcluster2[cds_mic$cluster==5] <- 'In7'
cds_mic$Subcluster2[cds_mic$cluster==6] <- 'In4'
cds_mic$Subcluster2[cds_mic$cluster==7] <- 'In8'
cds_mic$Subcluster2[cds_mic$cluster==8] <- 'In2'
cds_mic$Subcluster2[cds_mic$cluster==9] <- 'In9'
cds_mic$Subcluster2[cds_mic$cluster==10] <- 'In10'
cds_mic$Subcluster2[cds_mic$cluster==11] <- 'In6'


#extract the new subcluster labels
In_subs <- data.frame(
  TAG = colData(cds_mic)$TAG,
  Subcluster2 = colData(cds_mic)$Subcluster2
)
write.csv(In_subs, file="~/scAD_analysis/In_subs.csv")

plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))


################ find top genes within partitions in each broad cell type to classify subclusters #####################
#Astrocytes

cds_mic = cds[,cds$broad.cell.type2=='Ast']
plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#cds_mic <- clear_cds_slots(cds_mic)
#cds_mic = preprocess_cds(cds_mic, num_dim = 30,method="PCA", norm_method="none")
#adjust for batch effects
cds_mic$batch <- as.character(cds_mic$batch)
cds_mic$batch <- as.factor(cds_mic$batch)
cds_mic <- align_cds(cds_mic, 
                     preprocess_method="PCA",
                     alignment_group="batch",
                     residual_model_formula_str="~educ+pmi")
cds_mic = reduce_dimension(cds_mic)
cds_mic = cluster_cells(cds_mic, cluster_method="louvain", k=500)
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=1,show_trajectory_graph=FALSE)#+theme(
plot_cells(cds_mic, color_cells_by="ros_ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#Choose k=500, which yields 5 Astrocyte clusters
cds_mic$cluster <- clusters(cds_mic)
table(cds_mic$cluster)
cds_mic$Subcluster2[cds_mic$cluster==1] <- 'Ast1'
cds_mic$Subcluster2[cds_mic$cluster==2] <- 'Ast0'
cds_mic$Subcluster2[cds_mic$cluster==3] <- 'Ast2'
cds_mic$Subcluster2[cds_mic$cluster==4] <- 'Ast3'
cds_mic$Subcluster2[cds_mic$cluster==5] <- 'Ast4'



#extract the new subcluster labels
Ast_subs <- data.frame(
  TAG = colData(cds_mic)$TAG,
  Subcluster2 = colData(cds_mic)$Subcluster2
)
write.csv(Ast_subs, file="~/scAD_analysis/Ast_subs.csv")

plot_cells(cds_mic, color_cells_by="Subcluster2",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))




################ find top genes within partitions in each broad cell type to classify subclusters #####################
#Oligodendrocytes

cds_mic = cds[,cds$broad.cell.type2=='Oli']
plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#cds_mic <- clear_cds_slots(cds_mic)
#cds_mic = preprocess_cds(cds_mic, num_dim = 30,method="PCA", norm_method="none")
#adjust for batch effects
cds_mic$batch <- as.character(cds_mic$batch)
cds_mic$batch <- as.factor(cds_mic$batch)
cds_mic <- align_cds(cds_mic, 
                     preprocess_method="PCA",
                     alignment_group="batch",
                     residual_model_formula_str="~educ+pmi")
cds_mic = reduce_dimension(cds_mic)
cds_mic = cluster_cells(cds_mic, cluster_method="louvain", k=800)
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=1,show_trajectory_graph=FALSE)#+theme(
plot_cells(cds_mic, color_cells_by="ros_ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#Choose k=800, which yields 8 Oligo clusters
cds_mic$cluster <- clusters(cds_mic)
table(cds_mic$cluster)
cds_mic$Subcluster2[cds_mic$cluster==1] <- 'Oli0'
cds_mic$Subcluster2[cds_mic$cluster==2] <- 'Oli1'
cds_mic$Subcluster2[cds_mic$cluster==3] <- 'Oli2'
cds_mic$Subcluster2[cds_mic$cluster==4] <- 'Oli4'
cds_mic$Subcluster2[cds_mic$cluster==5] <- 'Oli3'
cds_mic$Subcluster2[cds_mic$cluster==6] <- 'Oli5'
cds_mic$Subcluster2[cds_mic$cluster==7] <- 'Oli6'
cds_mic$Subcluster2[cds_mic$cluster==8] <- 'Oli7'



#extract the new subcluster labels
Oli_subs <- data.frame(
  TAG = colData(cds_mic)$TAG,
  Subcluster2 = colData(cds_mic)$Subcluster2
)
write.csv(Oli_subs, file="~/scAD_analysis/Oli_subs.csv")

plot_cells(cds_mic, color_cells_by="Subcluster2",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))




################ find top genes within partitions in each broad cell type to classify subclusters #####################
#Opcs

cds_mic = cds[,cds$broad.cell.type2=='Opc']
plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#cds_mic <- clear_cds_slots(cds_mic)
#cds_mic = preprocess_cds(cds_mic, num_dim = 30,method="PCA", norm_method="none")
#adjust for batch effects
cds_mic$batch <- as.character(cds_mic$batch)
cds_mic$batch <- as.factor(cds_mic$batch)
cds_mic <- align_cds(cds_mic, 
                     preprocess_method="PCA",
                     alignment_group="batch",
                     residual_model_formula_str="~educ+pmi")
cds_mic = reduce_dimension(cds_mic)
cds_mic = cluster_cells(cds_mic, cluster_method="louvain", k=500)
plot_cells(cds_mic, color_cells_by="cluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=1,show_trajectory_graph=FALSE)#+theme(
plot_cells(cds_mic, color_cells_by="ros_ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#Choose k=500, which yields 4 Opc clusters
cds_mic$cluster <- clusters(cds_mic)
table(cds_mic$cluster)
cds_mic$Subcluster2[cds_mic$cluster==1] <- 'Opc0'
cds_mic$Subcluster2[cds_mic$cluster==2] <- 'Opc1'
cds_mic$Subcluster2[cds_mic$cluster==3] <- 'Opc2'
cds_mic$Subcluster2[cds_mic$cluster==4] <- 'Opc3'

table(cds_mic$Subcluster2)


#extract the new subcluster labels
Opc_subs <- data.frame(
  TAG = colData(cds_mic)$TAG,
  Subcluster2 = colData(cds_mic)$Subcluster2
)
write.csv(Opc_subs, file="~/scAD_analysis/Opc_subs.csv")

plot_cells(cds_mic, color_cells_by="Subcluster2",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))


################ find top genes within partitions in each broad cell type to classify subclusters #####################
#Endothelial cells and pericytes

cds_mic = cds[,cds$broad.cell.type2=='End']

cds_mic$batch <- as.character(cds_mic$batch)
cds_mic$batch <- as.factor(cds_mic$batch)
cds_mic <- align_cds(cds_mic, 
                     preprocess_method="PCA",
                     alignment_group="batch",
                     residual_model_formula_str="~educ+pmi")
cds_mic = reduce_dimension(cds_mic)
cds_mic = cluster_cells(cds_mic, cluster_method="louvain")
plot_cells(cds_mic, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

cds_mic$partition <- partitions(cds_mic)
table(cds_mic$partition)
cds_mic$Subcluster2[cds_mic$partition==1] <- 'End1'
cds_mic$Subcluster2[cds_mic$partition==2] <- 'End2'
table(cds_mic$Subcluster2)


#extract the new subcluster labels
End_subs <- data.frame(
  TAG = colData(cds_mic)$TAG,
  Subcluster2 = colData(cds_mic)$Subcluster2
)
write.csv(End_subs, file="~/scAD_analysis/End_subs.csv")


cds_mic = cds[,cds$broad.cell.type2=='Per']

cds_mic$batch <- as.character(cds_mic$batch)
cds_mic$batch <- as.factor(cds_mic$batch)
cds_mic <- align_cds(cds_mic, 
                     preprocess_method="PCA",
                     alignment_group="batch",
                     residual_model_formula_str="~educ+pmi")
cds_mic = reduce_dimension(cds_mic)
cds_mic = cluster_cells(cds_mic, cluster_method="louvain")
plot_cells(cds_mic, color_cells_by="partition",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

cds_mic$Subcluster2 <- 'Per'
table(cds_mic$Subcluster2)

#extract the new subcluster labels
Per_subs <- data.frame(
  TAG = colData(cds_mic)$TAG,
  Subcluster2 = colData(cds_mic)$Subcluster2
)
write.csv(Per_subs, file="~/scAD_analysis/Per_subs.csv")

cds_unknowns = cds[,cds$broad.cell.type2=='unknown']
Unknown_subs <- data.frame(
  TAG = colData(cds_unknowns)$TAG,
  Subcluster2 = colData(cds_unknowns)$broad.cell.type2
)
write.csv(Unknown_subs, file="~/scAD_analysis/Unknown_subs.csv")

# need to append the original cds with the new subcluster assignments
Subcluster2 <- rbind(Ast_subs, Mic_subs, In_subs, Ex_subs, Oli_subs, Opc_subs, End_subs, Per_subs, Unknown_subs)


rownames(Subcluster2) = Subcluster2$TAG
#reorder Subcluster2 df to match cds
tagnames_ordered <- match(colnames(cds), rownames(Subcluster2))
tagnames_ordered 
Subcluster2 <- Subcluster2[tagnames_ordered,]


cds$Subcluster2<-Subcluster2$Subcluster2


head(colData(cds))
subclusters <- data.frame(
  Subcluster = colData(cds)$Subcluster,
  Subcluster2 = colData(cds)$Subcluster2,
  broad.cell.type2 = colData(cds)$broad.cell.type2,
  ros_ids = colData(cds)$ros_ids,
  TAG = colData(cds)$TAG
)


saveRDS(cds, file="~/scAD_analysis/Mathys_scrannorm_newsubclusters_cds.rds")
#keep a spare full cds here:
cds_all <- cds




##plot diagnosis differences over pseudotime by sex:
#make sex-specific cds and clear the cds slots and reprocess the data
cdsF = cds[,cds$sex=='female']
dim(cdsF)
cdsF <- clear_cds_slots(cdsF)
cdsF$batch <- as.character(cdsF$batch)
cdsF$batch <- as.factor(cdsF$batch)
cdsF = preprocess_cds(cdsF, num_dim = 30,method="PCA",norm_method="none")
cdsF <- align_cds(cdsF, 
                  preprocess_method="PCA",
                  alignment_group="batch",
                  residual_model_formula_str="~educ+pmi")
cdsF = reduce_dimension(cdsF)
plot_cells(cdsF, color_cells_by="Subcluster2")
cdsF$Diagnosis[cdsF$Diagnosis=='Control']='Cont'
cdsF$Diagnosis[cdsF$Diagnosis=='Late-path']='Late'
cdsF$Diagnosis[cdsF$Diagnosis=='Early-path']='Early'

#Male cds 
cdsM = cds[,cds$sex=='male']
dim(cdsM)
cdsM <- clear_cds_slots(cdsM)
cdsM$batch <- as.character(cdsM$batch)
cdsM$batch <- as.factor(cdsM$batch)
cdsM = preprocess_cds(cdsM, num_dim = 30,method="PCA",norm_method="none")
cdsM <- align_cds(cdsM, 
                  preprocess_method="PCA",
                  alignment_group="batch",
                  residual_model_formula_str="~educ+pmi")
cdsM = reduce_dimension(cdsM)
plot_cells(cdsM, color_cells_by="Subcluster2")
cdsM$Diagnosis[cdsM$Diagnosis=='Control']='Cont'
cdsM$Diagnosis[cdsM$Diagnosis=='Late-path']='Late'
cdsM$Diagnosis[cdsM$Diagnosis=='Early-path']='Early'




# ####### Plot all trajectories for each new subcluster ####################
# #Calculate association between diagnosis (0/1) and pseudotime for each subcluster with >100 cells and plot, by sex:
# cds <- cdsF
# #cds <- cdsM
# plot_list<- list()
# 
# index=0
# for (subcluster in sort(unique(cds$Subcluster2))){
#   if (length(cds$Subcluster2[cds$Subcluster2==subcluster])>=100){
#     index=index+1
#     if (length(cds$Subcluster2[cds$Subcluster2==subcluster])<1000){
#       cell_size=1}
#     else if (length(cds$Subcluster2[cds$Subcluster2==subcluster])<3000){
#       cell_size=.5}
#     else {
#       cell_size=.5
#     }
#     print(subcluster)
#     cds_subset <- cds
#     cds_subset = cds_subset[,cds_subset$Subcluster2==subcluster]
#     cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
#     cds_subset$batch <- as.character(cds_subset$batch)
#     cds_subset$batch <- as.factor(cds_subset$batch)
#     cds_subset <- align_cds(cds_subset, 
#                             preprocess_method="PCA",
#                             alignment_group="batch",
#                             residual_model_formula_str="~educ+pmi")
#     cds_subset = reduce_dimension(cds_subset)
#     cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
#     
#     if (length(unique(cds_subset@clusters$UMAP$partitions))>1){
#       m<-0
#       p<-"1"
#       for (partition in unique(cds_subset@clusters$UMAP$partitions)){
#         if (length(cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions==partition])>m){
#           m<-length(cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions==partition])
#           p<-partition
#         }
#       }
#       cds_subset = cds_subset[,cds_subset@clusters$UMAP$partitions==p]
#       cds_subset = reduce_dimension(cds_subset)
#       cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
#     }
#     
#     pval_data <- data.frame(
#       batch = colData(cds_subset)$batch,
#       pmi = colData(cds_subset)$pmi,
#       educ = colData(cds_subset)$educ,
#       Diagnosis = cds_subset$pmi,
#       ros_ids = cds_subset$ros_ids
#     )
#     pval_data$Diagnosis[cds_subset$Diagnosis=='Cont']=0
#     pval_data$Diagnosis[cds_subset$Diagnosis=='Early']=1
#     pval_data$Diagnosis[cds_subset$Diagnosis=='Late']=1   
#     
#     cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions == "2"] <- "1"
#     cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions == "3"] <- "1"
#     cds_subset <- learn_graph(cds_subset)
#     inds = which(cds_subset$Diagnosis=='Cont')
#     inds = which(cds_subset$cogdx[inds]<2)
#     inds <- sample(inds) # shuffle indices
#     
#     pvals = c()
#     logfcs = c()
#     best_ind = 1
#     best_ind2 = 1
#     best_dif = -Inf
#     i=0
#     for (ind in inds[1:min(length(inds),50)]){
#       i=i+1
#       cds_subset = order_cells(cds_subset,root_cells=row.names(colData(cds_subset))[ind])
#       mc = median(pseudotime(cds_subset)[cds_subset$Diagnosis=='Cont'])
#       ma = median(pseudotime(cds_subset)[cds_subset$Diagnosis!='Cont'])
#       if (mc > ma){
#         pval = 1
#         logfc = 0
#       }
#       else{
#         cds_subset$pseudotime = pseudotime(cds_subset)
#         pval_data$pseudotime = cds_subset$pseudotime
#         pval_data$Diagnosis = as.numeric(pval_data$Diagnosis)
#         dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
#         fit1 <- lmFit(pval_data$Diagnosis,dm)
#         fit2 <- eBayes(fit1)
#         pval=topTable(fit2,coef=2)$adj.P.Val   
#         logfc = topTable(fit2,coef=2)$logFC
#       }
#       if (ma-mc > best_dif){
#         best_dif = ma-mc
#         best_ind = i
#         best_ind2 = ind
#       }
#       pvals <- c(pvals,pval)
#       logfcs <- c(logfcs,logfc)
#     }
#     pval = pvals[best_ind]
#     logfc = mean(logfcs)
#     print(logfc)
#     print(pval)
#     cds_subset = order_cells(cds_subset,root_cells=row.names(colData(cds_subset))[best_ind2]) # set root to best ind    
#     cds_subset$pseudotime = pseudotime(cds_subset)
#     pval_data$pseudotime = cds_subset$pseudotime
#     stars=''
#     if (pval<.01){
#       if (pval < .01){
#         stars = '*'
#       }
#       #         if (logfc > .02){
#       #             stars = '**'
#       #         }
#       #         if (logfc > .03){
#       #             stars = '***'
#       #}
#     }
#     pval_data$Diagnosis = as.factor(cds_subset$Diagnosis)
#     dp <- ggplot(pval_data, aes(x=Diagnosis, y=pseudotime, fill=Diagnosis)) + 
#       geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
#     p1<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#       theme(axis.title=element_text(size=8))+
#       ggtitle(paste0('p=',formatC(pval, format = "e", digits = 2)))+
#       theme(plot.title = element_text(size = 10))
#     
#     p2<-plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=cell_size,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)+
#       theme(legend.position = "none")+
#       theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#       theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
#       ggtitle(paste0('logFC=',formatC(logfc, format = "f", digits = 3),stars))+
#       theme(plot.title = element_text(size = 10))
#     p3<-arrangeGrob(p1,p2,ncol=2,top = textGrob(paste0(subcluster,stars),gp=gpar(fontsize=10,font=2)))
#     plot_list[[index]]<-p3
#   }}
# 
# #plot_list[sapply(plot_list, is.null)] <- NULL
# 
# 
# pdf(paste0("~/scAD_analysis/Figures/subclusters_pstime_M_part1",".pdf"))
# do.call(grid.arrange, c(plot_list[1:18], ncol=3))
# dev.off()
# 
# pdf(paste0("~/scAD_analysis/Figures/subclusters_pstime_M_part2",".pdf"))
# do.call(grid.arrange, c(plot_list[19:36], ncol=3))
# dev.off()
# 
# pdf(paste0("~/scAD_analysis/Figures/subclusters_pstime_M_part3",".pdf"))
# do.call(grid.arrange, c(plot_list[37:length(plot_list)], ncol=3))
# dev.off()








#Calculate association between diagnosis (0/1) and pseudotime for each subcluster with >100 cells and plot, by sex:
cds <- cdsF
#cds <- cdsM
plot_list<- list()
index=0
for (subcluster in sort(unique(cds$Subcluster2))){
  if (length(cds$Subcluster2[cds$Subcluster2==subcluster])>=100){
    index=index+1
    if (length(cds$Subcluster2[cds$Subcluster2==subcluster])<1000){
      cell_size=1}
    else if (length(cds$Subcluster2[cds$Subcluster2==subcluster])<3000){
      cell_size=.5}
    else {
      cell_size=.5
    }
    print(subcluster)
    cds_subset <- cds
    cds_subset = cds_subset[,cds_subset$Subcluster2==subcluster]
    cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
    cds_subset$batch <- as.character(cds_subset$batch)
    cds_subset$batch <- as.factor(cds_subset$batch)
    cds_subset <- align_cds(cds_subset,
                            preprocess_method="PCA",
                            alignment_group="batch",
                            residual_model_formula_str="~educ+pmi")
    cds_subset = reduce_dimension(cds_subset)
    cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
    cds_subset <- learn_graph(cds_subset)
    
    cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions == "2"] <- "1"
    cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions == "3"] <- "1"
    cds_subset <- learn_graph(cds_subset)

    get_earliest_principal_node <- function(cds_subset, diagnosis="Cont"){
      cell_ids <- which(colData(cds_subset)[, "Diagnosis"]==diagnosis)

      closest_vertex <-
        cds_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds_subset), ])
      root_pr_nodes <-
        igraph::V(principal_graph(cds_subset)[["UMAP"]])$name[as.numeric(names
                                                                         (which.max(table(closest_vertex[cell_ids,]))))]
      }
      cds_subset <- order_cells(cds_subset, root_pr_nodes=get_earliest_principal_node(cds_subset))
      cds_subset$pseudotime = pseudotime(cds_subset)


    pval_data <- data.frame(
      batch = colData(cds_subset)$batch,
      pmi = colData(cds_subset)$pmi,
      educ = colData(cds_subset)$educ,
      Diagnosis = cds_subset$pmi,
      ros_ids = cds_subset$ros_ids,
      pseudotime = colData(cds_subset)$pseudotime
    )
    pval_data$Diagnosis[cds_subset$Diagnosis=='Cont']=0
    pval_data$Diagnosis[cds_subset$Diagnosis=='Early']=1
    pval_data$Diagnosis[cds_subset$Diagnosis=='Late']=1

    dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
    fit1 <- lmFit(pval_data$Diagnosis,dm)
    fit2 <- eBayes(fit1)
    pval=topTable(fit2,coef=2)$adj.P.Val
    logfc = topTable(fit2,coef=2)$logFC
    stars=''
    if (pval<.01){
      if (pval < .01){
        stars = '*'
      }
      #         if (logfc > .02){
      #             stars = '**'
      #         }
      #         if (logfc > .03){
      #             stars = '***'
      #}
    }
  print(pval)
  print(logfc)

    pval_data$Diagnosis = as.factor(cds_subset$Diagnosis)
    dp <- ggplot(pval_data, aes(x=Diagnosis, y=pseudotime, fill=Diagnosis)) +
      geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
    p1<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
      theme(axis.title=element_text(size=8))+
      ggtitle(paste0('p=',formatC(pval, format = "e", digits = 2)))+
      theme(plot.title = element_text(size = 10))

    p2<-plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=cell_size,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)+
      theme(legend.position = "none")+
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      ggtitle(paste0('logFC=',formatC(logfc, format = "f", digits = 3),stars))+
      theme(plot.title = element_text(size = 10))
    p3<-arrangeGrob(p1,p2,ncol=2,top = textGrob(paste0(subcluster,stars),gp=gpar(fontsize=10,font=2)))
    plot_list[[index]]<-p3
  }}

pdf(paste0("~/scAD_analysis/Figures/subclusters_pstime_F_part1A",".pdf"))
do.call(grid.arrange, c(plot_list[1:18], ncol=3))
dev.off()

pdf(paste0("~/scAD_analysis/Figures/subclusters_pstime_F_part2A",".pdf"))
do.call(grid.arrange, c(plot_list[19:36], ncol=3))
dev.off()

pdf(paste0("~/scAD_analysis/Figures/subclusters_pstime_F_part3A",".pdf"))
do.call(grid.arrange, c(plot_list[37:length(plot_list)], ncol=3))
dev.off()













