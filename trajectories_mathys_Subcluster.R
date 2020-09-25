### Load packages ###
library(Matrix)

synapser::synLogin()


#get counts matrix, batch info, and cell tags with mathys broad cell types and subcluster info
synapser::synLogin()

p <- synapser::synGet('syn18686381')
counts <- Matrix::readMM(p$path)
syn18642934 <- synapser::synGet(entity='syn18642934') 
batches <- read.csv(syn18642934$path,stringsAsFactors = F)
p2 <- synapser::synGet('syn18686372')
Labels <- read.delim(p2$path,stringsAsFactors = F)
p3 <- synapser::synGet('syn18686382')
rownames(counts) <- readLines(p3$path) 



# get metadata files and rosmap IDs to match sample/patient/study IDs to each other
p3 <- synapser::synGet('syn18694015')
ids <- read.csv(p3$path,stringsAsFactors = F)
metadata <- read.csv('~/scAD_analysis/Mathys_supplement_5.csv')
p4 <- synapser::synGet(entity='syn3191087') 
metadata2 <- read.csv(p4$path)
metadata2 <- plyr::match_df(metadata2, ids, on="projid")
paste('Imputing PMI to:',median(metadata2$pmi[!is.na(metadata2$pmi)]))


metadata2$pmi[is.na(metadata2$pmi)] <- 7
metadata2$apoe_genotype[is.na(metadata2$apoe_genotype)] = 33
head(metadata2)

#Create columns to add to the cell metadata dataframe, and map rosmap ids (ros_ids) to sample ids (fastq) and project ids (projid). We will ultimately add the metadata and batch information to the Labels file so that each cell is appended with appropriate metadata.

sex = c()
m = as.character(metadata$msex)
fileName = c()
batch = c()
ros_ids = c()
projid = c()
Diagnosis = c()
cogdx = c()
ceradsc = c()
braaksc = c()
tangles = c()
apoe_genotype = c()
pmi = c()
educ=c()
race=c()
AOD=c()

#Need to be able to harmonize labels/counts/metadata by the different identifiers (projid, Subject, rosids)
for(i in 1:length(rownames(Labels))){
  ros_ids = c(ros_ids,ids$Subject[c(which(ids$projid==Labels$projid[i])[1])])
  fileName = c(fileName,ids$fastq[c(which(ids$projid==Labels$projid[i])[1])])
}
Labels$ros_ids = ros_ids
Labels$fileName = fileName
head(Labels)

#map metadata to cells
for(i in 1:length(rownames(Labels))){
  batch = c(batch,batches$sequencingBatch[c(which(batches$fileName==Labels$fileName[i])[1])])
  Diagnosis = c(Diagnosis,metadata$pathology.group[c(which(metadata$Subject==Labels$ros_id[i])[1])])
  cogdx = c(cogdx,metadata$cogdx[c(which(metadata$Subject==Labels$ros_id[i])[1])])
  sex = c(sex,m[c(which(metadata$Subject==Labels$ros_id[i])[1])])
  ceradsc = c(ceradsc,metadata$ceradsc[c(which(metadata$Subject==Labels$ros_id[i])[1])])
  braaksc = c(braaksc,metadata$braaksc[c(which(metadata$Subject==Labels$ros_id[i])[1])])
  tangles = c(tangles,metadata$tangles[c(which(metadata$Subject==Labels$ros_id[i])[1])])
  apoe_genotype = c(apoe_genotype,metadata2$apoe_genotype[c(which(metadata2$projid==Labels$projid[i])[1])])
  pmi = c(pmi,metadata2$pmi[c(which(metadata2$projid==Labels$projid[i])[1])])
  educ = c(educ,metadata2$educ[c(which(metadata2$projid==Labels$projid[i])[1])])
  race = c(race,metadata2$race[c(which(metadata2$projid==Labels$projid[i])[1])])
  AOD = c(AOD,metadata2$age_first_ad_dx[c(which(metadata2$projid==Labels$projid[i])[1])])
}


Labels$batch = batch
Labels$Diagnosis = Diagnosis
Labels$cogdx = cogdx
Labels$sex = sex
Labels$ceradsc = ceradsc
Labels$braaksc = braaksc
Labels$tangles = tangles
Labels$apoe_genotype = apoe_genotype
Labels$pmi = pmi
Labels$educ=educ
Labels$race=race
Labels$AOD=AOD

Labels$Diagnosis[Labels$Diagnosis=='late-pathology'] <- "Late-path"
Labels$Diagnosis[Labels$Diagnosis=='no-pathology'] <- 'Control'
Labels$Diagnosis[Labels$Diagnosis=='early-pathology'] <- 'Early-path'
Labels$simpleDiagnosis = Labels$Diagnosis
Labels$simpleDiagnosis[Labels$simpleDiagnosis!='Control'] <- "AD"
Labels$Diagnosis[Labels$Diagnosis==1]='Early-path'
Labels$Diagnosis[Labels$Diagnosis==2]='Late-path'
Labels$Diagnosis[Labels$Diagnosis==3]='Control'
head(Labels)



####Preprocessing#####

#detach("package:synapser", unload=TRUE)
#unloadNamespace("PythonEmbedInR") # must unload synapser because causes multiple definitions of S4Vectors

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
colnames(counts) <- Labels[,1]
rownames(Labels) <- Labels[,1]
gene_short_name <- data.frame(rownames(counts))
rownames(gene_short_name) <- rownames(counts)
gene_short_name$gene_short_name <- rownames(counts)
head(gene_short_name)


#Save the components of the cds so you don't need to recreate them again later
saveRDS(counts, file="~/scAD_analysis/mathys_counts_matrix.rds")
saveRDS(Labels, file="~/scAD_analysis/mathys_Labels.rds")
saveRDS(gene_short_name, file="~/scAD_analysis/mathys_gene_short_name.rds")
#counts <- readRDS(file="~/scAD_analysis/mathys_counts_matrix.rds")
#Labels <- readRDS(file="~/scAD_analysis/mathys_Labels.rds")
#gene_short_name <- readRDS(file="~/scAD_analysis/mathys_gene_short_name.rds")

cds <- new_cell_data_set(counts,
                         cell_metadata = Labels,
                         gene_metadata = gene_short_name)


#Clean up some variables before processing and plotting:

cds$educ = as.numeric(cds$educ)
cds$Education = cds$educ
cds$Sex = cds$sex
cds$CERAD = cds$ceradsc
cds$Braak = cds$braaksc
cds$APOE_genotype = cds$apoe_genotype
cds$batch = as.factor(cds$batch)

#Preprocessing the data using PCA
cds = preprocess_cds(cds, num_dim = 30,method="PCA",norm_method="size_only")
plot_pc_variance_explained(cds)

#adjust for batch effects
cds <- align_cds(cds, 
                 preprocess_method="PCA",
                 alignment_group="batch",
                 residual_model_formula_str="~educ+pmi")
cds = reduce_dimension(cds)
cds = cluster_cells(cds)
plot_cells(cds, color_cells_by="broad.cell.type",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

### save progress here to avoid rebuilding monocle object if R is terminated
saveRDS(cds, file = "Mathys_cds.rds")
#cds = readRDS(file = "~/scAD_analysis/Mathys_cds.rds")

plot_cells(cds, color_cells_by="Diagnosis",cell_size=.4,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds, color_cells_by="Sex",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cds, color_cells_by="Subcluster",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))



##plot diagnosis differences over pseudotime
library(qvalue)
library(metaseqR)

#save the entire cds so it doesn't get written over
cds_processed <- cds




##plot diagnosis differences over pseudotime by sex:
#make sex-specific cds and clear the cds slots and rescale data
cdsF = cds[,cds$sex=='female']
dim(cdsF)
cdsF <- clear_cds_slots(cdsF)
cdsF$batch <- as.character(cdsF$batch)
cdsF$batch <- as.factor(cdsF$batch)
cdsF = preprocess_cds(cdsF, num_dim = 30,method="PCA", norm_method = "size_only")
plot_pc_variance_explained(cdsF)
cdsF <- align_cds(cdsF, 
                        preprocess_method="PCA",
                        alignment_group="batch",
                        residual_model_formula_str="~educ+pmi")
cdsF = reduce_dimension(cdsF)
plot_cells(cdsF, color_cells_by="Subcluster")
cdsF$Diagnosis[cdsF$Diagnosis=='Control']='Cont'
cdsF$Diagnosis[cdsF$Diagnosis=='Late-path']='Late'
cdsF$Diagnosis[cdsF$Diagnosis=='Early-path']='Early'

#Male cds 
cdsM = cds[,cds$sex=='male']
dim(cdsM)
cdsM <- clear_cds_slots(cdsM)
cdsM$batch <- as.character(cdsM$batch)
cdsM$batch <- as.factor(cdsM$batch)
cdsM = preprocess_cds(cdsM, num_dim = 30,method="PCA", norm_method = "size_only")
plot_pc_variance_explained(cdsM)
cdsM <- align_cds(cdsM, 
                  preprocess_method="PCA",
                  alignment_group="batch",
                  residual_model_formula_str="~educ+pmi")
cdsM = reduce_dimension(cdsM)
plot_cells(cdsM, color_cells_by="Subcluster")
cdsM$Diagnosis[cdsM$Diagnosis=='Control']='Cont'
cdsM$Diagnosis[cdsM$Diagnosis=='Late-path']='Late'
cdsM$Diagnosis[cdsM$Diagnosis=='Early-path']='Early'




#Calculate association between diagnosis (0/1) and pseudotime for each subcluster with >100 cells and plot, by sex:
# cds <- cdsF
# # #cds <- cdsM
# plot_list<- list()
# index=0
# for (subcluster in sort(unique(cds$Subcluster))){#unique(cds$Subcluster)
#   if (length(cds$Subcluster[cds$Subcluster==subcluster])>=100){
#     index=index+1
#     if (length(cds$Subcluster[cds$Subcluster==subcluster])<1000){
#       cell_size=1}
#     else if (length(cds$Subcluster[cds$Subcluster==subcluster])<3000){
#       cell_size=.5}
#     else {
#       cell_size=.5
#     }
#     print(subcluster)
#     cds_subset <- cds
#     cds_subset = cds_subset[,cds_subset$Subcluster==subcluster]
#     cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
#     cds_subset$batch <- as.character(cds_subset$batch)
#     cds_subset$batch <- as.factor(cds_subset$batch)
#     cds_subset = preprocess_cds(cds_subset, num_dim = 30,method="PCA",norm_method="size_only")
#     cds_subset <- align_cds(cds_subset,
#                             preprocess_method="PCA",
#                             alignment_group="batch",
#                             residual_model_formula_str="~educ+pmi")
#     cds_subset = reduce_dimension(cds_subset)
#     cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
#     cds_subset <- learn_graph(cds_subset)
# 
# 
#     get_earliest_principal_node <- function(cds_subset, diagnosis="Cont"){
#       cell_ids <- which(colData(cds_subset)[, "Diagnosis"]==diagnosis)
# 
#       closest_vertex <-
#         cds_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#       closest_vertex <- as.matrix(closest_vertex[colnames(cds_subset), ])
#       root_pr_nodes <-
#         igraph::V(principal_graph(cds_subset)[["UMAP"]])$name[as.numeric(names
#                                                                          (which.max(table(closest_vertex[cell_ids,]))))]
# 
#       root_pr_nodes
#       }
#       cds_subset <- order_cells(cds_subset, root_pr_nodes=get_earliest_principal_node(cds_subset))
#       cds_subset$pseudotime = pseudotime(cds_subset)
# 
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
#     pval_data$pseudotime = cds_subset$pseudotime
# 
#     dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
#     fit1 <- lmFit(pval_data$Diagnosis,dm)
#     fit2 <- eBayes(fit1)
#     pval=topTable(fit2,coef=2)$adj.P.Val
#     logfc = topTable(fit2,coef=2)$logFC
#     stars=''
#     if (logfc>.02){
#       if (logfc > .02){
#         stars = '*'
#       }
#       #         if (logfc > .02){
#       #             stars = '**'
#       #         }
#       #         if (logfc > .03){
#       #             stars = '***'
#       #}
#     }
#   print(pval)
#   print(logfc)
#   
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
# # 
# # #plot_list[sapply(plot_list, is.null)] <- NULL
# # 
# pdf(paste0("~/scAD_analysis/Figures_LH/pstime_dx_part1",".pdf"))
# do.call(grid.arrange, c(plot_list[1:11], ncol=3))
# dev.off()
# 
# 
# #pdf(paste0("~/scAD_analysis/figures2/mathys_sn_diagnosis_LHpart2_F_batched",".pdf"))
# pdf(paste0("~/scAD_analysis/figures2/mathys_sn_diagnosis_LHpart2_M_batched",".pdf"))
# do.call(grid.arrange, c(plot_list[18:length(plot_list)], ncol=3))
# dev.off()


#Calculate association between diagnosis (0/1) and pseudotime for each subcluster with >100 cells and plot, by sex:
cds <- cdsF
#cds <- cdsM
plot_list<- list()

index=0
for (subcluster in sort(unique(cds$Subcluster))){#unique(cds$Subcluster)
  if (length(cds$Subcluster[cds$Subcluster==subcluster])>=100){
    index=index+1
    if (length(cds$Subcluster[cds$Subcluster==subcluster])<1000){
      cell_size=1}
    else if (length(cds$Subcluster[cds$Subcluster==subcluster])<3000){
      cell_size=.5}
    else {
      cell_size=.5
    }
    print(subcluster)
    cds_subset <- cds
    cds_subset = cds_subset[,cds_subset$Subcluster==subcluster]
    cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
    cds_subset$batch <- as.character(cds_subset$batch)
    cds_subset$batch <- as.factor(cds_subset$batch)
    cds_subset = preprocess_cds(cds_subset, num_dim = 30,method="PCA",norm_method="size_only")
    cds_subset <- align_cds(cds_subset, 
                            preprocess_method="PCA",
                            alignment_group="batch",
                            residual_model_formula_str="~educ+pmi")
    cds_subset = reduce_dimension(cds_subset)
    cds_subset = cluster_cells(cds_subset, cluster_method="louvain")

    if (length(unique(cds_subset@clusters$UMAP$partitions))>1){
      m<-0
      p<-"1"
      for (partition in unique(cds_subset@clusters$UMAP$partitions)){
        if (length(cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions==partition])>m){
          m<-length(cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions==partition])
          p<-partition
        }
      }
      cds_subset = cds_subset[,cds_subset@clusters$UMAP$partitions==p]
      cds_subset = reduce_dimension(cds_subset)
      cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
    }
    
    pval_data <- data.frame(
      batch = colData(cds_subset)$batch,
      pmi = colData(cds_subset)$pmi,
      educ = colData(cds_subset)$educ,
      Diagnosis = cds_subset$pmi,
      ros_ids = cds_subset$ros_ids
    )
    pval_data$Diagnosis[cds_subset$Diagnosis=='Cont']=0
    pval_data$Diagnosis[cds_subset$Diagnosis=='Early']=1
    pval_data$Diagnosis[cds_subset$Diagnosis=='Late']=1   
    
    cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions == "2"] <- "1"
    cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions == "3"] <- "1"
    cds_subset <- learn_graph(cds_subset)
    inds = which(cds_subset$Diagnosis=='Cont')
    inds = which(cds_subset$cogdx[inds]<2)
    inds <- sample(inds) # shuffle indices
    
    pvals = c()
    logfcs = c()
    best_ind = 1
    best_ind2 = 1
    best_dif = -Inf
    i=0
    for (ind in inds[1:min(length(inds),50)]){
      i=i+1
      cds_subset = order_cells(cds_subset,root_cells=row.names(colData(cds_subset))[ind])
      mc = median(pseudotime(cds_subset)[cds_subset$Diagnosis=='Cont'])
      ma = median(pseudotime(cds_subset)[cds_subset$Diagnosis!='Cont'])
      if (mc > ma){
        pval = 1
        logfc = 0
      }
      else{
        cds_subset$pseudotime = pseudotime(cds_subset)
        pval_data$pseudotime = cds_subset$pseudotime
        pval_data$Diagnosis = as.numeric(pval_data$Diagnosis)
        dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
        fit1 <- lmFit(pval_data$Diagnosis,dm)
        fit2 <- eBayes(fit1)
        pval=topTable(fit2,coef=2)$adj.P.Val   
        logfc = topTable(fit2,coef=2)$logFC
      }
      if (ma-mc > best_dif){
        best_dif = ma-mc
        best_ind = i
        best_ind2 = ind
      }
      pvals <- c(pvals,pval)
      logfcs <- c(logfcs,logfc)
    }
    pval = pvals[best_ind]
    logfc = mean(logfcs)
    print(logfc)
    print(pval)
    cds_subset = order_cells(cds_subset,root_cells=row.names(colData(cds_subset))[best_ind2]) # set root to best ind    
    cds_subset$pseudotime = pseudotime(cds_subset)
    pval_data$pseudotime = cds_subset$pseudotime
    stars=''
    if (logfc>.02){
      if (logfc > .02){
        stars = '*'
      }
      #         if (logfc > .02){
      #             stars = '**'
      #         }
      #         if (logfc > .03){
      #             stars = '***'
      #}
    }
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

#plot_list[sapply(plot_list, is.null)] <- NULL

#pdf(paste0("~/scAD_analysis/figures2/mathys_sn_diagnosis_LHpart1_F_batched",".pdf"))
pdf(paste0("~/scAD_analysis/Figures_LH/pstime_dx_F_part1_size",".pdf"))
do.call(grid.arrange, c(plot_list[1:18], ncol=3))
dev.off()


#pdf(paste0("~/scAD_analysis/figures2/mathys_sn_diagnosis_LHpart2_F_batched",".pdf"))
pdf(paste0("~/scAD_analysis/Figures_LH/pstime_dx_F_part2_size",".pdf"))
do.call(grid.arrange, c(plot_list[19:length(plot_list)], ncol=3))
dev.off()


#plot diagnostic criterial (CERAD, Braak, cogdx, apoe genotype) for subtypes which show variation in diagnosis across pseudotime
cds <- cdsF
#cds <- cdsM
plot_list<- list()
index=0
#for (subcluster in c('Oli0','Oli1', 'Oli4', 'Oli5')){#,'Oli0','Oli1'
for (subcluster in c('Ast0','Ex5', 'In0', 'Opc1')){#,'Oli0','Oli1'
  index=index+1
  if (length(cds$Subcluster[cds$Subcluster==subcluster])<1000){
    cell_size=1}
  else if (length(cds$Subcluster[cds$Subcluster==subcluster])<3000){
    cell_size=.5}
  else {
    cell_size=.5
  }
  print(subcluster)
  cds_subset = cds
  cds_subset = cds_subset[,cds_subset$broad.cell.type==gsub('[[:digit:]]+', '', subcluster)]
  cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
  cds_subset$batch <- as.character(cds_subset$batch)
  cds_subset$batch <- as.factor(cds_subset$batch)
  cds_subset = preprocess_cds(cds_subset, num_dim = 30,method="PCA",norm_method="size_only")
  cds_subset <- align_cds(cds_subset, 
                          preprocess_method="PCA",
                          alignment_group="batch",
                          residual_model_formula_str="~educ+pmi")
  cds_subset = cds_subset[,cds_subset$Subcluster==subcluster] # mic3 subcluster
  cds_subset = reduce_dimension(cds_subset)
  cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
  if (length(unique(cds_subset@clusters$UMAP$partitions))>1){
    m<-0
    p<-"1"
    for (partition in unique(cds_subset@clusters$UMAP$partitions)){
      if (length(cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions==partition])>m){
        m<-length(cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions==partition])
        p<-partition
      }
    }
    cds_subset = cds_subset[,cds_subset@clusters$UMAP$partitions==p]
    cds_subset = reduce_dimension(cds_subset)
    cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
  }
  cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions == "2"] <- "1"
  cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions == "3"] <- "1"
  cds_subset <- learn_graph(cds_subset)
  
  inds = which(cds_subset$Diagnosis=='Cont')
  inds = which(cds_subset$CERAD[inds]>3)
  best_ind = 1
  best_ind2 = 1
  best_dif = -Inf
  i=0
  for (ind in inds[1:min(length(inds),50)]){
    i=i+1
    cds_subset = order_cells(cds_subset,root_cells=row.names(colData(cds_subset))[ind])
    mc = median(pseudotime(cds_subset)[cds_subset$Diagnosis=='Cont'])
    ma = median(pseudotime(cds_subset)[cds_subset$Diagnosis!='Cont'])
    if (ma-mc > best_dif){
      best_dif = ma-mc
      best_ind = i
      best_ind2 = ind
    }
  }
  cds_subset = order_cells(cds_subset,root_cells=row.names(colData(cds_subset))[best_ind2]) # set root to best ind    
  cds_subset$pseudotime = pseudotime(cds_subset)
  pval_data <- data.frame(
    pseudotime = cds_subset$pseudotime, 
    batch = colData(cds_subset)$batch,
    pmi = colData(cds_subset)$pmi,
    educ = colData(cds_subset)$educ,
    CERAD = cds_subset$CERAD,
    Braak = cds_subset$Braak,
    COGDX = cds_subset$cogdx,
    APOE_genotype = cds_subset$pmi,
    ros_ids = cds_subset$ros_ids
  )
  pval_data$pseudotime = cds_subset$pseudotime
  pval_data$APOE_genotype[cds_subset$APOE_genotype=='23']=0
  pval_data$APOE_genotype[cds_subset$APOE_genotype=='33']=1
  pval_data$APOE_genotype[cds_subset$APOE_genotype=='34']=2
  pval_data$APOE_genotype[cds_subset$APOE_genotype=='44']=3
  pval_data$APOE_genotype = as.numeric(pval_data$APOE_genotype)
  ### CERAD ###
  dm <- model.matrix(~pseudotime+pmi+batch+educ, data=pval_data)
  fit1 <- lmFit(pval_data$CERAD,dm)
  fit2 <- eBayes(fit1)
  pval=topTable(fit2,coef=2)$adj.P.Val
  stars=''
  if (pval < .01){
    stars = '*'
  }
  if (pval < .001){
    stars = '**'
  }
  if (pval < .0001){
    stars = '***'
  }
  pval_data$CERAD = as.factor(cds_subset$CERAD)
  dp <- ggplot(pval_data, aes(x=CERAD, y=pseudotime, fill=CERAD)) + 
    geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
  p1<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    ggtitle(paste('p=',formatC(pval, format = "e", digits = 2),stars))+theme(axis.title=element_text(size=10))
  ### Braak ###
  dm <- model.matrix(~pseudotime+pmi+batch+educ, data=pval_data)
  fit1 <- lmFit(pval_data$Braak,dm)
  fit2 <- eBayes(fit1)
  pval=topTable(fit2,coef=2)$adj.P.Val
  stars=''
  if (pval < .01){
    stars = '*'
  }
  if (pval < .001){
    stars = '**'
  }
  if (pval < .0001){
    stars = '***'
  }
  pval_data$Braak = as.factor(cds_subset$Braak)
  dp <- ggplot(pval_data, aes(x=Braak, y=pseudotime, fill=Braak)) + 
    geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
  p2<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    ggtitle(paste('p=',formatC(pval, format = "e", digits = 2),stars))+theme(axis.title=element_text(size=10))
  ### COGDX ###
  dm <- model.matrix(~pseudotime+pmi+batch+educ, data=pval_data)
  fit1 <- lmFit(pval_data$COGDX,dm)
  fit2 <- eBayes(fit1)
  pval=topTable(fit2,coef=2)$adj.P.Val
  stars=''
  if (pval < .01){
    stars = '*'
  }
  if (pval < .001){
    stars = '**'
  }
  if (pval < .0001){
    stars = '***'
  }
  pval_data$COGDX = as.factor(cds_subset$cogdx)
  dp <- ggplot(pval_data, aes(x=COGDX, y=pseudotime, fill=COGDX)) + 
    geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
  p3<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    ggtitle(paste('p=',formatC(pval, format = "e", digits = 2),stars))+theme(axis.title=element_text(size=10))
  ### APOE_genotype ###
  dm <- model.matrix(~pseudotime+pmi+batch+educ, data=pval_data)
  fit1 <- lmFit(pval_data$APOE_genotype,dm)
  fit2 <- eBayes(fit1)
  pval=topTable(fit2,coef=2)$adj.P.Val
  stars=''
  if (pval < .01){
    stars = '*'
  }
  if (pval < .001){
    stars = '**'
  }
  if (pval < .0001){
    stars = '***'
  }
  pval_data$APOE_genotype = as.factor(cds_subset$APOE_genotype)
  dp <- ggplot(pval_data, aes(x=APOE_genotype, y=pseudotime, fill=APOE_genotype)) + 
    geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
  p4<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    ggtitle(paste('p=',formatC(pval, format = "e", digits = 2),stars))+theme(axis.title=element_text(size=10))    
  
  p3<-arrangeGrob(p1,p2,p3,p4,ncol=4,top = textGrob(subcluster,gp=gpar(fontsize=15,font=2)))
  plot_list[[index]]<-p3
}
#pdf(paste0("~/scAD_analysis/Figures_LH/neuropath_Olis_female",".pdf"))
pdf(paste0("~/scAD_analysis/Figures_LH/neuropath_Ast0_Ex5_Mic1_In0_female",".pdf"))
pdf(paste0("~/scAD_analysis/Figures_LH/neuropath_Opc1_female",".pdf"))
do.call(grid.arrange, c(plot_list, ncol=1))
dev.off()


#################################### TRAJECTORIES ######################

#Explore a specific cell subcluster of your choice.
#Previous studies and preliminary analyses indicate important differences in gene expression in different cell types by sex. We will examine trajectories in female samples only.

#After subsetting the cds to include female subcluster cells only:
#  1) repeat preprocess_cds to normalize by sequencing depth in the new smaller data set. 
#2) only keep genes non-zero in at least 20 cells
#3) change batch into a character and back to factor (to delete batch bins with 0 counts, it is a weird glitch)
#4) adjust for batch/pmi/education, reduce dimensions & cluster, this time using the louvain method of clustering (default is leiden, which generally produces less complex trajectories). 
#5) Call the learn_graph function to plot a trajectory, then plot by diagnosis (control, early-pathology AD, or late-pathology AD).
cds_ex = cds[,cds$broad.cell.type=='Ex']
cds_subset <- cds_ex[,cds_ex$sex=='female'] 
cds_subset <- cds_subset[,cds_subset$Subcluster=='Ex5']
#cds_subset <- clear_cds_slots(cds_subset)
cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,]
#cds_subset <- preprocess_cds(cds_subset, num_dim = 30, method = "PCA", norm_method="log")
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

plot_cells(cds_subset, color_cells_by="ros_ids",cell_size=2,label_cell_groups=0)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))


#The function below will pick the root at the "beginning" of the trajectory by first grouping the cells according to which trajectory graph node they are nearest to. It calcuates what fraction of the cells at each node come from Control samples. Then it picks the node that is most heavily occupied by Control cells and returns that as the root.

#After setting the root, call the order_cells function to calculate pseudotime based on order, and plot by pseudotime:
  

get_earliest_principal_node <- function(cds_subset, diagnosis="Cont"){
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
#Create a dataframe that will contain variables necessary to regress: diagnosis ~ pseudotime + educ + pmi
#Note: we are not adjusting for batch because the align_cds function already effectively subtracted out batch effects and we do not want to overcorrect. the align_cds residual_model_formula_str function only miniimally affects the clustering of the cells, so these variables should be included.


pval_data <- data.frame(
  pmi = colData(cds_subset)$pmi,
  educ = colData(cds_subset)$educ,
  Diagnosis = cds_subset$Diagnosis,
  ros_ids = cds_subset$ros_ids)
pval_data$diagnosis2 <- ifelse(pval_data$Diagnosis=='Cont',0,1)
pval_data$pseudotime = cds_subset$pseudotime
fit <- lm(diagnosis2~pseudotime+pmi+educ, data=pval_data)
summary(fit)


#Use limma to pull the pseudotime coefficient and p-value, and assemble two plots: violin plot of pseudotime by diagnosis category alongside the trajectory plot.


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
pdf(paste0("~/scAD_analysis/Figures_LH/F_Ast0_dx_pseudotime.pdf"))
grid.arrange(arrangeGrob(p1,p2,p4,ncol=3,top = textGrob(paste0('Ast0 Female',stars),gp=gpar(fontsize=18,font=2))))
dev.off()












table(cds$broad.cell.type)
table(cds$Subcluster)

##### Microglia trajectories ####
cds_subset <- cds[,cds$broad.cell.type=='Mic']
cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
dim(cds_subset)
cds_subset <- preprocess_cds(cds_subset, num_dim = 30,residual_model_formula_str="~batch+educ+pmi")
cds_subset = reduce_dimension(cds_subset, reduction_method="PCA")
cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
plot_cells(cds_subset, color_cells_by="ros_ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cds_subset, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cds_subset, color_cells_by="partition",cell_size=1,label_cell_groups=1,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)

#exclude partition 1 cells
cds_subset$partition <- partitions(cds_subset)
cds_subset = cds_subset[,cds_subset$partition==2]
dim(cds_subset)
plot_cells(cds_subset, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
#there are two outliers. track a trajectory then pull these out by their pseudotimes
cds_subset = reduce_dimension(cds_subset, reduction_method="PCA")
cds_subset = cluster_cells(cds_subset, cluster_method="louvain")

cds_subset <- learn_graph(cds_subset)
cds_subset <- order_cells(cds_subset)
cds_subset$pseudotimes <- pseudotime(cds_subset)
summary(cds_subset$pseudotimes)
plot_cells(cds_subset, color_cells_by="pseudotime",cell_size=1,label_cell_groups=1,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
cds_subset <- cds_subset[,cds_subset$pseudotimes>7]
dim(cds_subset)

cds_subset = reduce_dimension(cds_subset, reduction_method="PCA")
cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
plot_cells(cds_subset, color_cells_by="ros_ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cds_subset, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)

##### Microglia trajectories ####
cds_subset<-0
cds_subset <- cds[,cds$Subcluster=='Mic1']
cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
dim(cds_subset)
cds_subset <- preprocess_cds(cds_subset, num_dim = 30,norm_method="size_only", residual_model_formula_str="~batch+educ+pmi")
cds_subset = reduce_dimension(cds_subset, reduction_method="PCA")
cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
plot_cells(cds_subset, color_cells_by="ros_ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cds_subset, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cds_subset, color_cells_by="partition",cell_size=1,label_cell_groups=1,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)

cds_subset <- learn_graph(cds_subset)
cds_subset <- order_cells(cds_subset)
cds_subset$pseudotimes <- pseudotime(cds_subset)
summary(cds_subset$pseudotimes)
plot_cells(cds_subset, color_cells_by="pseudotime",cell_size=1,label_cell_groups=1,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
cds_subset <- cds_subset[,cds_subset$pseudotimes>2]
dim(cds_subset)

cds_subset <- preprocess_cds(cds_subset, num_dim = 30,norm_method="size_only", residual_model_formula_str="~batch+educ+pmi")
cds_subset$batch <- as.factor(cds_subset$batch)
cds_subset <- align_cds(cds_subset, 
                    preprocess_method="PCA",
                    alignment_group="batch",
                    residual_model_formula_str="~educ+pmi",
                    verbose=TRUE)
cds_subset = reduce_dimension(cds_subset)
cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
plot_cells(cds_subset, color_cells_by="ros_ids",cell_size=2,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
plot_cells(cds_subset, color_cells_by="Subcluster",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE)


