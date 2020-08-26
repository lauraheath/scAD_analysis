### Load packages ###
library(Matrix)
library(synapser)
#devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
library(monocle3)
library(Seurat)
#BiocManager::install("variancePartition")
synLogin()

# get data: this is the "filtered_count_matrix.mtx" file (snRNAseqPFC_BA10->Gene Expression(RNAseq)->Processed)
p <- synapser::synGet('syn18686381')
counts <- readMM(p$path)

# get sample QC/batch data: snRNAseqPFC_BA10_assay_scRNAseq_metadata.csv
syn18642934 <- synGet(entity='syn18642934') 
batches <- read.csv(syn18642934$path,stringsAsFactors = F)

#get tags and celltype data: filtered_column_metadata.txt
p2 <- synapser::synGet('syn18686372')
Labels <- read.delim(p2$path,stringsAsFactors = F)

#get short gene names list and make them into rownames on counts file: filtered_gene_row_names.txt
p3 <- synapser::synGet('syn18686382')
rownames(counts) <- readLines(p3$path)   



# get rosmap id mappings: snRNAseqPFC_BA10_id_mapping.csv
p3 <- synapser::synGet('syn18694015')
ids <- read.csv(p3$path,stringsAsFactors = F)

####upload a csv version of supplemtary table 5 from mathys et al. into R to get neuropath data
metadata <- read.csv('~/scAD_analysis/ROSMAP_metadata2.csv')
#ROSMAP_clinical.csv
p4 <- synGet(entity='syn3191087') 
metadata2 <- read.csv(p4$path)

#extract the mathys patients only and impute the missing PMI as the median PMI of the dataset, add that to 
#patient missing pmi in metadata file. Also impute missing apoe genotypes to 3/3
metadata3 <- metadata2
metadata3 <- plyr::match_df(metadata3, ids, on="projid")

paste('Imputing PMI to:',median(metadata3$pmi[!is.na(metadata3$pmi)]))
#add this back into the metadata file
metadata2$pmi[is.na(metadata2$pmi)] <- 7
#impute missing apoe genotypes as '33'
metadata2$apoe_genotype[is.na(metadata2$apoe_genotype)] = 33
head(metadata2)


#create list categories for the features data frame
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

Labels$Diagnosis[Labels$Diagnosis=='late-pathology'] <- "Late-pathology AD"
Labels$Diagnosis[Labels$Diagnosis=='no-pathology'] <- 'Control'
Labels$Diagnosis[Labels$Diagnosis=='early-pathology'] <- 'Early-pathology AD'
Labels$simpleDiagnosis = Labels$Diagnosis
Labels$simpleDiagnosis[Labels$simpleDiagnosis!='Control'] <- "AD"
Labels$Diagnosis[Labels$Diagnosis==1]='Early-pathology AD'
Labels$Diagnosis[Labels$Diagnosis==2]='Late-pathology AD'
Labels$Diagnosis[Labels$Diagnosis==3]='Control'
head(Labels)

colnames(counts) <- Labels[,1]
rownames(Labels) <- Labels[,1]

####Preprocessing#####

detach("package:synapser", unload=TRUE)
unloadNamespace("PythonEmbedInR") # must unload synapser because causes multiple definitions of S4Vectors

library(scater)
library(monocle3)

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

gene_short_name <- data.frame(rownames(counts))
rownames(gene_short_name) <- rownames(counts)
gene_short_name$gene_short_name <- rownames(counts)
head(gene_short_name)
length(gene_short_name)

# inds = rowSums(counts)>=250
# length(inds)
# counts = counts[inds,] # only keep genes non-zero in at least 250 cells

#Save the components of the cds so you don't need to recreate them again later
saveRDS(counts, file="~/scAD_analysis/mathys_counts_matrix.rds")
saveRDS(Labels, file="~/scAD_analysis/mathys_Labels.rds")
saveRDS(gene_short_name, file="~/scAD_analysis/mathys_gene_short_name.rds")
counts <- readRDS(file="~/scAD_analysis/mathys_counts_matrix.rds")
Labels <- readRDS(file="~/scAD_analysis/mathys_Labels.rds")
gene_short_name <- readRDS(file="~/scAD_analysis/mathys_gene_short_name.rds")

cds <- new_cell_data_set(counts,
                         cell_metadata = Labels,
                         gene_metadata = gene_short_name)





### preprocessing and reduce dimensionality
### preprocessing includes library size normalization and regressing out batch, education, and pmi
### note: data is not log-normalized

#upload into Rstudio a csv version of the table with gene names, adj.p, IFC, and subpopulations from supplementary table 8 in mathys et al.
mathy_marker_genes <- read.csv('~/scAD_analysis/mathys_marker_genes.csv', header=TRUE)
#mathy_marker_genes$adj.pvals = as.numeric(as.character(mathy_marker_genes$adj.pvals))

genes2<-c()
for (gene in unique(c(as.vector(mathy_marker_genes$gene.name),c("SYT1","SNAP25","GRIN1","GAD1","GAD2","SLC17A7","CAMK2A","NRGN","AQP4",
                                                                "GFAP","MBP","MOBP","PLP1","PDGFRA","VCAN","CD74","CSF1R","C3","FLT1","CLDN5")))){
  if (gene %in% rownames(cds)){
    genes2 <- c(genes2,which(rownames(cds)==gene))
  }
}
length(genes2)
cds_mathys <- cds[genes2,]

library(ggplot2)
library(gridExtra)
library(MASS)
library(limma)
library(edgeR)
library(grid)
cds_mathys$educ = as.numeric(cds_mathys$educ)
cds_mathys$Education = cds_mathys$educ
cds_mathys$Sex = cds_mathys$sex
cds_mathys$CERAD = cds_mathys$ceradsc
cds_mathys$Braak = cds_mathys$braaksc
cds_mathys$APOE_genotype = cds_mathys$apoe_genotype
cds_mathys$batch = as.factor(cds_mathys$batch)

cds_mathys$Diagnosis[cds_mathys$Diagnosis=='Early-pathology AD']='Early-path AD'
cds_mathys$Diagnosis[cds_mathys$Diagnosis=='Late-pathology AD']='Late-path AD'





#dimension reduction and normalization
cds_mathys = preprocess_cds(cds_mathys, num_dim = 30,method="PCA",norm_method="size_only")
plot_pc_variance_explained(cds_mathys)
cds_mathys <- align_cds(cds_mathys, 
                        preprocess_method="PCA",
                        alignment_group="batch",
                        residual_model_formula_str="~educ+pmi")
cds_mathys = reduce_dimension(cds_mathys)
plot_cells(cds_mathys, color_cells_by="Subcluster")


### save progress here to avoid rebuilding monocle object if R is terminated
saveRDS(cds_mathys, file = "Mathys_monocle_preprocessed_cds.rds")
saveRDS(cds, file = "Mathys_allgenes_monocle_cds.rds")

cds_mathys = readRDS(file = "~/scAD_analysis/Mathys_monocle_preprocessed_cds.rds")
cds = readRDS(file = "~/scAD_analysis/Mathys_allgenes_monocle_cds.rds")


p1<-plot_cells(cds_mathys, color_cells_by="Education",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10))
p2<-plot_cells(cds_mathys, color_cells_by="Sex",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10))
p3<-plot_cells(cds_mathys, color_cells_by="Diagnosis",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8))
p4<-plot_cells(cds_mathys, color_cells_by="ros_ids",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(legend.position = "none")
p6<-plot_cells(cds_mathys, color_cells_by="Subcluster",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm") )
#pdf("~/mathys_sn_summary.pdf")
grid.arrange(arrangeGrob(p1,p4, ncol=2),arrangeGrob(p3,p2,ncol=2),p6, heights=c(2,2,4), ncol=1)
#dev.off()

cds <- cds_mathys



##plot diagnosis differences over pseudotime
library(qvalue)
library(metaseqR)
plot_list<- list()

index=0
#'Ex0','Ex3','Ex4','Ex6','Ex7','Mic1','Oli0','Oli1','Opc0','Opc1'
for (subcluster in sort(unique(cds$Subcluster))){#unique(cds$Subcluster)
  if (length(cds$Subcluster[cds$Subcluster==subcluster])>=500){
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
    #cds_subset = cds[,cds$sex=='female']
    #cds_subset = cds[,cds$sex=='male']
    cds_subset = cds_subset[,cds_subset$Subcluster==subcluster] # mic3 subcluster
    cds_subset = cds_subset[rowSums(counts(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
    cds_subset$batch <- as.factor(cds_subset$batch)
    cds_subset <- align_cds(cds_subset, 
                            preprocess_method="PCA",
                            alignment_group="batch",
                            residual_model_formula_str="~educ+pmi")
    cds_subset = reduce_dimension(cds_subset)
    cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
    cds_subset$Diagnosis[cds_subset$Diagnosis=='Control']='Cont'
    cds_subset$Diagnosis[cds_subset$Diagnosis=='Late-path AD']='Late'
    cds_subset$Diagnosis[cds_subset$Diagnosis=='Early-path AD']='Early'
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

pdf(paste0("~/scAD_analysis/figures2/mathys_sn_diagnosis_LHpart1_louvain_batched",".pdf"))
do.call(grid.arrange, c(plot_list[1:15], ncol=3))
dev.off()


pdf(paste0("~/scAD_analysis/figures2/mathys_sn_diagnosis_LHpart2_louvain_batched",".pdf"))
do.call(grid.arrange, c(plot_list[16:length(plot_list)], ncol=3))
dev.off()



##plot diagnosis differences over pseudotime by sex:
#make sex-specific cds and clear the cds slots and rescale data
cdsF = cds_mathys[,cds_mathys$sex=='female']
dim(cdsF)
cdsF <- clear_cds_slots(cdsF)
cdsF$batch <- as.character(cdsF$batch)
cdsF$batch <- as.factor(cdsF$batch)
cdsF = preprocess_cds(cdsF, num_dim = 30,method="PCA",norm_method="size_only")
cdsF <- align_cds(cdsF, 
                        preprocess_method="PCA",
                        alignment_group="batch",
                        residual_model_formula_str="~educ+pmi")
cdsF = reduce_dimension(cdsF)
plot_cells(cdsF, color_cells_by="Subcluster")
cdsF$Diagnosis[cdsF$Diagnosis=='Control']='Cont'
cdsF$Diagnosis[cdsF$Diagnosis=='Late-path AD']='Late'
cdsF$Diagnosis[cdsF$Diagnosis=='Early-path AD']='Early'

#Male cds 
cdsM = cds_mathys[,cds_mathys$sex=='male']
dim(cdsM)
cdsM <- clear_cds_slots(cdsM)
cdsM$batch <- as.character(cdsM$batch)
cdsM$batch <- as.factor(cdsM$batch)
cdsM = preprocess_cds(cdsM, num_dim = 30,method="PCA",norm_method="size_only")
cdsM <- align_cds(cdsM, 
                  preprocess_method="PCA",
                  alignment_group="batch",
                  residual_model_formula_str="~educ+pmi")
cdsM = reduce_dimension(cdsM)
plot_cells(cdsM, color_cells_by="Subcluster")
cdsM$Diagnosis[cdsM$Diagnosis=='Control']='Cont'
cdsM$Diagnosis[cdsM$Diagnosis=='Late-path AD']='Late'
cdsM$Diagnosis[cdsM$Diagnosis=='Early-path AD']='Early'





#Calculate association between diagnosis (0/1) and pseudotime for each subcluster with >100 cells and plot, by sex:
#cds <- cdsF
cds <- cdsM
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
    #cds_subset = cds[,cds$sex=='female']
    #cds_subset = cds[,cds$sex=='male']
    cds_subset = cds_subset[,cds_subset$Subcluster==subcluster] # mic3 subcluster
    cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
    cds_subset$batch <- as.character(cds_subset$batch)
    cds_subset$batch <- as.factor(cds_subset$batch)
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
pdf(paste0("~/scAD_analysis/figures2/mathys_sn_diagnosis_LHpart1_M_batched",".pdf"))
do.call(grid.arrange, c(plot_list[1:18], ncol=3))
dev.off()


#pdf(paste0("~/scAD_analysis/figures2/mathys_sn_diagnosis_LHpart2_F_batched",".pdf"))
pdf(paste0("~/scAD_analysis/figures2/mathys_sn_diagnosis_LHpart2_M_batched",".pdf"))
do.call(grid.arrange, c(plot_list[18:length(plot_list)], ncol=3))
dev.off()


#plot diagnostic criterial (CERAD, Braak, cogdx, apoe genotype) for subtypes which show variation in diagnosis across pseudotime
plot_list<- list()
cds$APOE_genotype[is.na(cds$APOE_genotype)]='33'
index=0
for (subcluster in c('Mic1','Oli0','Oli1')){#,'Oli0','Oli1'
  index=index+1
  if (length(cds$Subcluster[cds$Subcluster==subcluster])<1000){
    cell_size=1}
  else if (length(cds$Subcluster[cds$Subcluster==subcluster])<3000){
    cell_size=.5}
  else {
    cell_size=.5
  }
  print(subcluster)
  cds_subset = cds[,cds$sex=='male']
  cds_subset = cds_subset[,cds_subset$broad.cell.type==gsub('[[:digit:]]+', '', subcluster)] # mic3 subcluster
  cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells
  #cds_subset = preprocess_cds(cds_subset, num_dim = 30,method="LSI",norm_method="size_only")
  cds_subset = cds_subset[,cds_subset$Subcluster==subcluster] # mic3 subcluster
  cds_subset = reduce_dimension(cds_subset)
  cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
  cds_subset$Diagnosis[cds_subset$Diagnosis=='Control']='Cont'
  cds_subset$Diagnosis[cds_subset$Diagnosis=='Late-path AD']='Late'
  cds_subset$Diagnosis[cds_subset$Diagnosis=='Early-path AD']='Early'
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
pdf(paste0("~/scAD_analysis/figures2/mathys_sn_diagnostic_criteria_male",".pdf"))
do.call(grid.arrange, c(plot_list, ncol=1))
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


