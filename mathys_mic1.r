#read in cell data set object that is preprocessed with mathys genes and celltypes
cds_mathys = readRDS(file = "~/scAD_analysis/Mathys_monocle_preprocessed_cds.rds")

library(ggplot2)
library(gridExtra)
library(MASS)
library(limma)
library(edgeR)
library(grid)
library(monocle3)
library(qvalue)
library(metaseqR)
cds_mathys$educ = as.numeric(cds_mathys$educ)
cds_mathys$Education = cds_mathys$educ
cds_mathys$Sex = cds_mathys$sex
cds_mathys$CERAD = cds_mathys$ceradsc
cds_mathys$Braak = cds_mathys$braaksc
cds_mathys$APOE_genotype = cds_mathys$apoe_genotype
cds_mathys$batch = as.factor(cds_mathys$batch)
cds_mathys$Diagnosis[cds_mathys$Diagnosis=='Control']='Cont'
cds_mathys$Diagnosis[cds_mathys$Diagnosis=='Late-path AD']='Late'
cds_mathys$Diagnosis[cds_mathys$Diagnosis=='Early-path AD']='Early'


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


plot_cells(cdsF, genes=c("MSN"),cell_size=0.1,label_cell_groups=0,show_trajectory_graph=FALSE)
plot_cells(cdsM, genes=c("MSN"),cell_size=0.1,label_cell_groups=0,show_trajectory_graph=FALSE)





###get pseudotime for one subcluster at a time, starting with Mic1 (starting with all-female or all-male cds generated in monocle_mathys_data_LH.R)

cds_subset = cdsF[,cdsF$Subcluster=='Oli5'] 
#cds_subset = cdsM[,cdsM$Subcluster=='Ast1'] 


cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells

dim(cds_subset)
#cds_subset$batch <- as.character(cds_subset$batch)
#cds_subset$batch <- as.factor(cds_subset$batch)
#cds_subset = preprocess_cds(cds_subset, num_dim = 30,method="PCA",norm_method="size_only")
cds_subset <- align_cds(cds_subset, 
                        preprocess_method="PCA",
                        alignment_group="batch",
                        residual_model_formula_str="~educ+pmi")
cds_subset = reduce_dimension(cds_subset)
cds_subset = cluster_cells(cds_subset, cluster_method="louvain")

#see what the mic1 plot looks like if we do not batch correct:
#cds_test <- cdsF[,cdsF$Subcluster=='Ast1']
cds_test <- cdsM[,cdsM$Subcluster=='Ast1']
cds_test = cds_test[rowSums(exprs(cds_subset) != 0) >= 20,]
cds_test <- clear_cds_slots(cds_test)
cds_test = preprocess_cds(cds_test, num_dim = 30,method="PCA",norm_method="size_only", residual_model_formula_str="~pmi+batch+educ")
cds_test = reduce_dimension(cds_test)
cds_test = cluster_cells(cds_test, cluster_method="louvain")

plot_cells(cds_test, color_cells_by="ros_ids",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

#for males, tweak clustering algorithm for continuous pseudotime tree
#cds_subset = cluster_cells(cds_subset, cluster_method="louvain", k=35)
plot_cells(cds_subset, color_cells_by="ros_ids",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

#There is a weird tiny cluster off by itself, need to remove:
# plot_cells(cds_subset, color_cells_by="cluster",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1)+theme(
#   legend.title = element_text(size = 10),
#   legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
#   legend.key.width = unit(0.1,"cm"))

#it is cluster 1
# cds_subset$cluster <- clusters(cds_subset, reduction_method="UMAP")
# #cds_subset <- cds_subset[,cds_subset$cluster!=1]
# plot_cells(cds_subset, color_cells_by="ros_ids")
# cds_subset <- align_cds(cds_subset, 
#                         preprocess_method="PCA",
#                         alignment_group="batch",
#                         residual_model_formula_str="~educ+pmi")
# cds_subset = reduce_dimension(cds_subset)
# cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
# plot_cells(cds_subset, color_cells_by="ros_ids",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1)+theme(
#   legend.title = element_text(size = 10),
#   legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
#   legend.key.width = unit(0.1,"cm"))

cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

cds_subset<-order_cells(cds_subset)
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))
   
cds_subset$pseudotime = pseudotime(cds_subset)
plot_cells(cds_subset, color_cells_by="pseudotime",cell_size=2)

pval_data <- data.frame(
  batch = colData(cds_subset)$batch,
  pmi = colData(cds_subset)$pmi,
  educ = colData(cds_subset)$educ,
  Diagnosis = cds_subset$Diagnosis,
  ros_ids = cds_subset$ros_ids)
pval_data$diagnosis2 <- ifelse(pval_data$Diagnosis=='Cont',0,1)
pval_data$pseudotime = cds_subset$pseudotime

fit <- lm(diagnosis2~pseudotime+pmi+educ, data=pval_data)
#fit <- lm(pseudotime~diagnosis2+pmi+educ, data=pval_data)

summary(fit)

#fit_lmm <- lmer(diagnosis2~pseudotime + pmi + educ + (1 | ros_ids), data=pval_data)
lmmfit <- lmer(pseudotime~diagnosis2 + pmi + educ + (1|ros_ids), data=pval_data)
#lmmfit <- lmer(pseudotime~diagnosis2 +(1|ros_ids), data=pval_data)

summary(lmmfit)






dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
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
pdf(paste0("~/scAD_analysis/figures2/mathys_Ast1_M_dx_pseudotime.pdf"))
p3<-grid.arrange(arrangeGrob(p1,p2,ncol=2,top = textGrob(paste0('Ast1 Males',stars),gp=gpar(fontsize=18,font=2))))
dev.off()



##### plot diagnostic criterial (CERAD, Braak, cogdx, apoe genotype) for subtypes which show variation in diagnosis across pseudotime
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
dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
fit1 <- lmFit(pval_data$CERAD,dm)
fit2 <- eBayes(fit1)
pval=topTable(fit2,coef=2)$adj.P.Val
stars=''
pval_data$CERAD = as.factor(cds_subset$CERAD)
dp <- ggplot(pval_data, aes(x=CERAD, y=pseudotime, fill=CERAD)) + 
    geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
p1<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    ggtitle(paste('p=',formatC(pval, format = "e", digits = 2),stars))+theme(axis.title=element_text(size=10))
  
  ### Braak ###
dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
fit1 <- lmFit(pval_data$Braak,dm)
fit2 <- eBayes(fit1)
pval=topTable(fit2,coef=2)$adj.P.Val
 
pval_data$Braak = as.factor(cds_subset$Braak)
dp <- ggplot(pval_data, aes(x=Braak, y=pseudotime, fill=Braak)) + 
    geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
p2<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    ggtitle(paste('p=',formatC(pval, format = "e", digits = 2),stars))+theme(axis.title=element_text(size=10))
 
   ### COGDX ###
dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
fit1 <- lmFit(pval_data$COGDX,dm)
fit2 <- eBayes(fit1)
pval=topTable(fit2,coef=2)$adj.P.Val
  
pval_data$COGDX = as.factor(cds_subset$cogdx)
dp <- ggplot(pval_data, aes(x=COGDX, y=pseudotime, fill=COGDX)) + 
    geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
p3<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    ggtitle(paste('p=',formatC(pval, format = "e", digits = 2),stars))+theme(axis.title=element_text(size=10))
 
   ### APOE_genotype ###
dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
fit1 <- lmFit(pval_data$APOE_genotype,dm)
fit2 <- eBayes(fit1)
pval=topTable(fit2,coef=2)$adj.P.Val

pval_data$APOE_genotype = as.factor(cds_subset$APOE_genotype)
dp <- ggplot(pval_data, aes(x=APOE_genotype, y=pseudotime, fill=APOE_genotype)) + 
    geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
p4<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    ggtitle(paste('p=',formatC(pval, format = "e", digits = 2),stars))+theme(axis.title=element_text(size=10))    
  
#pdf(paste0("~/scAD_analysis/figures2/mathys_Mic1_F_All_neuropath_pseudotime.pdf"))
#tiff(file="~/scAD_analysis/figures2/mathys_Mic0_F_neuropath_pseudotime.tiff", height=85,width=100,units='mm',res=300)

p5<-grid.arrange(arrangeGrob(p1,p2,p3,p4,ncol=4,top = textGrob(paste0('Ast1 M',stars),gp=gpar(fontsize=15,font=2))))
#dev.off()
p5

#Differentially expressed genes: genes associated with pseuodtime

gene_fits = fit_models(cds_subset, model_formula_str = "~pseudotime+pmi+educ")
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
write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='~/scAD_analysis/figures2/Ast1_M_degs.csv')




#Download csv file with significant gene x pseudotime associations and fill in genes of interest for plotting
#pdf("~/scAD_analysis/figures2/Oli5_All_top4DEGS.pdf")
p <- plot_cells(cds_subset, genes=c("PVALB", "COX7C", "GAD1", "PEBP1", "ATP1B1", "HSP90AA1", "MDH1", "CKB", "SPARCL1"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
#dev.off()
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=1,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)

top_downgenes <- c("PVALB", "COX7C", "GAD1", "PEBP1", "ATP1B1", "HSP90AA1", "MDH1", "CKB", "SPARCL1")
top_downgenes <- c("CKB", "COX4I1", "PEBP1", "NDUFA4", "GAPDH", "CYCS", "HSP90AA1", "MDH1", "CALM1")

Oli5_downs <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_downgenes,]
 #                              colData(cds_subset)$Subcluster %in% c("Mic1")]
pdf("~/scAD_analysis/figures2/In0_M_downgenes_x_pstime.pdf")
p <- plot_genes_in_pseudotime(Oli5_downs, color_cells_by="ros_ids",min_expr=0.5)
p
dev.off()


p <- plot_cells(cds_subset, genes=c("CNTNAP2","RBFOX1","SGCZ","ROBO2","DLG2","LRP1B", "DLGAP1", "GRIK1","NRXN3"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p <- plot_cells(cds_subset, genes=c("CNTNAP2","DLGAP1","ERBB4","NRG3","RBFOX1","DLG2", "TENM2", "DPP10","ROBO2"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)

p
dev.off()
top_upgenes <- c("CNTNAP2","RBFOX1","SGCZ","ROBO2","DLG2","LRP1B", "DLGAP1", "GRIK1","NRXN3")
top_upgenes <- c("CNTNAP2","DLGAP1","ERBB4","NRG3","RBFOX1","DLG2", "TENM2", "DPP10","ROBO2")

Oli5_ups <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_upgenes,]
 #                              colData(cds_subset)$Subcluster %in% c("Mic1")]
pdf("~/scAD_analysis/figures2/In0_M_upgenes_x_pseudotime.pdf")
p <- plot_genes_in_pseudotime(Oli5_ups, color_cells_by="ros_ids",min_expr=0.5)
p
dev.off()


#### finding genes that are differentially expressed on different paths through the trajectory ###
mic1_graphtest <- graph_test(cds_subset, neighbor_graph="principal_graph")
mic1_deg_ids <- subset(mic1_graphtest, q_value < .05)
write.csv(mic1_deg_ids,file='~/scAD_analysis/figures2/mic1_F_all_graphtestdegs.csv')
plot_cells(cds_subset, genes=c("PLEKHA7","RPL11","RPL7",'RPS4X'),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)

#cds_top4 <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_genes,]
pdf("~/scAD_analysis/figures2/Mic1_All_top4_graphtestDEGS.pdf")
p <- plot_cells(cds_subset, genes=c("FAM65B","C1QC", "RPL19", "APOE"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
dev.off()

top_genes1 <- c("FTL", "RPS20", "PLEKHA7", "RPL19", "FAM65B", "C1QC", "RPL23A", "APOE")
Mic1_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_genes1,
                               colData(cds_subset)$Subcluster %in% c("Mic1")]
pdf("~/scAD_analysis/figures2/Mic1_All_genesXpseudotime.pdf")
p <- plot_genes_in_pseudotime(Mic1_lineage_cds, color_cells_by="Diagnosis",min_expr=0.5)
p
dev.off()

pdf("~/scAD_analysis/figures2/Mic1_All_Diagnosis.pdf")
p <- plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
p
dev.off()

plot_cells(cds_subset, color_cells_by='pseudotime', cell_size=2)












#see what genes are differentially expressed in early pseudotime or late pseudotime
#calculating DE of genes in mic1 compared to all other cell pops
cds_early <- choose_cells(cds_subset)
cds_late <- choose_cells(cds_subset)

cdsF$Mic1 <- ifelse(cdsF$Subcluster=='Mic1',1,0)
gene_fits = fit_models(cds_subset, model_formula_str = "~pseudotime+pmi+educ")
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
write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='~/scAD_analysis/figures2/mic1F_all_degs.csv')




plot_cells(cds_subset, genes=c("MSN"),cell_size=2,label_cell_groups=0,show_trajectory_graph=TRUE)
gene_fits = fit_models(cds_subset, model_formula_str = "~pseudotime+pmi+batch+educ")
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
write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='~/scAD_analysis/figures2/mic1_F_all_degs.csv')

pdf("~/scAD_analysis/figures2/Mic1_All_top4DEGS.pdf")
p <- plot_cells(cds_subset, genes=c("SRGN","RPL27A","HIST1H2AC",'RPS6'),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
dev.off()
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=1,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)



dim(cds_subset)


#### finding genes that are differentially expressed on different paths through the trajectory ###
mic1_graphtest <- graph_test(cds_subset, neighbor_graph="principal_graph")
mic1_deg_ids <- subset(mic1_graphtest, q_value < .05)
write.csv(mic1_deg_ids,file='~/scAD_analysis/figures2/mic1_all_graphtestdegs.csv')
plot_cells(cds_subset, genes=c("FAM65B","C1QC","RPL23A",'RPL32'),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)










cds_subset = cdsF[,cdsF$Subcluster=='Oli5'] 

cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells

dim(cds_subset)
cds_subset$batch <- as.character(cds_subset$batch)
cds_subset$batch <- as.factor(cds_subset$batch)
cds_subset <- align_cds(cds_subset, 
                        preprocess_method="PCA",
                        alignment_group="batch",
                        residual_model_formula_str="~educ+pmi")
cds_subset = reduce_dimension(cds_subset)
cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
plot_cells(cds_subset, genes="FTH1",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)











cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells

dim(cds_subset)
cds_subset$batch <- as.character(cds_subset$batch)
cds_subset$batch <- as.factor(cds_subset$batch)
cds_subset <- align_cds(cds_subset, 
                        preprocess_method="PCA",
                        alignment_group="batch",
                        residual_model_formula_str="~educ+pmi")
cds_subset = reduce_dimension(cds_subset)
cds_subset = cluster_cells(cds_subset, cluster_method="louvain")

#see what the mic1 plot looks like if we do not batch correct:
# cds_test <- cdsF[,cdsF$Subcluster=='Mic1'] 
# cds_test <- clear_cds_slots(cds_test)
# cds_test = preprocess_cds(cds_test, num_dim = 30,method="PCA",norm_method="size_only", residual_model_formula_str="~pmi+batch+educ")
# cds_test = reduce_dimension(cds_test)
# cds_test = cluster_cells(cds_test, cluster_method="louvain")
# 
# plot_cells(cds_test, color_cells_by="ros_ids",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1)+theme(
#   legend.title = element_text(size = 10),
#   legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
#   legend.key.width = unit(0.1,"cm"))
# 
#for males, tweak clustering algorithm for continuous pseudotime tree
#cds_subset = cluster_cells(cds_subset, cluster_method="louvain", k=35)
plot_cells(cds_subset, color_cells_by="ros_ids",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

cds_subset<-order_cells(cds_subset)
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),legend.key.size = unit(.1, "cm"),
  legend.key.width = unit(0.1,"cm"))

cds_subset$pseudotime = pseudotime(cds_subset)
plot_cells(cds_subset, color_cells_by="pseudotime",cell_size=2)

pval_data <- data.frame(
  batch = colData(cds_subset)$batch,
  pmi = colData(cds_subset)$pmi,
  educ = colData(cds_subset)$educ,
  Diagnosis = cds_subset$Diagnosis,
  ros_ids = cds_subset$ros_ids)
pval_data$diagnosis2 <- ifelse(pval_data$Diagnosis=='Cont',0,1)
pval_data$pseudotime = cds_subset$pseudotime

#fit_lmm <- lmer(diagnosis2~pseudotime + (pseudotime | ros_ids), data=pval_data)
fit <- lm(diagnosis2~pseudotime+pmi+educ, data=pval_data)

summary(fit)


dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
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
pdf(paste0("~/scAD_analysis/figures2/mathys_Oli5_F_dx_pseudotime.pdf"))
p3<-grid.arrange(arrangeGrob(p1,p2,ncol=2,top = textGrob(paste0('Oli5 Females',stars),gp=gpar(fontsize=18,font=2))))
dev.off()



##### plot diagnostic criterial (CERAD, Braak, cogdx, apoe genotype) for subtypes which show variation in diagnosis across pseudotime
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
dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
fit1 <- lmFit(pval_data$CERAD,dm)
fit2 <- eBayes(fit1)
pval=topTable(fit2,coef=2)$adj.P.Val
stars=''
pval_data$CERAD = as.factor(cds_subset$CERAD)
dp <- ggplot(pval_data, aes(x=CERAD, y=pseudotime, fill=CERAD)) + 
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
p1<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  ggtitle(paste('p=',formatC(pval, format = "e", digits = 2),stars))+theme(axis.title=element_text(size=10))

### Braak ###
dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
fit1 <- lmFit(pval_data$Braak,dm)
fit2 <- eBayes(fit1)
pval=topTable(fit2,coef=2)$adj.P.Val

pval_data$Braak = as.factor(cds_subset$Braak)
dp <- ggplot(pval_data, aes(x=Braak, y=pseudotime, fill=Braak)) + 
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
p2<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  ggtitle(paste('p=',formatC(pval, format = "e", digits = 2),stars))+theme(axis.title=element_text(size=10))

### COGDX ###
dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
fit1 <- lmFit(pval_data$COGDX,dm)
fit2 <- eBayes(fit1)
pval=topTable(fit2,coef=2)$adj.P.Val

pval_data$COGDX = as.factor(cds_subset$cogdx)
dp <- ggplot(pval_data, aes(x=COGDX, y=pseudotime, fill=COGDX)) + 
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
p3<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  ggtitle(paste('p=',formatC(pval, format = "e", digits = 2),stars))+theme(axis.title=element_text(size=10))

### APOE_genotype ###
dm <- model.matrix(~pseudotime+pmi+educ, data=pval_data)
fit1 <- lmFit(pval_data$APOE_genotype,dm)
fit2 <- eBayes(fit1)
pval=topTable(fit2,coef=2)$adj.P.Val

pval_data$APOE_genotype = as.factor(cds_subset$APOE_genotype)
dp <- ggplot(pval_data, aes(x=APOE_genotype, y=pseudotime, fill=APOE_genotype)) + 
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")
p4<-dp+theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  ggtitle(paste('p=',formatC(pval, format = "e", digits = 2),stars))+theme(axis.title=element_text(size=10))    

#pdf(paste0("~/scAD_analysis/figures2/mathys_Mic1_F_All_neuropath_pseudotime.pdf"))
tiff(file="~/scAD_analysis/figures2/mathys_Mic1_F_All_neuropath_pseudotime.tiff", height=85,width=100,units='mm',res=300)

p5<-grid.arrange(arrangeGrob(p1,p2,p3,p4,ncol=4,top = textGrob(paste0('Oli5 All',stars),gp=gpar(fontsize=15,font=2))))
dev.off()


#Differentially expressed genes: genes associated with pseuodtime

gene_fits = fit_models(cds_subset, model_formula_str = "~pseudotime+pmi+educ")
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
write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='~/scAD_analysis/figures2/Oli5F_all_degs.csv')

#pdf("~/scAD_analysis/figures2/Oli5_All_top4DEGS.pdf")
p <- plot_cells(cds_subset, genes=c("DPYD", "FRMD5", "SLC24A2", "CTNNA3", "FRMD4B"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
dev.off()
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=1,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)

top_downgenes <- c("DPYD", "FRMD5", "SLC24A2", "CTNNA3", "FRMD4B")
Oli5_downs <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_downgenes,]
#                              colData(cds_subset)$Subcluster %in% c("Mic1")]
pdf("~/scAD_analysis/figures2/Oli5_F_downgenes_x_pstime.pdf")
p <- plot_genes_in_pseudotime(Oli5_downs, color_cells_by="ros_ids",min_expr=0.5)
p
dev.off()


p <- plot_cells(cds_subset, genes=c("FTH1","PTMA","FTL","QDPR","EEF1A1","PLP1", "TSC22D4", "MARCKSL1","CNP"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
dev.off()
top_upgenes <- c("FTH1","PTMA","FTL","QDPR","EEF1A1","PLP1", "TSC22D4", "MARCKSL1","CNP")
Oli5_ups <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_upgenes,]
#                              colData(cds_subset)$Subcluster %in% c("Mic1")]
pdf("~/scAD_analysis/figures2/Oli5_All_genesXpseudotime.pdf")
p <- plot_genes_in_pseudotime(Oli5_ups, color_cells_by="ros_ids",min_expr=0.5)
p
dev.off()


#### finding genes that are differentially expressed on different paths through the trajectory ###
mic1_graphtest <- graph_test(cds_subset, neighbor_graph="principal_graph")
mic1_deg_ids <- subset(mic1_graphtest, q_value < .05)
write.csv(mic1_deg_ids,file='~/scAD_analysis/figures2/mic1_F_all_graphtestdegs.csv')
plot_cells(cds_subset, genes=c("PLEKHA7","RPL11","RPL7",'RPS4X'),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)

#cds_top4 <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_genes,]
pdf("~/scAD_analysis/figures2/Mic1_All_top4_graphtestDEGS.pdf")
p <- plot_cells(cds_subset, genes=c("FAM65B","C1QC", "RPL19", "APOE"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
dev.off()

top_genes1 <- c("FTL", "RPS20", "PLEKHA7", "RPL19", "FAM65B", "C1QC", "RPL23A", "APOE")
Mic1_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% top_genes1,
                               colData(cds_subset)$Subcluster %in% c("Mic1")]
pdf("~/scAD_analysis/figures2/Mic1_All_genesXpseudotime.pdf")
p <- plot_genes_in_pseudotime(Mic1_lineage_cds, color_cells_by="Diagnosis",min_expr=0.5)
p
dev.off()

pdf("~/scAD_analysis/figures2/Mic1_All_Diagnosis.pdf")
p <- plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=2,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
p
dev.off()

plot_cells(cds_subset, color_cells_by='pseudotime', cell_size=2)












#see what genes are differentially expressed in early pseudotime or late pseudotime
#calculating DE of genes in mic1 compared to all other cell pops
cds_early <- choose_cells(cds_subset)
cds_late <- choose_cells(cds_subset)

cdsF$Mic1 <- ifelse(cdsF$Subcluster=='Mic1',1,0)
gene_fits = fit_models(cds_subset, model_formula_str = "~pseudotime+pmi+educ")
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
write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='~/scAD_analysis/figures2/mic1F_all_degs.csv')




plot_cells(cds_subset, genes=c("MSN"),cell_size=2,label_cell_groups=0,show_trajectory_graph=TRUE)
gene_fits = fit_models(cds_subset, model_formula_str = "~pseudotime+pmi+batch+educ")
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
write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='~/scAD_analysis/figures2/mic1_F_all_degs.csv')

pdf("~/scAD_analysis/figures2/Mic1_All_top4DEGS.pdf")
p <- plot_cells(cds_subset, genes=c("SRGN","RPL27A","HIST1H2AC",'RPS6'),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p
dev.off()
plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=1,label_cell_groups=FALSE, label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)



dim(cds_subset)


#### finding genes that are differentially expressed on different paths through the trajectory ###
mic1_graphtest <- graph_test(cds_subset, neighbor_graph="principal_graph")
mic1_deg_ids <- subset(mic1_graphtest, q_value < .05)
write.csv(mic1_deg_ids,file='~/scAD_analysis/figures2/mic1_all_graphtestdegs.csv')
plot_cells(cds_subset, genes=c("FAM65B","C1QC","RPL23A",'RPL32'),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)










cds_subset = cdsF[,cdsF$Subcluster=='Oli5'] 

cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 20 cells

dim(cds_subset)
cds_subset$batch <- as.character(cds_subset$batch)
cds_subset$batch <- as.factor(cds_subset$batch)
cds_subset <- align_cds(cds_subset, 
                        preprocess_method="PCA",
                        alignment_group="batch",
                        residual_model_formula_str="~educ+pmi")
cds_subset = reduce_dimension(cds_subset)
cds_subset = cluster_cells(cds_subset, cluster_method="louvain")
plot_cells(cds_subset, genes="FTH1",cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)