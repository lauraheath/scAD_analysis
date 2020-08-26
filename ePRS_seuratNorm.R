#starting with the counts matrix derived in snRNAhackathon_UWdata_filtered.R:
library(Seurat)
library(scater)
library(batchelor)

counts <- dat[rowSums(dat != 0) >= 500,]
dim(counts)
#write.csv(as.matrix(counts),file='/Users/relyanow/Desktop/mount2/users/relyanow/test_data/UW_sn_AD.csv')

dat2 <- CreateSeuratObject(counts = as.matrix(counts), project = "HighLow", min.cells = 1, min.features = 3)
counts<-0
dat2

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
dat2[["percent.mt"]] <- PercentageFeatureSet(dat2, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(dat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(dat2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dat2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

dat3 <- subset(dat2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 15)
dat3 <- NormalizeData(dat3, normalization.method = "LogNormalize", scale.factor = 10000)
str(dat3)
dat3 <- FindVariableFeatures(object = dat3,nfeatures = 10000)
str(dat3@assays$RNA@var.features)

dat3 <- ScaleData(dat3,features=dat3@assays$RNA@var.features)

dat3 <- RunPCA(dat3,npcs=50)

# Examine and visualize PCA results a few different ways
print(dat3[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(dat3, dims = 1:2, reduction = "pca")

DimPlot(dat3, reduction = "pca")

ElbowPlot(dat3,50)

dat3 <- FindNeighbors(dat, dims = 1:25)
dat <- FindClusters(dat, resolution = 0.25)






#Apply sctransform normalization
#Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.