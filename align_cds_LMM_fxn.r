#' Align cells from different groups within a cds
#'
#' @description Data sets that contain cells from different groups often
#' benefit from alignment to subtract differences between them. Alignment
#' can be used to remove batch effects, subtract the effects of treatments,
#' or even potentially compare across species.
#' \code{align_cds} executes alignment and stores these adjusted coordinates.
#'
#' This function can be used to subtract both continuous and discrete batch
#' effects. For continuous effects, \code{align_cds} fits a linear model to the
#' cells' PCA or LSI coordinates and subtracts them using Limma. For discrete
#' effects, you must provide a grouping of the cells, and then these groups are
#' aligned using Batchelor, a "mutual nearest neighbor" algorithm described in:
#'
#' Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018). "Batch effects in
#' single-cell RNA-sequencing data are corrected by matching mutual nearest
#' neighbors." Nat. Biotechnol., 36(5), 421-427. doi: 10.1038/nbt.4091
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param preprocess_method a string specifying the low-dimensional space
#'   in which to perform alignment, currently either PCA or LSI. Default is
#'   "PCA".
#' @param residual_model_formula_str NULL or a string model formula specifying
#'   any effects to subtract from the data before dimensionality reduction.
#'   Uses a linear model to subtract effects. For non-linear effects, use
#'   alignment_group. Default is NULL.
#' @param alignment_group String specifying a column of colData to use for
#'  aligning groups of cells. The column specified must be a factor.
#'  Alignment can be used to subtract batch effects in a non-linear way.
#'  For correcting continuous effects, use residual_model_formula_str.
#'  Default is NULL.
#' @param alignment_k The value of k used in mutual nearest neighbor alignment
#' @param verbose Whether to emit verbose output during dimensionality
#'   reduction
#' @param ... additional arguments to pass to limma::lmFit if
#'   residual_model_formula is not NULL
#' @return an updated cell_data_set object
#' @export
align_cds_LMM <- function(cds,
                      preprocess_method = c("PCA", "LSI"),
                      alignment_group=NULL,
                      alignment_k=20,
                      residual_model_formula_str=NULL,
                      verbose=FALSE,
                      ...){
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(preprocess_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "preprocess_method must be one of 'PCA' or 'LSI'")
  preprocess_method <- match.arg(preprocess_method)
  
  preproc_res <- reducedDims(cds)[[preprocess_method]]
  assertthat::assert_that(!is.null(preproc_res),
                          msg = paste0("Preprocessing for '",
                                       preprocess_method, "' does not exist. ",
                                       "Please make sure you have run ",
                                       "preprocess_cds with",
                                       "preprocess_method = '",
                                       preprocess_method,
                                       "' before calling align_cds."))
  
  if (!is.null(residual_model_formula_str)) {
    if (verbose) message("Removing residual effects")
    X.model_mat <- Matrix::sparse.model.matrix(
      stats::as.formula(residual_model_formula_str),
      data = colData(cds),
      drop.unused.levels = TRUE)
    

    #fit <- lme4::lmer(..., data=Matrix::t(preproc_res))
    #beta <- coef(summary(fit))[, -1, drop=FALSE]
    corfit <- limma::duplicateCorrelation(Matrix::t(preproc_res), X.model_mat, block=cds$UWID)
    fit <- limma::lmFit(Matrix::t(preproc_res), X.model_mat, block = cds$UWID, correlation=corfit$consensus, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    preproc_res <- Matrix::t(as.matrix(Matrix::t(preproc_res)) -
                               beta %*% Matrix::t(X.model_mat[, -1]))
    cds@preprocess_aux$beta = beta
  }
  
  if(!is.null(alignment_group)) {
    message(paste("Aligning cells from different batches using Batchelor.",
                  "\nPlease remember to cite:\n\t Haghverdi L, Lun ATL,",
                  "Morgan MD, Marioni JC (2018). 'Batch effects in",
                  "single-cell RNA-sequencing data are corrected by matching",
                  "mutual nearest neighbors.' Nat. Biotechnol., 36(5),",
                  "421-427. doi: 10.1038/nbt.4091"))
    corrected_PCA = batchelor::fastMNN(as.matrix(preproc_res),
                                       batch=colData(cds)[,alignment_group],
                                       k=alignment_k,
                                       cos.norm=FALSE,
                                       pc.input = TRUE)
    preproc_res = corrected_PCA$corrected
    cds <- add_citation(cds, "MNN_correct")
  }
  
  reducedDims(cds)[["Aligned"]] <- as.matrix(preproc_res)
  
  cds
}






#' Preprocess a cds to prepare for trajectory inference
#'
#' @description Most analyses (including trajectory inference, and clustering)
#' in Monocle3, require various normalization and preprocessing steps.
#' \code{preprocess_cds} executes and stores these preprocessing steps.
#'
#' Specifically, depending on the options selected, \code{preprocess_cds} first
#' normalizes the data by log and size factor to address depth differences, or
#' by size factor only. Next, \code{preprocess_cds} calculates a lower
#' dimensional space that will be used as the input for further dimensionality
#' reduction like tSNE and UMAP.
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param method a string specifying the initial dimension method to use,
#'   currently either PCA or LSI. For LSI (latent semantic indexing), it
#'   converts the (sparse) expression matrix into tf-idf matrix and then
#'   performs SVD to decompose the gene expression / cells into certain
#'   modules / topics. Default is "PCA".
#' @param num_dim the dimensionality of the reduced space.
#' @param norm_method Determines how to transform expression values prior to
#'   reducing dimensionality. Options are "log", "size_only", and "none".
#'   Default is "log". Users should only use "none" if they are confident that
#'   their data is already normalized.
#' @param use_genes NULL or a list of gene IDs. If a list of gene IDs, only
#'   this subset of genes is used for dimensionality reduction. Default is
#'   NULL.
#' @param residual_model_formula_str NULL or a string model formula specifying
#'   any effects to subtract from the data before dimensionality reduction.
#'   Uses a linear model to subtract effects. For non-linear effects, use
#'   alignment_group. Default is NULL.
#' @param alignment_group String specifying a column of colData to use for
#'  aligning groups of cells. The column specified must be a factor.
#'  Alignment can be used to subtract batch effects in a non-linear way.
#'  For correcting continuous effects, use residual_model_formula_str.
#'  Default is NULL.
#' @param pseudo_count NULL or the amount to increase expression values before
#'   normalization and dimensionality reduction. If NULL (default), a
#'   pseudo_count of 1 is added for log normalization and 0 is added for size
#'   factor only normalization.
#' @param scaling When this argument is set to TRUE (default), it will scale
#'   each gene before running trajectory reconstruction. Relevant for
#'   method = PCA only.
#' @param verbose Whether to emit verbose output during dimensionality
#'   reduction
#' @param ... additional arguments to pass to limma::lmFit if
#'   residual_model_formula is not NULL
#' @return an updated cell_data_set object
#' @export
preprocess_cds <- function(cds, method = c('PCA', "LSI"),
                           num_dim=50,
                           norm_method = c("log", "size_only", "none"),
                           use_genes = NULL,
                           residual_model_formula_str=NULL,
                           alignment_group=NULL,
                           pseudo_count=NULL,
                           scaling = TRUE,
                           verbose=FALSE,
                           ...) {
  
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be one of 'PCA' or 'LSI'")
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(norm_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "norm_method must be one of 'log', 'size_only' or 'none'")
  assertthat::assert_that(assertthat::is.count(num_dim))
  if(!is.null(use_genes)) {
    assertthat::assert_that(is.character(use_genes))
    assertthat::assert_that(all(use_genes %in% row.names(rowData(cds))),
                            msg = paste("use_genes must be NULL, or all must",
                                        "be present in the row.names of rowData(cds)"))
  }
  assertthat::assert_that(!is.null(size_factors(cds)),
                          msg = paste("You must call estimate_size_factors before calling",
                                      "preprocess_cds."))
  assertthat::assert_that(sum(is.na(size_factors(cds))) == 0,
                          msg = paste("One or more cells has a size factor of",
                                      "NA."))
  
  method <- match.arg(method)
  norm_method <- match.arg(norm_method)
  
  #ensure results from RNG sensitive algorithms are the same on all calls
  set.seed(2016)
  FM <- normalize_expr_data(cds, norm_method, pseudo_count)
  
  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }
  
  if (!is.null(use_genes)) {
    FM <- FM[use_genes, ]
  }
  
  fm_rowsums = Matrix::rowSums(FM)
  FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]
  
  if(method == 'PCA') {
    if (verbose) message("Remove noise by PCA ...")
    
    irlba_res <- sparse_prcomp_irlba(Matrix::t(FM),
                                     n = min(num_dim,min(dim(FM)) - 1),
                                     center = scaling, scale. = scaling)
    preproc_res <- irlba_res$x
    row.names(preproc_res) <- colnames(cds)
    
    irlba_rotation <- irlba_res$rotation
    row.names(irlba_rotation) <- rownames(FM)
    cds@preprocess_aux$gene_loadings <- irlba_rotation
    cds@preprocess_aux$prop_var_expl <- irlba_res$sdev^2 / sum(irlba_res$sdev^2)
    
  } else if(method == "LSI") {
    
    preproc_res <- tfidf(FM)
    irlba_res <- irlba::irlba(Matrix::t(preproc_res),
                              nv = min(num_dim,min(dim(FM)) - 1))
    
    preproc_res <- irlba_res$u %*% diag(irlba_res$d)
    row.names(preproc_res) <- colnames(cds)
    
    irlba_rotation = irlba_res$v
    row.names(irlba_rotation) = rownames(FM)
    cds@preprocess_aux$gene_loadings = irlba_rotation
    
  }
  
  row.names(preproc_res) <- colnames(cds)
  
  reducedDims(cds)[[method]] <- as.matrix(preproc_res)
  cds@preprocess_aux$beta = NULL
  
  cds
}


# Helper function to normalize the expression data prior to dimensionality
# reduction
normalize_expr_data <- function(cds,
                                norm_method = c("log", "size_only", "none"),
                                pseudo_count = NULL) {
  norm_method <- match.arg(norm_method)
  
  FM <- SingleCellExperiment::counts(cds)
  
  # If we're going to be using log, and the user hasn't given us a
  # pseudocount set it to 1 by default.
  if (is.null(pseudo_count)){
    if(norm_method == "log")
      pseudo_count <- 1
    else
      pseudo_count <- 0
  }
  
  if (norm_method == "log") {
    # If we are using log, normalize by size factor before log-transforming
    
    FM <- Matrix::t(Matrix::t(FM)/size_factors(cds))
    
    if (pseudo_count != 1 || is_sparse_matrix(SingleCellExperiment::counts(cds)) == FALSE){
      FM <- FM + pseudo_count
      FM <- log2(FM)
    } else {
      FM@x = log2(FM@x + 1)
    }
    
  } else if (norm_method == "size_only") {
    FM <- Matrix::t(Matrix::t(FM)/size_factors(cds))
    FM <- FM + pseudo_count
  }
  return (FM)
}

# Andrew's tfidf
tfidf <- function(count_matrix, frequencies=TRUE, log_scale_tf=TRUE,
                  scale_factor=100000, block_size=2000e6) {
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf <- Matrix::t(Matrix::t(count_matrix) / Matrix::colSums(count_matrix))
  } else {
    # "raw count" method
    tf <- count_matrix
  }
  
  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }
  
  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(count_matrix) / Matrix::rowSums(count_matrix > 0))
  
  # Try to just to the multiplication and fall back on delayed array
  # TODO hopefully this actually falls back and not get jobs killed in SGE
  tf_idf_counts = tryCatch({
    tf_idf_counts = tf * idf
    tf_idf_counts
  }, error = function(e) {
    print(paste("TF*IDF multiplication too large for in-memory, falling back",
                "on DelayedArray."))
    options(DelayedArray.block.size=block_size)
    DelayedArray:::set_verbose_block_processing(TRUE)
    
    tf = DelayedArray::DelayedArray(tf)
    idf = as.matrix(idf)
    
    tf_idf_counts = tf * idf
    tf_idf_counts
  })
  
  rownames(tf_idf_counts) = rownames(count_matrix)
  colnames(tf_idf_counts) = colnames(count_matrix)
  tf_idf_counts = methods::as(tf_idf_counts, "sparseMatrix")
  return(tf_idf_counts)
}