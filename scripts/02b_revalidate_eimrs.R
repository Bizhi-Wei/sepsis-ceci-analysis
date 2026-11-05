#!/usr/bin/env Rscript
# Revalidate the full EIMRS development procedure without held-out screening.
suppressPackageStartupMessages({library(limma);library(survival);library(glmnet);library(pROC);library(jsonlite)})
if (!nzchar(Sys.getenv("CECI_WORK_DIR"))) stop("Set CECI_WORK_DIR to a private working folder.")
base <- normalizePath(Sys.getenv("CECI_WORK_DIR"), mustWork = TRUE)
repair <- file.path(base, '03_results/eimrs_validation')
dir.create(repair, recursive=TRUE, showWarnings=FALSE)
script_arg <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub('^--file=', '', script_arg[1]), mustWork = TRUE))
source(file.path(script_dir, '02a_eimrs_validation.R'))
set.seed(20261101)
geo <- readRDS(file.path(base,'01_data/processed/geo_annotated_cohorts.rds'))
locked <- read.csv(file.path(base, '03_results/tables/Table_4_Locked_EIMRS_Signature_Genes.csv'))
universe <- read.csv(file.path(base,'03_results/tables/Table_S8_EGIH_Univariate_Cox_GSE65682.csv'))$gene
ii <- which(geo$pheno65$cohort_group=='Sepsis_28d_Prognostic')
hh <- which(geo$pheno65$cohort_group=='Healthy_Control')
xx <- geo$expr65_gene[geo$shared_genes, ii, drop=FALSE]
healthy <- geo$expr65_gene[geo$shared_genes, hh, drop=FALSE]
time <- geo$pheno65$time_28d[ii]; event <- geo$pheno65$mort_28d[ii]
n <- length(event)
cat('Full derivation fit with all screening and scaling repeated\n');flush.console()
full <- eimrs_fit_pipeline(xx,healthy,time,event,seed=20261026,universe=universe)
lp_full <- eimrs_predict_pipeline(full,xx)
locked_beta <- setNames(locked$lasso_beta_per_SD,locked$gene)
weight_match <- setequal(names(locked_beta),names(full$beta)) && max(abs(locked_beta[names(full$beta)]-full$beta))<1e-8
if (!weight_match) stop('Full-data fit changed locked coefficients; update final model and downstream validation coherently before publishing')
concord <- function(tt,ee,lp) survival::concordance(survival::Surv(tt,ee)~lp,reverse=TRUE)$concordance
c_apparent <- concord(time,event,lp_full)
outer <- eimrs_stratified_folds(event,10,seed=20261101)
oof <- rep(NA_real_,n);fold_records <- list()
for (k in 1:10) {
  tr <- which(outer!=k);te <- which(outer==k)
  fit <- eimrs_fit_pipeline(xx[,tr,drop=FALSE],healthy,time[tr],event[tr],
                            seed=20261101+k,original_ids=tr,universe=universe)
  oof[te] <- eimrs_predict_pipeline(fit,xx[,te,drop=FALSE])
  fold_records[[k]] <- data.frame(Fold=k,N_train=length(tr),N_test=length(te),Events_test=sum(event[te]),
                                Candidates=fit$n_candidates,Selected=length(fit$beta),Lambda_rule=fit$lambda_rule)
  cat(sprintf('Outer fold %d/10: %d training-only candidates, %d genes\n',k,fit$n_candidates,length(fit$beta)));flush.console()
}
stopifnot(all(is.finite(oof)))
oof_auc <- as.numeric(pROC::auc(pROC::roc(event,oof,direction='<',quiet=TRUE)))
oof_c <- concord(time,event,oof)
write.csv(do.call(rbind,fold_records),file.path(repair,'outer_fold_audit.csv'),row.names=FALSE)
saveRDS(list(full=full,oof=oof,outer=outer,oof_auc=oof_auc,oof_c=oof_c),file.path(repair,'nested_cv_results.rds'))
cat(sprintf('Nested OOF AUC %.6f; C %.6f\n',oof_auc,oof_c));flush.console()
B <- as.integer(Sys.getenv('EIMRS_BOOTSTRAPS','1000'))
workers <- min(8L,as.integer(Sys.getenv('EIMRS_WORKERS','8')))
cl <- parallel::makePSOCKcluster(workers,outfile=file.path(repair,'bootstrap_workers.log'))
on.exit(try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
parallel::clusterEvalQ(cl,{Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1');suppressPackageStartupMessages({library(limma);library(survival);library(glmnet)});NULL})
parallel::clusterExport(cl,c('xx','healthy','time','event','n','universe','eimrs_fit_pipeline','eimrs_predict_pipeline','eimrs_stratified_folds','concord'),envir=environment())
boot_fun <- function(b) {
  set.seed(70261101+b)
  ids <- sample.int(n,n,replace=TRUE);hi <- sample.int(ncol(healthy),ncol(healthy),replace=TRUE)
  tryCatch({
    fit <- eimrs_fit_pipeline(xx[,ids,drop=FALSE],healthy[,hi,drop=FALSE],time[ids],event[ids],
                              seed=90261101+b,original_ids=ids,universe=universe)
    train_lp <- eimrs_predict_pipeline(fit,xx[,ids,drop=FALSE]);original_lp <- eimrs_predict_pipeline(fit,xx)
    data.frame(Replicate=b,Optimism=concord(time[ids],event[ids],train_lp)-concord(time,event,original_lp),
               Candidates=fit$n_candidates,Selected=length(fit$beta),Error='')
  },error=function(e)data.frame(Replicate=b,Optimism=NA_real_,Candidates=NA_integer_,Selected=NA_integer_,Error=conditionMessage(e)))
}
all_boot <- list();batch_size <- 40L
for (start in seq(1L,B,by=batch_size)) {
  ids <- start:min(B,start+batch_size-1L)
  result <- parallel::parLapplyLB(cl,ids,boot_fun)
  all_boot <- c(all_boot,result)
  boot <- do.call(rbind,all_boot)
  write.csv(boot,file.path(repair,'bootstrap_checkpoint.csv'),row.names=FALSE)
  cat(sprintf('Full-pipeline bootstrap %d/%d completed, failures=%d\n',nrow(boot),B,sum(!is.finite(boot$Optimism))));flush.console()
}
parallel::stopCluster(cl);cl <- NULL
valid <- is.finite(boot$Optimism)
if (sum(valid)<.95*B) stop('Too many bootstrap failures; retained checkpoints require investigation')
optimism <- mean(boot$Optimism[valid]);corrected <- c_apparent-optimism
boot_se <- stats::sd(boot$Optimism[valid])/sqrt(sum(valid))
summary <- list(Version='20261101',N=n,Events=sum(event),Full_model_genes=length(full$beta),
                Full_model_matches_locked_coefficients=weight_match,Full_model_rule=full$lambda_rule,
                Outer_folds=10,Inner_folds=10,Screening='All outcome-dependent screening repeated in each training resample',
                Bootstrap_replicates_requested=B,Bootstrap_replicates_valid=sum(valid),Bootstrap_failures=sum(!valid),
                Bootstrap_duplicate_policy='All copies of a patient stay in the same inner tuning fold',
                Healthy_reference_policy='Fixed healthy reference for outer CV; independently resampled healthy controls for bootstrap',
                Apparent_C=c_apparent,Mean_optimism=optimism,Corrected_C=corrected,Optimism_MCSE=boot_se,
                OOF_AUC=oof_auc,OOF_C=oof_c)
jsonlite::write_json(summary,file.path(repair,'validation_summary.json'),auto_unbox=TRUE,pretty=TRUE,digits=15)
input_md5 <- unname(tools::md5sum(file.path(base,'01_data/processed/geo_annotated_cohorts.rds')))
saveRDS(list(summary=summary,full_fit=full,oof=oof,outer_fold=outer,bootstrap=boot,input_md5=input_md5),
        file.path(base,'01_data/processed/eimrs_full_pipeline_validation.rds'))
cat(jsonlite::toJSON(summary,auto_unbox=TRUE,pretty=TRUE,digits=8),'\n')
