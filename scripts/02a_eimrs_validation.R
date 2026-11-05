# Training-only EIMRS screening, scaling, tuning and refitting.
# The caller supplies the septic training rows and the fixed healthy reference;
# no held-out septic expression or outcome enters these operations.
eimrs_stratified_folds <- function(event, k = 10L, seed = 20261026L, ids = seq_along(event)) {
  set.seed(seed)
  out <- integer(length(event))
  unique_ids <- unique(ids)
  ev <- event[match(unique_ids, ids)]
  assigned <- integer(length(unique_ids))
  for (value in c(0, 1)) {
    pos <- which(ev == value)
    assigned[pos] <- sample(rep(seq_len(k), length.out = length(pos)))
  }
  out[] <- assigned[match(ids, unique_ids)]
  out
}

eimrs_fit_pipeline <- function(case_expression, healthy_expression, time, event,
                               seed = 20261026L, inner_folds = 10L,
                               original_ids = seq_along(event),
                               universe = rownames(case_expression)) {
  stopifnot(ncol(case_expression) == length(time), length(time) == length(event),
            identical(rownames(case_expression), rownames(healthy_expression)))
  genes <- intersect(universe, rownames(case_expression))
  x_raw <- t(case_expression[genes, , drop = FALSE])
  mu <- colMeans(x_raw)
  sd <- apply(x_raw, 2, stats::sd)
  valid <- is.finite(sd) & sd > 0
  x_raw <- x_raw[, valid, drop = FALSE]
  genes <- genes[valid]; mu <- mu[valid]; sd <- sd[valid]
  x <- scale(x_raw, center = mu, scale = sd)
  y <- survival::Surv(time, event)
  p_cox <- vapply(seq_len(ncol(x)), function(j) {
    fit <- tryCatch(suppressWarnings(survival::coxph(y ~ x[, j])), error = function(e) NULL)
    if (is.null(fit)) return(1)
    val <- summary(fit)$coefficients[1, 5]
    if (is.finite(val)) val else 1
  }, numeric(1))
  fdr_cox <- stats::p.adjust(p_cox, method = 'BH')
  de_matrix <- cbind(healthy_expression, case_expression)
  group <- factor(c(rep('Healthy', ncol(healthy_expression)), rep('Sepsis', ncol(case_expression))),
                  levels = c('Healthy', 'Sepsis'))
  de <- limma::eBayes(limma::lmFit(de_matrix, stats::model.matrix(~group)))
  # Sepsis-vs-healthy BH adjustment uses the entire measured shared gene universe,
  # as in the original development analysis; the caller may pass its full matrix.
  de_p <- de$p.value[, 2]
  de_fdr <- stats::p.adjust(de_p, method = 'BH')[match(genes, rownames(de_matrix))]
  lfc <- de$coefficients[match(genes, rownames(de_matrix)), 2]
  chosen <- which(de_fdr < .05 & abs(lfc) > .30 & fdr_cox < .05)
  fallback <- FALSE
  if (length(chosen) < 15) {
    chosen <- which(de_fdr < .05 & p_cox < .005)
    fallback <- TRUE
  }
  chosen <- chosen[order(p_cox[chosen])]
  if (length(chosen) < 2) {
    return(list(mu = mu, sd = sd, genes = character(), beta = numeric(),
                n_candidates = length(chosen), lambda_rule = 'no eligible model',
                fallback = fallback, constant = TRUE))
  }
  xx <- x[, chosen, drop = FALSE]
  foldid <- eimrs_stratified_folds(event, inner_folds, seed, original_ids)
  cv <- glmnet::cv.glmnet(xx, y, family = 'cox', alpha = 1, foldid = foldid,
                          type.measure = 'C', parallel = FALSE)
  nz_1se <- sum(as.matrix(stats::coef(cv, s = 'lambda.1se'))[, 1] != 0)
  nz_min <- sum(as.matrix(stats::coef(cv, s = 'lambda.min'))[, 1] != 0)
  if (nz_1se >= 6 && nz_1se <= 18) {
    lambda <- cv$lambda.1se; rule <- 'lambda.1se'
  } else if (nz_min >= 6 && nz_min <= 18) {
    lambda <- cv$lambda.min; rule <- 'lambda.min'
  } else {
    nz <- colSums(as.matrix(cv$glmnet.fit$beta) != 0)
    okay <- which(nz >= 10 & nz <= 15)
    if (length(okay)) {
      best <- okay[which.max(cv$cvm[okay])]
      lambda <- cv$lambda[best]; rule <- sprintf('CV choice among 10-15 genes (%d genes)', nz[best])
    } else { lambda <- cv$lambda.min; rule <- 'lambda.min' }
  }
  beta <- as.matrix(stats::coef(cv, s = lambda))[, 1]
  beta <- beta[beta != 0]
  list(mu = mu, sd = sd, genes = names(beta), beta = beta,
       n_candidates = length(chosen), lambda_rule = rule,
       fallback = fallback, constant = length(beta) == 0)
}

eimrs_predict_pipeline <- function(fit, new_expression) {
  if (!length(fit$genes)) return(rep(0, ncol(new_expression)))
  xx <- scale(t(new_expression[fit$genes, , drop = FALSE]),
              center = fit$mu[fit$genes], scale = fit$sd[fit$genes])
  as.numeric(xx %*% fit$beta)
}
