mvnchol_bamlss <- function(k, type = c("basic", "modified", "chol"), ...) {
  type <- match.arg(type)
  switch(type, basic = mvn_chol(k = k, ...), chol = mvn_chol(k = k, ...), modified = mvn_modchol(k = k, ...))
}
mvn_modchol <- function(k = 2, ...) {
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  k_phi <- k * (k - 1)/2
  links <- c(rep("identity", k), rep("log", k), rep("identity", k_phi))
  names(links) <- c(mu, innov, phi)
  rval <- list(family = "mvnmodchol", names = c(mu, innov, phi), links = links, d = function(y, par, log = FALSE) {
 d <- log_dmvnmodchol(y, par)
 if (!log) d <- exp(d)
 return(d)
  })
  mu_score_calls <- paste0("function(y, par, ...) {mu_score_mvnmodchol(y, par, j = ", seq_len(k), ")}")
  innov_score_calls <- paste0("function(y, par, ...) {innov_score_mvnmodchol(y, par, j = ", seq_len(k), ")}")
  phi_score_calls <- utils::combn(seq_len(k), 2, function(x) {
 paste0("function(y, par, ...) {phi_score_mvnmodchol(y, par, i = ", x[1], ", j = ", x[2], ")}")
  })
  scores <- list()
  for (j in seq_along(mu)) {
 scores[[mu[j]]] <- eval(parse(text = mu_score_calls[j]))
  }
  for (j in seq_along(innov)) {
 scores[[innov[j]]] <- eval(parse(text = innov_score_calls[j]))
  }
  for (j in seq_along(phi)) {
 scores[[phi[j]]] <- eval(parse(text = phi_score_calls[j]))
  }
  rval$score <- scores
  mu_hess_calls <- paste0("function(y, par, ...) {mu_hess_mvnmodchol(y, par, j = ", seq_len(k), ")}")
  innov_hess_calls <- paste0("function(y, par, ...) {innov_hess_mvnmodchol(y, par, j = ", seq_len(k), ")}")
  phi_hess_calls <- utils::combn(seq_len(k), 2, function(x) {
 paste0("function(y, par, ...) {phi_hess_mvnmodchol(y, par, i = ", x[1], ", j = ", x[2], ")}")
  })
  hesses <- list()
  for (j in seq_along(mu)) {
 hesses[[mu[j]]] <- eval(parse(text = mu_hess_calls[j]))
  }
  for (j in seq_along(innov)) {
 hesses[[innov[j]]] <- eval(parse(text = innov_hess_calls[j]))
  }
  for (j in seq_along(phi)) {
 hesses[[phi[j]]] <- eval(parse(text = phi_hess_calls[j]))
  }
  rval$hess <- hesses
  rval$precision <- function(par) {
 k <- length(grep("innov", names(par)))
 n <- length(par[[1]])
 Linvt <- list()
 Siginv <- list()
 for (i in 1:n) {
   Linvt[[i]] <- matrix(0, nrow = k, ncol = k)
   for (j in 1:k) {
  Linvt[[i]][j, j] <- 1/sqrt(par[[paste0("innov", j)]][i])
  if (j < k) {
  for (l in (j + 1):k) {
    Linvt[[i]][j, l] <- -par[[paste0("phi", j, l)]][i]/sqrt(par[[paste0("innov", l)]][i])
  }
  }
   }
   Siginv[[i]] <- Linvt[[i]] %*% t(Linvt[[i]])
 }
 return(Siginv)
  }
  rval$covariance <- function(par) {
 k <- length(grep("innov", names(par)))
 n <- length(par[[1]])
 Linvt <- list()
 Sig <- list()
 for (i in 1:n) {
   Linvt[[i]] <- matrix(0, nrow = k, ncol = k)
   for (j in 1:k) {
  Linvt[[i]][j, j] <- 1/sqrt(par[[paste0("innov", j)]][i])
  if (j < k) {
  for (l in (j + 1):k) {
    Linvt[[i]][j, l] <- -par[[paste0("phi", j, l)]][i]/sqrt(par[[paste0("innov", l)]][i])
  }
  }
   }
   L <- solve(Linvt[[i]])
   Sig[[i]] <- t(L) %*% L
 }
 return(Sig)
  }
  rval$r <- function(par, expand = 1) {
 n <- length(par[[1]])
 k <- length(grep("innov", names(par)))
 Sigs <- rval$covariance(par)
 simdat <- matrix(ncol = k, nrow = n * expand)
 for (i in 1:n) {
   mu_vec <- vector(length = k)
   for (j in 1:k) {
  mu_vec[j] <- par[[paste0("mu", j)]][i]
   }
   simdat[seq_len(expand) + (i - 1) * expand, ] <- mvtnorm::rmvnorm(n = expand, mean = mu_vec, sigma = Sigs[[i]])
 }
 return(simdat)
  }
  class(rval) <- "family.bamlss"
  rval
}
log_dmvnmodchol_ref <- function(y, par) {
  n <- nrow(y)
  k <- ncol(y)
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  term_1 <- -k/2 * log(2 * pi)
  term_2 <- -1/2 * apply(log(subset(par_full, select = innov)), 1, sum)
  term_3 <- vector(length = n)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  for (ni in 1:n) {
 y_tild <- y_til[ni, ]
 L_inv <- matrix(0, nrow = k, ncol = k)
 for (l in 1:length(phi)) {
   i <- utils::combn(k, 2)[1, l]
   j <- utils::combn(k, 2)[2, l]
   L_inv[j, i] <- -par[[paste0("phi", i, j)]][ni]/sqrt(par[[paste0("innov", j)]][ni])
 }
 for (i in 1:length(innov)) {
   L_inv[i, i] <- 1/sqrt(par[[paste0("innov", i)]][ni])
 }
 term_3[ni] <- -1/2 * norm(L_inv %*% y_tild, type = "2")^2
  }
  ll <- term_1 + term_2 + term_3
  return(ll)
}
log_dmvnmodchol_C <- function(y, par) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  .Call("log_dmvnmodcholC", y, par, n, k, PACKAGE = "bamlss")
}
log_dmvnmodchol <- log_dmvnmodchol_C
mu_score_mvnmodchol_ref <- function(y, par, j) {
  n <- nrow(y)
  k <- ncol(y)
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dmu <- vector(length = n)
  for (ni in 1:n) {
 y_tild <- y_til[ni, ]
 L_inv <- matrix(0, nrow = k, ncol = k)
 temp <- utils::combn(k, 2)
 for (l in 1:length(phi)) {
   ii <- temp[1, l]
   jj <- temp[2, l]
   L_inv[jj, ii] <- -par[[paste0("phi", ii, jj)]][ni]/sqrt(par[[paste0("innov", jj)]][ni])
 }
 for (ii in 1:length(innov)) {
   L_inv[ii, ii] <- 1/sqrt(par[[paste0("innov", ii)]][ni])
 }
 Sigma_inv <- t(L_inv) %*% L_inv
 dl_dmu[ni] <- sum(Sigma_inv[j, ] * y_tild)
  }
  return(dl_dmu)
}
mu_score_mvnmodchol_C <- function(y, par, j) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  j <- as.integer(j)
  .Call("mu_score_mvnmodcholC", y, par, n, k, j, PACKAGE = "bamlss")
}
mu_score_mvnmodchol <- mu_score_mvnmodchol_C
innov_score_mvnmodchol_ref <- function(y, par, j) {
  n <- nrow(y)
  k <- ncol(y)
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dinnov <- vector(length = n)
  for (ni in 1:n) {
 y_tild <- y_til[ni, ]
 Phimat <- matrix(0, nrow = k, ncol = k)
 diag(Phimat) <- -1
 temp <- utils::combn(k, 2)
 for (l in 1:length(phi)) {
   ii <- temp[1, l]
   jj <- temp[2, l]
   Phimat[jj, ii] <- par[[paste0("phi", ii, jj)]][ni]
 }
 dl_dinnov[ni] <- -1/2 + 1/(2 * par[[paste0("innov", j)]][ni]) * sum(y_tild[1:j] * Phimat[j, 1:j])^2
  }
  return(dl_dinnov)
}
innov_score_mvnmodchol_C <- function(y, par, j) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  j <- as.integer(j)
  .Call("innov_score_mvnmodcholC", y, par, n, k, j, PACKAGE = "bamlss")
}
innov_score_mvnmodchol <- innov_score_mvnmodchol_C
phi_score_mvnmodchol_ref <- function(y, par, i, j) {
  n <- nrow(y)
  k <- ncol(y)
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dphi <- vector(length = n)
  for (ni in 1:n) {
 y_tild <- y_til[ni, ]
 Tmat <- matrix(0, nrow = k, ncol = k)
 diag(Tmat) <- 1
 temp <- utils::combn(k, 2)
 for (l in 1:length(phi)) {
   ii <- temp[1, l]
   jj <- temp[2, l]
   Tmat[jj, ii] <- -par[[paste0("phi", ii, jj)]][ni]
 }
 dl_dphi[ni] <- y_tild[i]/par[[paste0("innov", j)]][ni] * sum(y_tild[1:j] * Tmat[j, 1:j])
  }
  return(dl_dphi)
}
phi_score_mvnmodchol_C <- function(y, par, i, j) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  i <- as.integer(i)
  j <- as.integer(j)
  .Call("phi_score_mvnmodcholC", y, par, n, k, i, j, PACKAGE = "bamlss")
}
phi_score_mvnmodchol <- phi_score_mvnmodchol_C
mu_hess_mvnmodchol_ref <- function(y, par, j) {
  n <- nrow(y)
  k <- ncol(y)
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  dl_dmu <- vector(length = n)
  for (ni in 1:n) {
 dl_dmu[ni] <- 1/par[[paste0("innov", j)]][ni]
 if (j < k) {
   for (jj in (j + 1):k) {
  dl_dmu[ni] <- dl_dmu[ni] + par[[paste0("phi", j, jj)]][ni]^2/par[[paste0("innov", jj)]][ni]
   }
 }
  }
  return(dl_dmu)
}
mu_hess_mvnmodchol <- mu_hess_mvnmodchol_ref
innov_hess_mvnmodchol_ref <- function(y, par, j) {
  n <- nrow(y)
  k <- ncol(y)
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dinnov <- vector(length = n)
  for (ni in 1:n) {
 y_tild <- y_til[ni, ]
 Phimat <- matrix(0, nrow = k, ncol = k)
 diag(Phimat) <- -1
 temp <- utils::combn(k, 2)
 for (l in 1:length(phi)) {
   ii <- temp[1, l]
   jj <- temp[2, l]
   Phimat[jj, ii] <- par[[paste0("phi", ii, jj)]][ni]
 }
 dl_dinnov[ni] <- 0.5/par[[paste0("innov", j)]][ni] * sum(y_tild[1:j] * Phimat[j, 1:j])^2
  }
  return(dl_dinnov)
}
innov_hess_mvnmodchol <- innov_hess_mvnmodchol_ref
phi_hess_mvnmodchol_ref <- function(y, par, i, j) {
  n <- nrow(y)
  k <- ncol(y)
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dphi <- vector(length = n)
  for (ni in 1:n) {
 dl_dphi[ni] <- y_til[ni, i]^2/par[[paste0("innov", j)]][ni]
  }
  return(dl_dphi)
}
phi_hess_mvnmodchol <- phi_hess_mvnmodchol_ref
make_formula <- function(formula, type = "basic") {
  FORM <- Formula::as.Formula(formula)
  l <- length(FORM)
  k <- l[1]
  seq_k <- seq_len(k)
  if (l[2] == 1) {
 FORM <- stats::update(FORM, . ~ . | 1 | 1)
  }
  if (l[2] == 2) {
 FORM <- stats::update(FORM, . ~ . | . | 1)
  }
  fam <- mvnchol_bamlss(k = k, type = type)
  nms <- fam$names
  rval <- c(lapply(seq_k, formula, x = FORM, rhs = 1), rep(list(formula(FORM, lhs = FALSE, rhs = 2)), k), rep(list(formula(FORM, lhs = FALSE, rhs = 3)), k * (k - 1)/2))
  names(rval) <- nms
  rval
}
dist_mvnchol <- function(k, r = k - 1, type = c("basic", "modified", "chol"), ...) {
  type <- match.arg(type)
  switch(type, basic = dist_mvn_chol(k = k, r = r, ...), chol = dist_mvn_chol(k = k, r = r, ...), modified = dist_mvn_modchol(k = k, r = r, ...))
}
dist_mvn_modchol <- function(k, r = k - 1, ...) {
  nms_mu <- paste0("mu_", seq_len(k))
  nms_innov <- paste0("innov_", seq_len(k))
  combns <- utils::combn(k, 2)
  combns_cut <- combns[, which((combns[2, ] - combns[1, ]) <= r)]
  nms_phi <- paste("phi", combns_cut[1, ], combns_cut[2, ], sep = "_")
  k_phi <- ncol(combns_cut)
  k_all <- k + k + k_phi
  nms_par <- c(nms_mu, nms_innov, nms_phi)
  nms_eta <- c(nms_mu, paste0("log(", nms_innov, ")"), nms_phi)
  ddist <- function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {
 n <- nrow(y)
 term_1 <- -k/2 * log(2 * pi)
 term_2 <- -1/2 * sum(eta[paste0("log(", nms_innov, ")")])
 term_3 <- numeric(n)
 y_til <- y - matrix(rep(eta[nms_mu], n), ncol = k, byrow = TRUE)
 for (ni in seq_len(n)) {
   y_tild <- y_til[ni, ]
   L_inv <- matrix(0, nrow = k, ncol = k)
   for (l in seq_len(k_phi)) {
  i <- combns_cut[1, l]
  j <- combns_cut[2, l]
  L_inv[j, i] <- -eta[[paste0("phi_", i, "_", j)]]/sqrt(exp(eta[[paste0("log(innov_", j, ")")]]))
   }
   for (i in seq_len(k)) {
  L_inv[i, i] <- 1/sqrt(exp(eta[[paste0("log(innov_", i, ")")]]))
   }
   term_3[ni] <- -1/2 * norm(L_inv %*% y_tild, type = "2")^2
 }
 ll <- term_1 + term_2 + term_3
 if (isTRUE(log)) {
   if (isFALSE(sum)) {
  return(ll)
   }
   else {
  return(sum(ll))
   }
 }
 else {
   if (isFALSE(sum)) {
  return(exp(ll))
   }
   else {
  message("You are summing the (non-log) likelihoods!!!")
  return(sum(exp(ll)))
   }
 }
  }
  sdist_ref <- function(y, eta, weights = NULL, sum = FALSE) {
 n <- nrow(y)
 scores_mat <- matrix(ncol = k_all, nrow = n)
 colnames(scores_mat) <- nms_eta
 y_til <- y - matrix(rep(eta[nms_mu], n), ncol = k, byrow = TRUE)
 for (ni in seq_len(n)) {
   L_inv <- matrix(0, nrow = k, ncol = k)
   Phi_mat <- matrix(0, nrow = k, ncol = k)
   diag(Phi_mat) <- -1
   for (l in seq_len(k_phi)) {
  i <- combns_cut[1, l]
  j <- combns_cut[2, l]
  L_inv[j, i] <- -eta[[paste0("phi_", i, "_", j)]]/sqrt(exp(eta[[paste0("log(innov_", j, ")")]]))
  Phi_mat[j, i] <- eta[[paste0("phi_", i, "_", j)]]
   }
   for (i in seq_len(k)) {
  L_inv[i, i] <- 1/sqrt(exp(eta[[paste0("log(innov_", i, ")")]]))
   }
   Sigma_inv <- t(L_inv) %*% L_inv
   for (j in seq_len(k)) {
  scores_mat[ni, j] <- sum(Sigma_inv[j, ] * y_til[ni, ])
   }
   for (j in seq_len(k)) {
  scores_mat[ni, j + k] <- -0.5 + 1/(2 * exp(eta[[paste0("log(innov_", j, ")")]])) * sum(y_til[ni, 1:j] * Phi_mat[j, 1:j])^2
   }
   for (l in seq_len(k_phi)) {
  i <- combns_cut[1, l]
  j <- combns_cut[2, l]
  scores_mat[ni, 2 * k + l] <- -y_til[ni, i]/exp(eta[[paste0("log(innov_", j, ")")]]) * sum(y_til[ni, 1:j] * Phi_mat[j, 1:j])
   }
 }
 if (isFALSE(sum)) {
   return(scores_mat)
 }
 else {
   return(colSums(scores_mat))
 }
  }
  sdist_C <- function(y, eta, weights = NULL, sum = FALSE) {
 y <- as.matrix(y)
 storage.mode(y) <- "numeric"
 n <- nrow(y)
 k <- ncol(y)
 scores_mat <- matrix(ncol = k_all, nrow = n)
 for (j in 1:k) {
   scores_mat[, nms_mu[j]] <- .Call("mu_score_mvnmodcholC", y, eta, n, k, j, PACKAGE = "bamlss")
 }
 for (j in 1:k) {
   scores_mat[, paste0("log(innov_", j, ")")] <- .Call("innov_score_mvnmodcholC", y, eta, n, k, j, PACKAGE = "bamlss")
 }
 if (isFALSE(sum)) {
   return(scores_mat)
 }
 else {
   return(colSums(scores_mat))
 }
  }
  sdist <- sdist_ref
  link <- c(rep("identity", k), rep("log", k), rep("identity", k_phi))
  linkfun <- function(par) {
 eta <- par
 eta[seq_len(k) + k] <- log(par[seq_len(k) + k])
 names(eta) <- nms_eta
 return(eta)
  }
  linkinv <- function(eta) {
 par <- eta
 par[seq_len(k) + k] <- exp(eta[seq_len(k) + k])
 names(par) <- nms_par
 return(par)
  }
  linkinvdr <- function(eta) {
 dpardeta <- rep(1, k_all)
 dpardeta[seq_len(k) + k] <- exp(eta[seq_len(k) + k])
 names(dpardeta) <- nms_par
 return(dpardeta)
  }
  startfun <- function(y, weights = NULL) {
 n <- nrow(y)
 if (is.null(weights) || (length(weights) == 0)) {
   starteta <- rep(-999, k_all)
   names(starteta) <- nms_eta
   starteta[nms_mu] <- colMeans(y)
   y_til <- y - matrix(rep(starteta[nms_mu], n), ncol = k, byrow = TRUE)
   starteta[["log(innov_1)"]] <- log(stats::var(y_til[, 1]))
   for (i in 2:k) {
  X <- y_til[, max(1, r - k):(i - 1)]
  beta_vec <- solve(t(X) %*% X) %*% t(X) %*% y_til[, i]
  starteta[paste0("phi_", max(1, r - k):(i - 1), "_", i)] <- beta_vec
  starteta[paste0("log(innov_", i, ")")] <- log(stats::var(y_til[, i] - as.vector(X %*% beta_vec)))
   }
 }
 else {
   stop("Weights not implemented for startfun")
 }
 return(starteta)
  }
  mle <- TRUE
  dist_list <- list(family.name = "MVN Cholesky", ddist = ddist, sdist = sdist, link = link, linkfun = linkfun, linkinv = linkinv, linkinvdr = linkinvdr, startfun = startfun, mle = mle, gamlssobj = FALSE, censored = FALSE)
  class(dist_list) <- "disttree.family"
  return(dist_list)
}
dist_mvn_chol <- function(k, ...) {
  stop("distree family for basic Cholesky is not implemented yet!")
}
mvn_chol <- function(k = 2, ...) {
  mu <- paste0("mu", seq_len(k))
  lamdiag <- paste0("lamdiag", seq_len(k))
  lambda <- utils::combn(seq_len(k), 2, function(x) paste0("lambda", x[1], x[2]))
  k_lambda <- k * (k - 1)/2
  links <- c(rep("identity", k), rep("log", k), rep("identity", k_lambda))
  names(links) <- c(mu, lamdiag, lambda)
  rval <- list(family = "mvnchol", names = c(mu, lamdiag, lambda), links = links, d = function(y, par, log = FALSE) {
 d <- log_dmvnchol(y, par)
 if (!log) d <- exp(d)
 return(d)
  })
  mu_score_calls <- paste0("function(y, par, ...) {mu_score_mvnchol(y, par, j = ", seq_len(k), ")}")
  lamdiag_score_calls <- paste0("function(y, par, ...) {lamdiag_score_mvnchol(y, par, j = ", seq_len(k), ")}")
  lambda_score_calls <- utils::combn(seq_len(k), 2, function(x) {
 paste0("function(y, par, ...) {lambda_score_mvnchol(y, par, i = ", x[1], ", j = ", x[2], ")}")
  })
  scores <- list()
  for (j in seq_along(mu)) {
 scores[[mu[j]]] <- eval(parse(text = mu_score_calls[j]))
  }
  for (j in seq_along(lamdiag)) {
 scores[[lamdiag[j]]] <- eval(parse(text = lamdiag_score_calls[j]))
  }
  for (j in seq_along(lambda)) {
 scores[[lambda[j]]] <- eval(parse(text = lambda_score_calls[j]))
  }
  rval$score <- scores
  rval$precision <- function(par) {
 k <- length(grep("lamdiag", names(par)))
 n <- length(par[[1]])
 Linvt <- list()
 Siginv <- list()
 for (i in 1:n) {
   Linvt[[i]] <- matrix(0, nrow = k, ncol = k)
   for (j in 1:k) {
  Linvt[[i]][j, j] <- par[[paste0("lamdiag", j)]][i]
  if (j < k) {
  for (l in (j + 1):k) {
    Linvt[[i]][j, l] <- par[[paste0("lambda", j, l)]][i]
  }
  }
   }
   Siginv[[i]] <- Linvt[[i]] %*% t(Linvt[[i]])
 }
 return(Siginv)
  }
  rval$covariance <- function(par) {
 k <- length(grep("lamdiag", names(par)))
 n <- length(par[[1]])
 Linvt <- list()
 Sig <- list()
 for (i in 1:n) {
   Linvt[[i]] <- matrix(0, nrow = k, ncol = k)
   for (j in 1:k) {
  Linvt[[i]][j, j] <- par[[paste0("lamdiag", j)]][i]
  if (j < k) {
  for (l in (j + 1):k) {
    Linvt[[i]][j, l] <- par[[paste0("lambda", j, l)]][i]
  }
  }
   }
   L <- solve(Linvt[[i]])
   Sig[[i]] <- t(L) %*% L
 }
 return(Sig)
  }
  rval$correlation <- function(par) {
 k <- length(grep("lamdiag", names(par)))
 n <- length(par[[1]])
 Linvt <- list()
 Cor <- list()
 for (i in 1:n) {
   Linvt[[i]] <- matrix(0, nrow = k, ncol = k)
   for (j in 1:k) {
  Linvt[[i]][j, j] <- par[[paste0("lamdiag", j)]][i]
  if (j < k) {
  for (l in (j + 1):k) {
    Linvt[[i]][j, l] <- par[[paste0("lambda", j, l)]][i]
  }
  }
   }
   L <- solve(Linvt[[i]])
   Sig <- t(L) %*% L
   inv_sds <- matrix(0, nrow = k, ncol = k)
   diag(inv_sds) <- 1/sqrt(diag(Sig))
   Cor[[i]] <- inv_sds %*% Sig %*% inv_sds
 }
 return(Cor)
  }
  rval$means <- function(par) {
 k <- length(grep("lamdiag", names(par)))
 n <- length(par[[1]])
 lst <- lapply(seq_len(n), function(x) lapply(par[seq_len(k)], `[[`, i = x))
 lapply(lst, do.call, what = "c")
  }
  rval$stdev <- function(par) {
 k <- length(grep("lamdiag", names(par)))
 n <- length(par[[1]])
 Linvt <- list()
 Sig <- list()
 for (i in 1:n) {
   Linvt[[i]] <- matrix(0, nrow = k, ncol = k)
   for (j in 1:k) {
  Linvt[[i]][j, j] <- par[[paste0("lamdiag", j)]][i]
  if (j < k) {
  for (l in (j + 1):k) {
    Linvt[[i]][j, l] <- par[[paste0("lambda", j, l)]][i]
  }
  }
   }
   L <- solve(Linvt[[i]])
   Sig[[i]] <- t(L) %*% L
 }
 lapply(Sig, function(x) stats::setNames(sqrt(diag(x)), nm = paste0("sig", seq_len(k))))
  }
  class(rval) <- "family.bamlss"
  rval
}
log_dmvnchol_ref <- function(y, par) {
  n <- nrow(y)
  k <- ncol(y)
  mu <- paste0("mu", seq_len(k))
  lamdiag <- paste0("lamdiag", seq_len(k))
  lambda <- utils::combn(seq_len(k), 2, function(x) paste0("lambda", x[1], x[2]))
  par_full <- do.call("cbind", par)
  term_1 <- -k/2 * log(2 * pi)
  term_2 <- apply(log(subset(par_full, select = lamdiag)), 1, sum)
  term_3 <- vector(length = n)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  for (ni in 1:n) {
 y_tild <- y_til[ni, ]
 L_inv <- matrix(0, nrow = k, ncol = k)
 for (l in 1:length(lambda)) {
   i <- utils::combn(k, 2)[1, l]
   j <- utils::combn(k, 2)[2, l]
   L_inv[j, i] <- par[[paste0("lambda", i, j)]][ni]
 }
 for (i in 1:length(lamdiag)) {
   L_inv[i, i] <- par[[paste0("lamdiag", i)]][ni]
 }
 term_3[ni] <- -1/2 * norm(L_inv %*% y_tild, type = "2")^2
  }
  ll <- term_1 + term_2 + term_3
  return(ll)
}
log_dmvnchol_C <- function(y, par) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  .Call("log_dmvncholC", y, par, n, k, PACKAGE = "bamlss")
}
log_dmvnchol <- log_dmvnchol_C
mu_score_mvnchol_ref <- function(y, par, j) {
  n <- nrow(y)
  k <- ncol(y)
  mu <- paste0("mu", seq_len(k))
  lamdiag <- paste0("lamdiag", seq_len(k))
  lambda <- utils::combn(seq_len(k), 2, function(x) paste0("lambda", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dmu <- vector(length = n)
  for (ni in 1:n) {
 y_tild <- y_til[ni, ]
 L_inv <- matrix(0, nrow = k, ncol = k)
 temp <- utils::combn(k, 2)
 for (l in 1:length(lambda)) {
   ii <- temp[1, l]
   jj <- temp[2, l]
   L_inv[jj, ii] <- par[[paste0("lambda", ii, jj)]][ni]
 }
 for (ii in 1:length(lamdiag)) {
   L_inv[ii, ii] <- par[[paste0("lamdiag", ii)]][ni]
 }
 Sigma_inv <- t(L_inv) %*% L_inv
 dl_dmu[ni] <- sum(Sigma_inv[j, ] * y_tild)
  }
  return(dl_dmu)
}
mu_score_mvnchol_C <- function(y, par, j) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  j <- as.integer(j)
  .Call("mu_score_mvncholC", y, par, n, k, j, PACKAGE = "bamlss")
}
mu_score_mvnchol <- mu_score_mvnchol_C
lamdiag_score_mvnchol_ref <- function(y, par, j) {
  n <- nrow(y)
  k <- ncol(y)
  mu <- paste0("mu", seq_len(k))
  lamdiag <- paste0("lamdiag", seq_len(k))
  lambda <- utils::combn(seq_len(k), 2, function(x) paste0("lambda", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dlamdiag <- vector(length = n)
  for (ni in 1:n) {
 y_tild <- y_til[ni, ]
 L_inv <- matrix(0, nrow = k, ncol = k)
 temp <- utils::combn(k, 2)
 for (l in 1:length(lambda)) {
   ii <- temp[1, l]
   jj <- temp[2, l]
   L_inv[jj, ii] <- par[[paste0("lambda", ii, jj)]][ni]
 }
 for (ii in 1:length(lamdiag)) {
   L_inv[ii, ii] <- par[[paste0("lamdiag", ii)]][ni]
 }
 dl_dlamdiag[ni] <- 1 - L_inv[j, j] * y_tild[j] * sum(y_tild[1:j] * L_inv[j, 1:j])
  }
  return(dl_dlamdiag)
}
lamdiag_score_mvnchol_C <- function(y, par, j) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  j <- as.integer(j)
  .Call("lamdiag_score_mvncholC", y, par, n, k, j, PACKAGE = "bamlss")
}
lamdiag_score_mvnchol <- lamdiag_score_mvnchol_C
lambda_score_mvnchol_ref <- function(y, par, i, j) {
  n <- nrow(y)
  k <- ncol(y)
  mu <- paste0("mu", seq_len(k))
  lamdiag <- paste0("lamdiag", seq_len(k))
  lambda <- utils::combn(seq_len(k), 2, function(x) paste0("lambda", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dlambda <- vector(length = n)
  for (ni in 1:n) {
 y_tild <- y_til[ni, ]
 L_inv <- matrix(0, nrow = k, ncol = k)
 temp <- utils::combn(k, 2)
 for (l in 1:length(lambda)) {
   ii <- temp[1, l]
   jj <- temp[2, l]
   L_inv[jj, ii] <- par[[paste0("lambda", ii, jj)]][ni]
 }
 for (ii in 1:length(lamdiag)) {
   L_inv[ii, ii] <- par[[paste0("lamdiag", ii)]][ni]
 }
 dl_dlambda[ni] <- -y_tild[i] * sum(y_tild[1:j] * L_inv[j, 1:j])
  }
  return(dl_dlambda)
}
lambda_score_mvnchol_C <- function(y, par, i, j) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  i <- as.integer(i)
  j <- as.integer(j)
  .Call("lambda_score_mvncholC", y, par, n, k, i, j, PACKAGE = "bamlss")
}
lambda_score_mvnchol <- lambda_score_mvnchol_C
