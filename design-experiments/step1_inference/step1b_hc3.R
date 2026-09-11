# Follow-up: heteroscedasticity/imbalance broke every candidate. Does a sandwich (HC3) fix it,
# and what does it cost elsewhere? No permutations needed -> fast.
suppressMessages({library(limma); library(parallel); library(sandwich)})
src <- readLines("step1_inference.R"); end <- grep("^# -+ per-taxon test engine", src) - 1
eval(parse(text = src[seq_len(end)]))          # settings + simulate_cell
n_rep_null <- 20; n_rep_alt <- 15; m <- 500
test_hc <- function(y, x, z) {
  ok <- is.finite(y); y <- y[ok]; xx <- x[ok]; zz <- z[ok]; n <- length(y); df <- n - 3
  if (df < 2 || length(unique(xx)) < 2) return(NULL)
  f <- lm(y ~ zz + xx); b <- unname(coef(f)["xx"])
  s2 <- sum(resid(f)^2) / df; vxx <- summary(f)$cov.unscaled["xx","xx"]
  se_hc3 <- sqrt(vcovHC(f, type = "HC3")["xx","xx"])
  c(beta = b, s2 = s2, df = df, vxx = vxx,
    p_ols = 2*pt(-abs(b/sqrt(s2*vxx)), df),
    p_hc3 = 2*pt(-abs(b/se_hc3), df))
}
run_rep <- function(n, kind, null_cell, rep_id) {
  sim <- simulate_cell(n, kind, null_cell)
  res <- lapply(seq_len(m), function(j) test_hc(sim$Y[, j], sim$x, sim$z))
  keep <- !vapply(res, is.null, logical(1)); R <- do.call(rbind, res[keep]); is_da <- sim$is_da[keep]
  sq <- squeezeVar(R[,"s2"], df = R[,"df"])
  p_modt <- 2*pt(-abs(R[,"beta"]/sqrt(sq$var.post*R[,"vxx"])), R[,"df"] + sq$df.prior)
  # moderated HC3: scale HC3 SE by sqrt(var.post/s2) -- shrink the robust SE toward the prior
  se_hc3 <- R[,"beta"] / qt(R[,"p_hc3"]/2, R[,"df"], lower.tail = FALSE) ; se_hc3 <- abs(se_hc3)
  p_modhc3 <- 2*pt(-abs(R[,"beta"]/(se_hc3*sqrt(sq$var.post/R[,"s2"]))), R[,"df"] + sq$df.prior)
  P <- cbind(ols = R[,"p_ols"], modt = p_modt, hc3 = R[,"p_hc3"], modhc3 = p_modhc3)
  do.call(rbind, lapply(colnames(P), function(meth) {
    pv <- P[,meth]; q <- p.adjust(pv,"BH"); rej <- q <= 0.05
    tp <- sum(rej & is_da); fp <- sum(rej & !is_da); pn <- pv[!is_da]
    data.frame(n=n, error=kind, cell=if(null_cell)"null" else "alt", rep=rep_id, method=meth,
               fdr=if(tp+fp>0) fp/(tp+fp) else 0, tpr=if(sum(is_da)>0) tp/sum(is_da) else NA,
               fpr05=mean(pn<0.05), n_rej=tp+fp)
  }))
}
grid <- expand.grid(n=c(10,20,40,100), error=c("gauss","t3","contam","skew","hetero"), cell=c("null","alt"), stringsAsFactors=FALSE)
tasks <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) data.frame(grid[i,], rep=seq_len(if(grid$cell[i]=="null") n_rep_null else n_rep_alt))))
res <- mclapply(seq_len(nrow(tasks)), function(i){ tk<-tasks[i,]; set.seed(5e6+i)
  tryCatch(run_rep(tk$n,tk$error,tk$cell=="null",tk$rep), error=function(e) NULL)}, mc.cores=2)
res <- do.call(rbind,res); write.csv(res,"step1b_results.csv",row.names=FALSE)
cat("\n== FPR at 0.05 (null cells) ==\n")
a <- aggregate(fpr05 ~ method+error+n, subset(res,cell=="null"), mean); t<-reshape(a,idvar=c("error","n"),timevar="method",direction="wide"); names(t)<-sub("fpr05\\.","",names(t)); print(t[order(t$error,t$n),],digits=3,row.names=FALSE)
cat("\n== FDR / TPR at BH 0.05 (alt cells) ==\n")
a <- aggregate(cbind(fdr,tpr) ~ method+error+n, subset(res,cell=="alt"), mean); a$fdr<-round(a$fdr,3); a$tpr<-round(a$tpr,3)
t<-reshape(a,idvar=c("error","n"),timevar="method",direction="wide"); print(t[order(t$error,t$n),],row.names=FALSE)
