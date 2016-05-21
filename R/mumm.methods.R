
#' @export
print.mumm <- function(fit) {

  cat("Multiplicative mixel model fit by ML \n Formula:", deparse(fit$call$formula,width.cutoff = 500),
      "\n Data:", fit$call$data, "\n")

  cat("Log-likelihood at convergence:", -fit$objective, "\n")

  cat("Random effects: \n")
  #table with variance components
  table_sd = data.frame(Groups = names(fit$sigmas), Std.Dev. = fit$sigmas)
  print(table_sd, right = FALSE, row.names = FALSE)

  cat("Number of obs:", fit$nobs, "\n")

  cat("Fixed Effects: \n")
  print.default(fit$par_fix)

}

#' @export
ranef <- function(fit) UseMethod("ranef")

#' @export
ranef.mumm <- function(fit) {

  rand_ef = names(fit$sigmas)
  names_random_par = names(fit$par_rand)
  nlevels_par = fit$nlevels_par_rand
  par_rand_list = list()

  index = 1
  for(i in 1:length(fit$nlevels_par_rand)) {
    par_rand_list[[as.character(rand_ef[i])]] = fit$par_rand[index:(index+nlevels_par[[i]]-1)]
    index = index + nlevels_par[[i]]
  }

  print.default(par_rand_list)

}

#' @export
confint.mumm <- function(fit, names = "all", level = 0.95){

  name_vector = c(names(fit$par_fix),names(fit$sigmas))

  if(names[1] == "all"){
    index = 1:length(name_vector)
  } else {
    index = match(names,name_vector)
  }

  confints = c()

  for(i in 1:length(index)){
    profile = tmbprofile(fit$obj, index[i] , trace = FALSE)
    c = confint(profile, level = level)
    dimnames(c)[[1]] = names[i]
    confints = rbind(confints,c)

  }
  print.default(confints)
}

#' @export
lrt <- function(fit1,fit2) UseMethod("lrt")

#' @export
lrt.mumm <- function(fit1, fit2) {

  cat("Data:", fit1$call$data, "\n")
  cat("Models: \n Object:", deparse(fit1$call$formula),
      "\n ...1:", deparse(fit2$call$formula), "\n" )

  #table
  loglik1 = -fit1$objective
  loglik2 = -fit2$objective
  df = fit1$df-fit2$df
  chi = 2*(loglik1-loglik2)
  pvalue = 1-pchisq(chi,df  = df)
  table_sd = data.frame(Df = c(fit1$df,fit2$df), logLik = c(loglik1,loglik2), Chisq = c(NA,chi),
                        ChiDf = c(NA,df), pvalue = c(NA,pvalue))
  row.names(table_sd) <- c("Object","..1")
  print(as.matrix(table_sd), right = FALSE, digits = 4, na.print = "")


}
