#' @export
print.mumm <- function(fit) {

  cat("Multiplicative mixel model fit by ML \n Formula:")

  for(i in 1:length(deparse(fit$call$formula))) {
    if(i>1){cat("      ")}
    cat(deparse(fit$call$formula)[i],"\n")
  }


  cat("Data:", fit$call$data, "\n")

  cat("Log-likelihood at convergence:", -fit$objective, "\n")

  cat("Random effects: \n")
  #table with variance components
  table_sd = data.frame(Groups = names(fit$sigmas), Std.Dev. = fit$sigmas)
  print(table_sd, right = FALSE, row.names = FALSE)

  cat("Correlation:", fit$est_cor, "\n")

  cat("Number of obs:", fit$nobs, "\n")

  cat("Fixed Effects: \n")
  print.default(fit$par_fix)

}

#' Extract Random Effects
#'
#' A function to extract the estimated random effects from a model object of class mumm.
#'
#' @usage ## S3 method for class mumm
#' ranef(fit)
#'
#' @param fit an object of class mumm
#'
#' @return A named list with the estimated random effects, where each element in the list is
#' a numeric vector consisting of the estimated random effect coefficients for a random factor in the model.
#'
#' @examples
#' set.seed(100)
#' sigma_e <- 1.5
#' sigma_a <- 0.8
#' sigma_b <- 0.5
#' sigma_d <- 0.7
#' nu <- c(8.2, 6.2, 2.3, 10.4, 7.5, 1.9)
#'
#' nA <- 15
#' nP <- 6
#' nR <- 5
#'
#' a <- rnorm(nA, mean = 0, sd = sigma_a)
#' b <- rnorm(nA, mean = 0, sd = sigma_b)
#' d <- rnorm(nA*nP, mean = 0, sd = sigma_d)
#' e <- rnorm(nA*nP*nR, mean = 0, sd = sigma_e)
#'
#' Assessor <- factor(rep(seq(1,nA),each = (nP*nR)))
#' Product <- factor(rep(rep(seq(1,nP),each = nR), nA))
#' AssessorProduct <- (Assessor:Product)
#'
#' y <- nu[Product] + a[Assessor] + b[Assessor]*(nu[Product]-mean(nu)) + d[AssessorProduct] + e
#'
#' sim_data <- data.frame(y, Assessor, Product)
#'
#' fit <- mumm(y ~ 1 + Product + (1|Assessor) + (1|Assessor:Product) +
#'              mp(Assessor,Product) ,data = sim_data)
#'
#' ranef(fit)
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

#' Confidence Intervals for Model Parameters
#'
#' Computes confidence intervals for the fixed effect parameters and the variance components
#' in a fitted multiplicative mixed model.
#'
#' @usage ## S3 method for class mumm
#' confint(fit, parm = "all", level = 0.95)
#'
#' @param fit an object of class mumm.
#'
#' @param parm a vector of names specifying which parameters to compute confidence intervals for. If missing,
#' confidence intervals will be computed for all of the fixed effect paramters and all of the variance components.
#'
#' @param level the confidence level.
#'
#' @details The confidence intervals are computed by the profile likelihood method.
#'
#' @return A matrix with the first column showing the lower confidence limit and the second column showing the
#' upper limit for each parameter.
#'
#' @examples
#' set.seed(100)
#' sigma_e <- 1.5
#' sigma_a <- 0.8
#' sigma_b <- 0.5
#' sigma_d <- 0.7
#' nu <- c(8.2, 6.2, 2.3, 10.4, 7.5, 1.9)
#'
#' nA <- 15
#' nP <- 6
#' nR <- 5
#'
#' a <- rnorm(nA, mean = 0, sd = sigma_a)
#' b <- rnorm(nA, mean = 0, sd = sigma_b)
#' d <- rnorm(nA*nP, mean = 0, sd = sigma_d)
#' e <- rnorm(nA*nP*nR, mean = 0, sd = sigma_e)
#'
#' Assessor <- factor(rep(seq(1,nA),each = (nP*nR)))
#' Product <- factor(rep(rep(seq(1,nP),each = nR), nA))
#' AssessorProduct <- (Assessor:Product)
#'
#' y <- nu[Product] + a[Assessor] + b[Assessor]*(nu[Product]-mean(nu)) + d[AssessorProduct] + e
#'
#' sim_data <- data.frame(y, Assessor, Product)
#'
#' fit <- mumm(y ~ 1 + Product + (1|Assessor) + (1|Assessor:Product) +
#'              mp(Assessor,Product) ,data = sim_data)
#'
#' confint(fit)
#' confint(fit, parm = c('Product3', 'Assessor', 'mp Assessor:Product'), level = 0.90)
#'
#' @export
confint.mumm <- function(fit, parm = "all", level = 0.95){

  name_vector = c(names(fit$par_fix),names(fit$sigmas))

  if(parm[1] == "all"){
    index = 1:length(name_vector)
  } else {
    index = match(parm,name_vector)
    name_vector <- name_vector[index]
  }

  confints = c()

  n_parfix = length(fit$par_fix)

  for(i in 1:length(index)){
    profile = tmbprofile(fit$obj, index[i] , trace = FALSE)
    c = confint(profile, level = level)
    dimnames(c)[[1]] = parm[i]

    #If parameter is a variance component
    if(index[i]>n_parfix){
      confints = rbind(confints,exp(c))
    }
    else{
      confints = rbind(confints,c)
    }

  }
  dimnames(confints)[[1]] <- name_vector
  print.default(confints)
}

#' Likelihood Ratio Test
#'
#' A function to perform the likelihood ratio test for testing two nested models against each other.
#'
#' @param fit1 a fitted model object of class mumm.
#'
#' @param fit2 a fitted model object of class mumm or a model object returned by \code{lm} or \code{lmer}.
#'
#' @details Performs the likelihood ratio test for testing two nested models against each other. The model in
#' \code{fit2} should be nested within the model in \code{fit1}.
#'
#' @return A matrix with the likelihood ratio test statistic and the corresponding p-value.
#'
#' @examples
#' set.seed(100)
#' sigma_e <- 1.5
#' sigma_a <- 0.8
#' sigma_b <- 0.5
#' sigma_d <- 0.7
#' nu <- c(8.2, 6.2, 2.3, 10.4, 7.5, 1.9)
#'
#' nA <- 15
#' nP <- 6
#' nR <- 5
#'
#' a <- rnorm(nA, mean = 0, sd = sigma_a)
#' b <- rnorm(nA, mean = 0, sd = sigma_b)
#' d <- rnorm(nA*nP, mean = 0, sd = sigma_d)
#' e <- rnorm(nA*nP*nR, mean = 0, sd = sigma_e)
#'
#' Assessor <- factor(rep(seq(1,nA),each = (nP*nR)))
#' Product <- factor(rep(rep(seq(1,nP),each = nR), nA))
#' AssessorProduct <- (Assessor:Product)
#'
#' y <- nu[Product] + a[Assessor] + b[Assessor]*(nu[Product]-mean(nu)) + d[AssessorProduct] + e
#'
#' sim_data <- data.frame(y, Assessor, Product)
#'
#' fit <- mumm(y ~ 1 + Product + (1|Assessor) + (1|Assessor:Product) +
#'              mp(Assessor,Product) ,data = sim_data)
#'
#' fit2 <- mumm(y ~ 1 + Product + (1|Assessor) + mp(Assessor,Product) ,data = sim_data)
#' lrt(fit,fit2)
#'
#' @export
lrt <- function(fit1,fit2) UseMethod("lrt")

#' @export
lrt.mumm <- function(fit1, fit2) {

  cat("Data:", fit1$call$data, "\n")
  cat("Models: \n Object:")

  for(i in 1:length(deparse(fit1$call$formula))) {
    if(i>1){cat("    ")}
    cat(deparse(fit1$call$formula)[i],"\n")
  }


  cat("...1   :")

  if(class(fit2)=="mumm") {

    loglik2 = -fit2$objective
    df2 = fit2$df

    for(i in 1:length(deparse(fit2$call$formula))) {
      if(i>1){cat("    ")}
      cat(deparse(fit2$call$formula)[i],"\n")
    }
  } else {
    loglikelihood = logLik(fit2, REML = F)
    loglik2 = loglikelihood[1]
    df2 = attr(loglikelihood,"df")

    if(class(fit2)=="lm"){

      for(i in 1:length(deparse(fit2$call$formula))) {
        if(i>1){cat("    ")}
        cat(deparse(fit2$call$formula)[i],"\n")
      }

    } else {

      for(i in 1:length(deparse(fit2@call$formula))) {
        if(i>1){cat("    ")}
        cat(deparse( fit2@call$formula)[i],"\n")
      }
    }

  }

  #table
  loglik1 = -fit1$objective
  df1 = fit1$df
  df = df1-df2
  chi = 2*(loglik1-loglik2)
  pvalue = 1-pchisq(chi,df  = df)
  table_sd = data.frame(Df = c(df1,df2), logLik = c(loglik1,loglik2), Chisq = c(NA,chi),
                        ChiDf = c(NA,df), pvalue = c(NA,pvalue))
  row.names(table_sd) <- c("Object","..1")
  print(as.matrix(table_sd), right = FALSE, digits = 4, na.print = "")


}

