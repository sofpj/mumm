
significance <- function(fit) {

  #----------------------- Testing the random effecs:..........................

  formula1 <- as.formula(fit$call$formula)
  terms1 <- terms(formula1)
  mp_string = 'mp\\('

  rand_terms = findbars(formula1)

  # Helping function for reformulate and update
  parens <- function(x) paste0(".-(",x,")")
  for_reform <- function(x) paste0(".-",x,"")
  for_reform_p <- function(x) paste0(".+",x,"")

  pvalue= c()
  count = 0
  #Going through all of the random effects (not mp terms)
  if (!is.null(rand_terms)) {
    for(i in 1:length(rand_terms)){

      reform = reformulate(parens(rand_terms[i]))
      formula2 = update(formula1,reform)
      fit2 <- mumm(formula2,data = fit$data)

      loglik1 = -fit$objective
      loglik2 = -fit2$objective
      df = fit1$df-fit2$df
      chi = 2*(loglik1-loglik2)
      count = count + 1
      pvalue[count] = 1-pchisq(chi,df  = df)

    }
  }


  fixedform <-nobars(formula1)
  terms_fix <- terms(fixedform)

  #Going through all of the multiplicative effects:
  if(length(fit$index_num) == 1) {
    term_mp = attr(terms_fix,"term.labels")[fit$index_num]

    reform = reformulate(for_reform(term_mp))
    formula2 = update(formula1,reform)
    if (!is.null(rand_terms)) {
      fit2 <- lmer(formula2,data = fit$data, REML = F)
    } else {
      fit2 <- lm(formula2,data = fit$data)
    }
    loglikelihood = logLik(fit2)
    loglik1 = -fit$objective
    loglik2 = loglikelihood[1]
    df = fit1$df-attr(loglikelihood,"df")
    chi = 2*(loglik1-loglik2)
    count = count + 1
    pvalue[count] = 1-pchisq(chi,df  = df)

  } else {

    for(i in 1:length(fit$index_num)) {
      term_mp = attr(terms_fix,"term.labels")[fit$index_num[i]]
      reform = reformulate(for_reform(term_mp))
      formula2 = update(formula1,reform)

      fit2 = mumm(formula2,data = fit$data)

      loglik1 = -fit$objective
      loglik2 = -fit2$objective
      df = fit1$df-fit2$df
      chi = 2*(loglik1-loglik2)
      count = count + 1
      pvalue[count] = 1-pchisq(chi,df  = df)

    }

  }

  #----------------------- Testing the fixed effecs:..........................

  num_fixed = length(attr(terms_fix,"term.labels"))-length(fit$index_num)

  terms_only_fix = drop.terms(terms_fix,fit$index_num)

  all_terms_vstring = attr(terms1,"term.labels")

  for(i in 1:num_fixed){
    term_fix =attr(terms_only_fix,"term.labels")[i]

    #NOTE: maybe a problem with the interactions!

    index_remove = grep(term_fix,all_terms_vstring)
    remove_terms_vstring = all_terms_vstring[index_remove]
    formula_small = formula1

    for(j in 1:length(index_remove)) {
    reform = reformulate(for_reform(remove_terms_vstring[j]))
    formula_small = update(formula_small,reform)
    }

    reform = reformulate(for_reform_p(term_fix))
    formula_big = update(formula_small,reform)

    small_terms_vstring = attr(terms(formula_small),"term.labels")

    #If there are mp terms in the formulas
    if(sum(grep(mp_string,small_terms_vstring))>0){

      fit_small = mumm(formula_small, data = fit$data)
      fit_big = mumm(formula_big, data = fit$data)

      loglik1 = -fit_big$objective
      loglik2 = -fit_small$objective
      df = fit1$df-fit2$df
      chi = 2*(loglik1-loglik2)
      count = count + 1
      pvalue[count] = 1-pchisq(chi,df  = df)

    #If there are NO mp terms in the formulas
    } else {

      if(is.null(findbars(formula_small))) {
        fit_small = lm(formula_small, data = fit$data)
        fit_big = lm(formula_big, data = fit$data)
      } else {
        fit_small = lmer(formula_small, data = fit$data)
        fit_big = lmer(formula_big, data = fit$data)
      }

      loglikelihood_big = logLik(fit_big)
      loglikelihood_small = logLik(fit_small)
      loglik1 = loglikelihood_big[1]
      loglik2 = loglikelihood_small[1]
      df = attr(loglikelihood_big,"df")-attr(loglikelihood_small,"df")
      chi = 2*(loglik1-loglik2)
      count = count + 1
      pvalue[count] = 1-pchisq(chi,df  = df)


    }


  }


}
