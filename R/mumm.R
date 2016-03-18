#' Fit Multiplicative Mixed Models with TMB
#'
#' Fit a multiplicative mixed-effects model to data with use of the Template Model Builder.
#'
#' @useDynLib mumm
#' @importFrom lme4 subbars findbars mkReTrms nobars
#' @importFrom TMB MakeADFun sdreport
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix t
#' @export
mumm <- function(formula, data, report = FALSE) {

  #Checiking input:


  #-------------Building fixed effect design matrix, X-----------------------

  fixedform <- formula   #formula_obj
  fixedform[[3]]<-nobars(fixedform[[3]])
  terms_fix_mult = attr(terms(fixedform),"term.labels")

  #Seperating the fixed effect terms and the multiplicative terms
  TFind = rep(FALSE,length(terms_fix_mult))

  for(i in 1:length(terms_fix_mult)) {
    TFind[i] = substr(terms_fix_mult[i],1,3) == "mp(";
  }
  NUMind = 1:length(terms_fix_mult);
  NUMind1 = NUMind[TFind]  #index for mp terms
  NUMind2 = NUMind[!TFind] #index for fixed terms

  terms_fix = drop.terms(terms(fixedform),NUMind1,keep.response = TRUE)
  formula_fix = formula(terms_fix)
  terms_mult = drop.terms(terms(fixedform),NUMind2,keep.response = TRUE)
  formula_mult = formula(terms_mult)

  #Finding the random and fixed effects that are part of the mulitiplicative terms
  randomef = vector(length = 0)
  fixedef = vector(length = 0)
  for(i in 1:length(attr(terms_mult,"term.labels"))) {
    mult_term = attr(terms_mult,"term.labels")[i];
    pos = stringr::str_locate_all(pattern = ",",mult_term)
    randomef[i] = substr(mult_term,4,pos[[1]][1]-1);
    fixedef[i] = substr(mult_term,pos[[1]][1]+1,nchar(mult_term)-1);

  }
  #Removing potential blank space from the strings
  randomef = stringr::str_trim(randomef)
  fixedef = stringr::str_trim(fixedef)


  #Seperating the fixed effects that is and isn't part of the multiplicative term
  remove = rep(FALSE,length(attr(terms_fix,"term.labels")))
  for(i in 1:length(fixedef)){
    remove = remove + (attr(terms_fix,"term.labels")==fixedef[i])
  }
  NUMremove = 1:length(attr(terms_fix,"term.labels"));
  NUMremove1 = NUMremove[as.logical(remove)];  #part of scaling
  NUMremove2 = NUMremove[!as.logical(remove)]; #not part of scaling



  #------------Building fixed effect design matrices, X and Xnu -------------------

  #If all of the fixed effects are part of the multiplicative terms
  if (length(attr(terms_fix,"term.labels")) == length(attr(terms_mult,"term.labels"))) {

    X = matrix(0, nrow = nrow(data), ncol = 0)
    formula_fixInMult = formula(terms_fix)
    terms_fixInMult = terms(formula_fixInMult)
    attr(terms_fixInMult,"intercept") <- 0
    Xnu = model.matrix(terms_fixInMult,data = data)

  } else{

    terms_fix_ed = drop.terms(terms_fix,NUMremove1,keep.response = TRUE)
    formula_fix_ed = formula(terms_fix_ed)

    terms_fixInMult = drop.terms(terms_fix,NUMremove2,keep.response = TRUE)
    formula_fixInMult = formula(terms_fixInMult)

    X = model.matrix(formula_fix_ed,data = data)
    X = X[,-1]    #Removing the intercept

    attr(terms_fixInMult,"intercept") <- 0
    Xnu = model.matrix(terms_fixInMult,data = data)

  }

  #Building random effect design matrix, Z
  randform <- formula

  if (is.null(findbars(randform[[3]]))) {

    Z = matrix(0, nrow = nrow(data), ncol = 0)
    Z = as(Z,"dgTMatrix")
    rterms = NULL
    npar= 0

  } else {

    rterms = mkReTrms(findbars(randform[[3]]),data)
    Z = t(rterms$Zt)
    npar = sapply(rterms$flist,nlevels)
  }


  data_named = data

  colnames(data_named)[colnames(data)==formula[[2]]] <- "y"

  nlevelsf = sapply(data_named[fixedef],nlevels);
  nlevelsr = sapply(data_named[randomef],nlevels);

  ffac = data.matrix(data_named[fixedef]);
  ffac = matrix(ffac,dimnames = NULL, ncol = length(fixedef))
  rfac = data.matrix(data_named[randomef]);
  rfac = matrix(rfac,dimnames = NULL, ncol = length(randomef))

  dataTMB = list(y = data_named[['y']],ffac = ffac, rfac = rfac,
                 X = X, Z = Z, Xnu = Xnu,  npar = npar, nlevelsf = nlevelsf, nlevelsr = nlevelsr)


  parameters = list(
    beta  = rep(0,ncol(dataTMB$X)),
    a = rep(0,ncol(dataTMB$Z)),
    b  = rep(0,sum(nlevelsr)),
    nu  = rep(0,ncol(Xnu)),
    log_sigma_b     = rep(0,ncol(data_named[randomef])),
    log_sigma_a     = rep(0,length(rterms$cnms)),
    log_sigma       = 0
  )


  obj <- MakeADFun(
    data=dataTMB,
    parameters= parameters,
    random = c("a","b"),
    DLL    = "mumm",
    silent = FALSE
  )

  opt = nlminb(obj$par,obj$fn,obj$gr);

  if(report == TRUE){
    sdr = sdreport(obj);
    return(list(opt = opt, sdr = sdr, obj = obj))
  } else {
    return(list(opt = opt, obj = obj))
  }

}
