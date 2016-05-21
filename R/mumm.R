#' Fit Multiplicative Mixed Models with TMB
#'
#' Fit a multiplicative mixed-effects model to data with use of the Template Model Builder.
#'
#' @useDynLib mumm
#' @importFrom lme4 subbars findbars mkReTrms nobars
#' @importFrom TMB MakeADFun sdreport tmbprofile
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

  #The number of multiplicative terms in the model
  n_mult = length(attr(terms_mult,"term.labels"))

  #Finding the random and fixed effects that are part of the mulitiplicative terms
  randomef = vector(length = 0)
  fixedef = vector(length = 0)
  for(i in 1:n_mult) {
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

    if(sum((attr(terms_fix,"term.labels")==fixedef[i]))==0) {
      stop(sprintf("%s in multipliactive term is not a fixed effect",fixedef[i]))
      }

    remove = remove + (attr(terms_fix,"term.labels")==fixedef[i])
  }
  NUMremove = 1:length(attr(terms_fix,"term.labels"));
  NUMremove1 = NUMremove[as.logical(remove)];  #part of scaling
  NUMremove2 = NUMremove[!as.logical(remove)]; #not part of scaling


  #------------Building fixed effect design matrices, X and Xnu -------------------

  #If all of the fixed effects are part of the multiplicative terms
  if (length(NUMremove2) == 0) {

    terms_fix_ed = NULL  #fixed effects that are not part of the multiplicative terms
    X = matrix(1, nrow = nrow(data), ncol = 1) #intercept
    colnames(X) = "(Intercept)"
    Xnu = model.matrix(terms_fix,data = data)
    t = as.data.frame(table(attr(Xnu,"assign")))
    Xnu = Xnu[,-1]    #Removing the intercept
    sizenu = t[t[["Var1"]]%in%NUMremove1,]$Freq

  } else {

    #fixed effects that are not part of the multiplicative terms
    #terms_fix_ed = drop.terms(terms_fix,NUMremove1,keep.response = TRUE)
    #formula_fix_ed = formula(terms_fix_ed)
    #X = model.matrix(formula_fix_ed,data = data) #with intercept

    terms_fixInMult = drop.terms(terms_fix,NUMremove2,keep.response = TRUE)
    formula_fixInMult = formula(terms_fixInMult)

    #Xnu = model.matrix(terms_fixInMult,data = data)
    #Xnu = Xnu[,-1]    #Removing the intercept

    Xbig = model.matrix(terms_fix,data = data)

    X = Xbig[,1, drop = F]
    Xindex = as.logical(attr(Xbig,"assign")%in%NUMremove2 + (attr(Xbig,"assign")==0))
    X = Xbig[,Xindex, drop= F]
    Xnu = Xbig[,!Xindex, drop = F]

    t = as.data.frame(table(attr(Xbig,"assign")))

    #The number of parameters to be estimated for each fixed effect in the multiplicative term
    sizenu = t[t[["Var1"]]%in%NUMremove1,]$Freq

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

  #number of random effects in the model (excluding random regression coef.)
  n_rand = length(rterms$cnms)

  data_named = data

  colnames(data_named)[colnames(data)==formula[[2]]] <- "y"

  attach(data_named)

  #Want to include mp interactions in data_named
  for(i in 1:length(fixedef)){
    if(regexpr(":",fixedef[i]) != -1){
      data_named[fixedef[i]] = eval(parse(text=fixedef[i]))
    }
  }

  detach(data_named)

  nlevelsf = sapply(data_named[fixedef],nlevels);
  nlevelsr = sapply(data_named[randomef],nlevels);

  ffac = data.matrix(data_named[fixedef]);
  ffac = matrix(ffac,dimnames = NULL, ncol = length(fixedef))
  rfac = data.matrix(data_named[randomef]);
  rfac = matrix(rfac,dimnames = NULL, ncol = length(randomef))

  dataTMB = list(y = data_named[['y']],ffac = ffac, rfac = rfac,
                 X = X, Z = Z, Xnu = Xnu,  npar = npar, nlevelsf = nlevelsf, nlevelsr = nlevelsr, sizenu = sizenu)


  parameters = list(
    beta  = rep(0,ncol(dataTMB$X)),
    a = rep(0,ncol(dataTMB$Z)),
    b  = rep(0,sum(nlevelsr)),
    nu  = rep(0,ncol(Xnu)),
    log_sigma_a     = rep(0,n_rand),
    log_sigma_b     = rep(0,n_mult),
    log_sigma       = 0
  )


  obj <- MakeADFun(
    data=dataTMB,
    parameters= parameters,
    random = c("a","b"),
    DLL    = "mumm",
    silent = TRUE
  )

  opt = nlminb(obj$par,obj$fn,obj$gr);

  sdr = sdreport(obj)


  #---------------------------- Finalizing output ---------------------------------------

  ## The estimated fixed effect coefficients
  par_fix = opt$par[1:(length(opt$par)-(n_rand+n_mult+1))]


  names_fixed_ef = c(colnames(X),colnames(Xnu))

  names(par_fix) = names_fixed_ef



  ##The estimated variance components
  sigmas = exp(opt$par[(length(opt$par)-(n_rand+n_mult)):length(opt$par)])

  ##Creating a vector with names for the variance components
  names_random_ef = c(names(rterms$flist),paste0("mp ",randomef,":",fixedef),"Residual")
  names(sigmas) = names_random_ef


  ##Random effects
  par_rand = sdr$par.random
  ## Creating a vector with names for the random parameters
  mp_rdata = data_named[randomef]

  if(is.null(findbars(randform[[3]]))){
    names_random_par = sapply(mp_rdata,levels)
  } else {
    names_random_par = c(rterms$Zt@Dimnames[[1]],sapply(mp_rdata,levels))
  }

  names(par_rand) = names_random_par

  nlevels_par_rand = c(sapply(rterms$flist,nlevels),sapply(mp_rdata,nlevels))



  res = list(par = opt$par, objective = opt$objective, convergence = opt$convergence,
             iterations = opt$iterations, evaluations = opt$evaluations, convmessage = opt$message,
             par_fix = par_fix, sigmas = sigmas , par_rand = par_rand, nlevels_par_rand =  nlevels_par_rand,
             call = match.call(), nobs = nrow(data), df = length(opt$par), sdreport = sdr, obj = obj)

  class(res) <- "mumm"
  return(res)


