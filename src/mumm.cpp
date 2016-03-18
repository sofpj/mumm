#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Data Section*/
  DATA_MATRIX(X);
  DATA_MATRIX(Xnu);
  DATA_SPARSE_MATRIX(Z);
  DATA_MATRIX(rfac); /* måske DATA_MATRIX og så konverter til Integer når de skal bruges som index */
  DATA_MATRIX(ffac);
  DATA_VECTOR(y);
  DATA_IVECTOR(npar); /* The number of levels in all of the random effects belonging to a */
  DATA_IVECTOR(nlevelsf); /* The number of levels in all of the fixed effects in the multiplicative terms (nu) */
  DATA_IVECTOR(nlevelsr); /* The number of levels in all of the random effects in the multiplicative terms* (b) */


  using CppAD::Integer;

  /* Parameter Section */
  PARAMETER_VECTOR(beta); /*related to X */
  PARAMETER_VECTOR(a);      /*related to Z */
  PARAMETER_VECTOR(b);      /*random scalling coef*/
  PARAMETER_VECTOR(nu);

  PARAMETER_VECTOR(log_sigma_a);
  PARAMETER_VECTOR(log_sigma_b);
  PARAMETER(log_sigma);

  Type sigma = exp(log_sigma);
  vector<Type> sigma_a = exp(log_sigma_a);
  vector<Type> sigma_b = exp(log_sigma_b);


  /* Including the estimates and standard errors on natural scale */
  ADREPORT(sigma);
  ADREPORT(sigma_a);
  ADREPORT(sigma_b);


  Type nll = 0;
  int nobs = y.size();

  /* The linear part of the model */
  vector<Type> linmodel = X*beta + Z*a + Xnu*nu;


  /* The multiplicative part of the model */

  vector<Type> mult(nobs);      /*"nu"*/
  int jumpf = 0;
  int jumpr = 0;
  int indx = 0;

  /*Going through all of the multiplicative terms*/

  /*The first multiplicative term j = 0 */

  vector<Type> nuj(nlevelsf[0]);
  vector<Type> bj(nlevelsr[0]);
  indx = 0;

  /*Build nu for the j'th multiplicative term*/
  for (int l = 0; l<nlevelsf[0]; l++) {
    indx = l+jumpf;
    nuj[l] = nu[indx];
  }

  indx = 0;
  /*Build b for the j'th multiplicative term*/
  for (int l=0; l<nlevelsr[0]; l++) {
    indx = l + jumpr;
    bj[l] = b[indx];
  }

  vector<Type> nu2 = nuj - nuj.sum()/nuj.size();

  /*Going through all of the observations */
  for (int i = 0; i<nobs; i++){

    mult[i] += bj[Integer(rfac(i,0))-1]*nu2[Integer(ffac(i,0))-1];

  }
  jumpf += nlevelsf[0];
  jumpr += nlevelsr[0];


  /* The rest of the multiplicative terms */

  for (int j = 1; j<nlevelsf.size(); j++){
    vector<Type> nuj(nlevelsf[j]);
    vector<Type> bj(nlevelsr[j]);
    indx = 0;

    /*Build nu for the j'th multiplicative term*/
    nuj[0] = 0;   /*to overcome identifiability issues*/
    for (int l = 1; l<nlevelsf[j]; l++) {
      indx = l+jumpf-1;
      nuj[l] = nu[indx];
    }

    indx = 0;
    /*Build b for the j'th multiplicative term*/
    for (int l=0; l<nlevelsr[j]; l++) {
      indx = l + jumpr;
      bj[l] = b[indx];
    }

    vector<Type> nu2 = nuj - nuj.sum()/nuj.size();

    /*Going through all of the observations */
    for (int i = 0; i<nobs; i++){

      mult[i] += bj[Integer(rfac(i,j))-1]*nu2[Integer(ffac(i,j))-1];

    }
    jumpf += nlevelsf[j];
    jumpr += nlevelsr[j];
  }



  /* The negative joint log-likelihood function */

  for(int i=0; i<nobs; i++){

    nll -= dnorm(y[i],
                 linmodel[i]+
                   mult[i],sigma,true);

  }

  int index_count = 0;
  int index = 0;

  /* Going through all of the random effects in a */
  for(int i=0; i<sigma_a.size(); i++){

    for(int j=0; j<npar[i]; j++){
      index = j+index_count;
      nll -= dnorm(a[index], Type(0), sigma_a[i], true);
    }

    index_count += npar[i];
  }


  index_count = 0;

  /* Going through all of the random effects in b */
  for(int i=0; i<sigma_b.size(); i++){

    for(int j=0; j<nlevelsr[i]; j++){
      index = j+index_count;
      nll -= dnorm(b[index], Type(0), sigma_b[i], true);
    }


    index_count += nlevelsr[i];
  }

  return nll;
}

