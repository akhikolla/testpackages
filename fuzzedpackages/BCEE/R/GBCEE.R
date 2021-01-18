GBCEE <- function(X, Y, U, omega, niter = 5000, family.X = "gaussian", family.Y = "gaussian", 
                   X1 = 1, X0 = 0, priorX = NA, priorY = NA, maxsize = NA, OR = 20, truncation = c(0.01, 0.99),
                   var.comp = "asymptotic", B = 200)
{
 # Error checks
 if(!is.numeric(X)){stop("X should be a numeric variable")};
 if(!is.numeric(Y)){stop("Y should be a numeric variable")};
 if(!is.matrix(U)){stop("U should be a matrix")};
 if(!is.numeric(U)){stop("U should contain numeric variables")};
 if(!is.numeric(omega)){stop("omega should be a numeric value")};
 if(length(omega) > 1){stop("omega should contain a single value")};
 if(omega <0){stop("omega should be a value between 0 and Inf")};
 if(!is.numeric(niter)){stop("niter should be a numeric value")};
 if(length(niter) > 1){stop("niter should contain a single value")};
 if(niter < 1){stop("niter should be >= 1")};
 if(niter%%1 != 0){niter = niter - niter%%1;
                   warning("niter was not an integer; it has been rounded down.");};
 if(!(family.X %in%c("gaussian", "binomial"))){stop("family.X should either be 'gaussian' or 'binomial'")};
 if(!(family.Y %in%c("gaussian", "binomial"))){stop("family.Y should either be 'gaussian' or 'binomial'")};
 if(!is.numeric(X1)){stop("X1 should be a numeric value")};
 if(length(X1) > 1){stop("X1 should be a single value")};
 if(!is.numeric(X0)){stop("X0 should be a numeric value")};
 if(length(X0) > 1){stop("X0 should be a single value")};
 if(!is.na(priorX[1])){
  if(!is.numeric(priorX)){stop("priorX should either be NA or a numeric vector of length ncol(U)")};
  if(length(priorX) != ncol(U)){stop("priorX should either be NA or a numeric vector of length ncol(U)")};
 }
 if(!is.na(priorY[1])){
  if(!is.numeric(priorY)){stop("priorY should either be NA or a numeric vector of length ncol(U)")};
  if(length(priorY) != ncol(U)){stop("priorY should either be NA or a numeric vector of length ncol(U)")};
 }
 if(!is.na(maxsize)){
  if(!is.numeric(maxsize)){stop("maxsize should either be NA a single numeric value")};
  if(length(maxsize) > 1){stop("maxsize should either be NA a single numeric value")};
  if(maxsize < 0){stop("maxsize should be >= 0")};
  if(maxsize%%1 != 0){maxsize = maxsize - maxsize%%1;
                      warning("maxsize was not an integer; it has been rounded down.")};
 }
 if(!is.numeric(OR)){stop("OR should be a numeric value")};
 if(length(OR) > 1){stop("OR should be a single value")};
 if(OR < 1){stop("OR should be >=1")};
 if(!is.numeric(truncation)){stop("truncation should be a numeric vector of length 2")};
 if(length(truncation) != 2){stop("truncation should be a numeric vector of length 2")};
 if(truncation[1] < 0 | truncation[1] > 1 | truncation[2] < 0 | truncation[2] > 1){
  stop("truncation should have values between 0 and 1")};
 if(!(var.comp %in% c("asymptotic", "bootstrap"))){stop("var.comp should either be 'asymptotic' or 'bootstrap'")};
 if(var.comp == "bootstrap"){
  if(!is.numeric(B)){stop("B should be a numeric value")};
  if(length(B) > 1){stop("B should be a single value")};
  if(B < 2){stop("At least 2 bootstrap samples are required to compute a standard deviation (but much more is desirable)")};
  if(B%%1 !=0){B = B - B%%1; warning("B was not an integer; it has been rounded down.");};
 }

 sy = sd(Y);
 if(family.Y == "binomial") sy = 1;
 su = apply(as.matrix(U), 2, sd);
 n_cov = ncol(as.matrix(U)); #Number of covariate, potential confounders
 n = length(Y); #Sample size

 norm.sample = rnorm(1000);

 if(is.na(priorX[1])) priorX = rep(0.5, n_cov)
 if(is.na(priorY[1])) priorY = rep(0.5, n_cov)

 #Define maximal model size if not supplied by the user
 if(is.na(maxsize))
 {
  maxsize = floor(n_cov); 
 }

 resultsX = summary(regsubsets(y = X, x = as.matrix(U), nbest = 1, really.big = T, nvmax = maxsize, method = "forward")); #use linear regression for first screening of models
 alpha_X = resultsX$which[which.min(resultsX$bic), -1]; #Inital exposure model is the best fitting linear regression
	
 ###Obtaining BIC for X
 if(sum(alpha_X) == 0) 
 {
  model_X = glm(X~1, family = family.X);
 } else
 {
  model_X = glm(X~U[,alpha_X == 1], family = family.X);
 }
 bic_init = bic_x = BIC(model_X);

 if(family.X == "binomial")
 {
  if(length(table(X)) > 2) stop("X should only have two levels");
  X = as.numeric(X==X1);
 }

 models_X = matrix(0, nrow = niter, ncol = n_cov); #objects that will contain results
 tested_models_X = matrix(-1, nrow = niter, ncol = 1); #objects that will contain information for models already tested in character form
 tested_models_X2 = matrix(NA, nrow = niter, ncol = n_cov); #objects that will contain information for models already tested in numeric form
 pX_tested = matrix(NA, nrow = niter, ncol = 1); #objects that will contain information for models already tested
 alpha_X = as.numeric(alpha_X);
 char_alpha_X = paste(alpha_X, collapse = "");

 tested_models_X[1] = char_alpha_X; #Recording info about first model
 tested_models_X2[1,] = alpha_X; #Recording info about first model
 pX = pX_tested[1] = exp(-bic_x/2 + bic_init/2)*prod((2*priorX)**alpha_X*(2*(1-priorX))**(1-alpha_X)); #Recording info about first model

 for(i in 2:niter)
 {
  alpha_X_0 = alpha_X; #Current alpha_X
  if(sum(alpha_X) < maxsize){
   alpha_X_1 = alpha_X; temp = sample(n_cov, 1); alpha_X_1[temp] = (alpha_X_1[temp] + 1)%%2; 
   #Candidate alpha_X
   #According to MC3, we select a model in the neighborhood of the candidate model,
   #that is, a model with one more variable or one variable fewer. Each model
   #have the same probability. Here, I randomly select one of the n_cov potential
   #confounder and either remove it or add it according to whether the variable
   #was already in the model or not, respectively.
  } else {
   alpha_X_1 = alpha_X; temp = sample((1:n_cov)[alpha_X == 1], 1); alpha_X_1[temp] = (alpha_X_1[temp] + 1)%%2;
   # if model is already at maxsize, sample a variable already in the model for removal
  } 
  pX_0 = pX;

  char_alpha_X_1 = paste(alpha_X_1, collapse = "");
  tested = which(tested_models_X == char_alpha_X_1);
  if(length(tested) == 1) #The candidate model was already tested
  {
   pX_1 = pX_tested[tested];
  } else #The candidate model was never tested
  {
   ###Obtaining elements to calculate B10 for X
   if(sum(alpha_X_1) == 0) 
   {
    model_X = glm(X~1, family = family.X);
   }
   else{
    model_X = glm(X~U[,alpha_X_1 == 1], family = family.X);
   }
   bic_x_1 = BIC(model_X);
   tested_models_X[i] = char_alpha_X_1; #Recording info about candidate model
   tested_models_X2[i,] = alpha_X_1;
   pX_1 = pX_tested[i] = exp(-bic_x_1/2 + bic_init/2)*prod((2*priorX)**alpha_X_1*(2*(1-priorX))**(1-alpha_X_1)); #Recording info about candidate model
  }
  ratio = pX_1/pX_0; #Marginal likelihood x prior probability
  if(sample(c(1,0), size = 1, prob = c(min(1, ratio), 1 - min(1, ratio))))
  {
   alpha_X = alpha_X_1;
   pX = pX_1;
  }
 } #end of loop
 tested_models_X2 = na.omit(tested_models_X2);
 pX_tested = as.numeric(na.omit(pX_tested));
 pX_tested[pX_tested/max(pX_tested) < 1/OR] = 0;
 pX_tested = pX_tested/sum(pX_tested);
 alpha_X = colSums(tested_models_X2*pX_tested);
 models.X = cbind(tested_models_X2,pX_tested);
 colnames(models.X) = c(1:n_cov, "PostProb");

 resultsY = summary(regsubsets(y = Y, x = cbind(X, U), nbest = 1, really.big = T, nvmax = maxsize, force.in = 1, method = "forward")); #use linear regression for first screening of models
 alpha_Y = resultsY$which[which.min(resultsY$bic), -c(1,2)];  #Inital outcome model is the best fitting linear regression

 if(sum(alpha_Y) == 0)
 {
  model_Y0 = glm(Y~X, family = family.Y);
  model_X0 = glm(X~1, family = family.X);
 } else
 {
  model_Y0 = glm(Y~X + U[ , alpha_Y == 1], family = family.Y);
  model_X0 = glm(X~1 + U[ , alpha_Y == 1], family = family.X);
 }

 bic_init = bic_y = BIC(model_Y0);

 if(family.Y == "binomial")
 {
  betas_t = matrix(NA, nrow = niter, ncol = 2); #exposure betas estimated with TMLE or AIPTW
  vars_t = matrix(NA, nrow = niter, ncol = 2); #exposure beta variances with TMLE or AIPTW
 } else #family.Y == "gaussian"
 {
  betas_t = matrix(NA, nrow = niter, ncol = 1); #exposure betas estimated with TMLE or AIPTW
  vars_t = matrix(NA, nrow = niter, ncol = 1); #exposure beta variances with TMLE or AIPTW
 }
 models_Y = matrix(NA, nrow = niter, ncol = n_cov); #Models in numeric form
 tested_models_Y = matrix(-1, nrow = niter, ncol = 1); #Models in alphanumeric form
 pY_tested = matrix(NA, nrow = niter, ncol = 1); #BIC
 pY_tested1 = matrix(NA, nrow = niter, ncol = n_cov); #betas of the U covariates of tested models
 pY_tested2 = matrix(NA, nrow = niter, ncol = n_cov); #sd of the betas of the U covariates
 pY_tested3 = matrix(NA, nrow = niter, ncol = 1); #ratio of posterior probability relative to initial model

 #Recording infos about initial model:
 if(family.Y == "gaussian" & family.X == "gaussian")
 {
  #TMLE
  B0 = coef(model_Y0)[2]; # Initial estimate of slope
  Q0 = predict(model_Y0); # Y hat
  r0 = Q0 - B0*X; 
  predX0 = predict(model_X0); # X hat
  #Perform truncation ...
  Hg = X - predX0; # Clever covariate
  epsilon = coef(lm(Y~-1+Hg, offset = Q0)); # Fluctuation of estimate
  B1 = betas_t[1] = B0 + epsilon; # Final TMLE estimate of slope
		
  if(var.comp == "asymptotic")
  {
   # Compute empirical influence curve for variance estimation
   r1 = r0 - epsilon*predX0; 
   Q1 = B1*X + r1; # Updated Y hat
   k = -mean(Hg*X); 
   D = Hg*(Y - Q1); # Score
   IC = 1/k*D; # Influence curve
   vars_t[1] = mean(IC**2)/n;
  } else #var.comp == "bootstrap"
  {
   bootf = function(ds, i)
   {
    Y = ds[i,1];
    X = ds[i,2];
    U = ds[i, -c(1,2)];
    if(sum(alpha_Y) == 0)
    {
     model_Y0 = glm(Y~X, family = family.Y);
     model_X0 = glm(X~1, family = family.X);
    } else
    {
     model_Y0 = glm(Y~X + U[ , alpha_Y == 1], family = family.Y);
     model_X0 = glm(X~1 + U[ , alpha_Y == 1], family = family.X);
    }
    B0 = coef(model_Y0)[2]; # Initial estimate of slope
    Q0 = predict(model_Y0); # Y hat
    r0 = Q0 - B0*X; 
    predX0 = predict(model_X0); # X hat
    #Perform truncation ...
    Hg = X - predX0; # Clever covariate
    epsilon = coef(lm(Y~-1+Hg, offset = Q0)); # Fluctuation of estimate
    B1 = B0 + epsilon; # Final TMLE estimate of slope
    return(B1);
   }
   vars_t[1] = var(boot(data = cbind(Y, X, U), statistic = bootf, R = B)$t, na.rm = TRUE);
  }
 } else if(family.Y == "gaussian" & family.X == "binomial")
 {
  pX1 = plogis(cbind(1, U[ , alpha_Y == 1])%*%coef(model_X0));
  bounds = quantile(pX1, truncation); 
  pX1 = pmin(pX1, bounds[2]);
  pX1 = pmax(pX1, bounds[1]);
  H1 = X/pX1; H0 = (1 - X)/(1 - pX1);
  w = H1 + H0;
  Q1_0 = cbind(1, 1, U[,alpha_Y == 1])%*%coef(model_Y0);
  Q0_0 = cbind(1, 0, U[,alpha_Y == 1])%*%coef(model_Y0);
  QX_0 = cbind(1, X, U[,alpha_Y == 1])%*%coef(model_Y0);

  #TMLE
  epsilon = coef(lm(Y~1+X+offset(QX_0), weights = w));
  Q1_1 = Q1_0 + sum(epsilon);
  Q0_1 = Q0_0 + epsilon[1];
  QX_1 = Q1_1*X + Q0_1*(1-X);
  betas_t[1] = mean(Q1_1 - Q0_1);
  if(var.comp == "asymptotic")
  {
   vars_t[1] = var((H1 - H0)*(Y - QX_1) + Q1_1 - Q0_1 - betas_t[1,1])/n;
  } else #var == "bootstrap"
  {
   bootf = function(ds, i)
   {
    Y = ds[i,1];
    X = ds[i,2];
    U = ds[i, -c(1,2)];
    if(sum(alpha_Y) == 0)
    {
     model_Y0 = glm(Y~X, family = family.Y);
     model_X0 = glm(X~1, family = family.X);
    } else
    {
     model_Y0 = glm(Y~X + U[ , alpha_Y == 1], family = family.Y);
     model_X0 = glm(X~1 + U[ , alpha_Y == 1], family = family.X);
    }
    pX1 = plogis(cbind(1, U[ , alpha_Y == 1])%*%coef(model_X0));
    bounds = quantile(pX1, truncation); 
    pX1 = pmin(pX1, bounds[2]);
    pX1 = pmax(pX1, bounds[1]);
    H1 = X/pX1; H0 = (1 - X)/(1 - pX1);
    w = H1 + H0;
    Q1_0 = cbind(1, 1, U[,alpha_Y == 1])%*%coef(model_Y0);
    Q0_0 = cbind(1, 0, U[,alpha_Y == 1])%*%coef(model_Y0);
    QX_0 = cbind(1, X, U[,alpha_Y == 1])%*%coef(model_Y0);

    #TMLE
    epsilon = coef(lm(Y~1+X+offset(QX_0), weights = w));
    Q1_1 = Q1_0 + sum(epsilon);
    Q0_1 = Q0_0 + epsilon[1];
    QX_1 = Q1_1*X + Q0_1*(1-X);
    B1 = mean(Q1_1 - Q0_1);
    return(B1);
   }
   vars_t[1] = var(boot(data = cbind(Y, X, U), statistic = bootf, R = B)$t, na.rm = TRUE);
  }
 } else if(family.Y == "binomial" & family.X == "binomial")
 {
  pX1 = plogis(cbind(1, U[ , alpha_Y == 1])%*%coef(model_X0));
  bounds = quantile(pX1, truncation); 
  pX1 = pmin(pX1, bounds[2]);
  pX1 = pmax(pX1, bounds[1]);
  H1 = X/pX1; H0 = (1 - X)/(1 - pX1);
  w = H1 + H0;

  Q1_0 = plogis(cbind(1, 1, U[,alpha_Y == 1])%*%coef(model_Y0));
  Q0_0 = plogis(cbind(1, 0, U[,alpha_Y == 1])%*%coef(model_Y0));
  QX_0 = plogis(cbind(1, X, U[,alpha_Y == 1])%*%coef(model_Y0));

  #TMLE
  epsilon = coef(glm(Y~X+offset(qlogis(QX_0)), family = "binomial", weights = w));
  Q1_1 = plogis(qlogis(Q1_0) + sum(epsilon));
  Q0_1 = plogis(qlogis(Q0_0) + epsilon[1]);
  QX_1 = Q1_1*X + Q0_1*(1-X);

  # Risk difference
  betas_t[1,1] = mean(Q1_1 - Q0_1);

  # Relative risk
  D1 = (H1*(Y - Q1_1) + Q1_1 - mean(Q1_1));
  D0 = (H0*(Y - Q0_1) + Q0_1 - mean(Q0_1));
  betas_t[1,2] = mean(Q1_1)/mean(Q0_1);

  if(var.comp == "asymptotic")
  {
   vars_t[1,1] = var((H1 - H0)*(Y - QX_1) + Q1_1 - Q0_1 - betas_t[1,1])/n;
   vars_t[1,2] = var(D1/mean(Q0_1) - D0*mean(Q1_1)/mean(Q0_1)**2)/n;
  } else #var == "bootstrap"
  {
   bootf = function(ds, i)
   {
    Y = ds[i,1];
    X = ds[i,2];
    U = ds[i, -c(1,2)];
    if(sum(alpha_Y) == 0)
    {
     model_Y0 = glm(Y~X, family = family.Y);
     model_X0 = glm(X~1, family = family.X);
    } else
    {
     model_Y0 = glm(Y~X + U[ , alpha_Y == 1], family = family.Y);
     model_X0 = glm(X~1 + U[ , alpha_Y == 1], family = family.X);
    }
    pX1 = plogis(cbind(1, U[ , alpha_Y == 1])%*%coef(model_X0));
    bounds = quantile(pX1, truncation); 
    pX1 = pmin(pX1, bounds[2]);
    pX1 = pmax(pX1, bounds[1]);
    H1 = X/pX1; H0 = (1 - X)/(1 - pX1);
    w = H1 + H0;

    Q1_0 = plogis(cbind(1, 1, U[,alpha_Y == 1])%*%coef(model_Y0));
    Q0_0 = plogis(cbind(1, 0, U[,alpha_Y == 1])%*%coef(model_Y0));
    QX_0 = plogis(cbind(1, X, U[,alpha_Y == 1])%*%coef(model_Y0));

    #TMLE
    epsilon = coef(glm(Y~X+offset(qlogis(QX_0)), family = "binomial", weights = w));
    Q1_1 = plogis(qlogis(Q1_0) + sum(epsilon));
    Q0_1 = plogis(qlogis(Q0_0) + epsilon[1]);
    QX_1 = Q1_1*X + Q0_1*(1-X);

    # Risk difference
    B1 = mean(Q1_1 - Q0_1);

    # Relative risk
    D1 = (H1*(Y - Q1_1) + Q1_1 - mean(Q1_1));
    D0 = (H0*(Y - Q0_1) + Q0_1 - mean(Q0_1));
    B2 = mean(Q1_1)/mean(Q0_1);

    return(c(B1, B2));
   }
   boot.samples = boot(data = cbind(Y, X, U), statistic = bootf, R = B)$t;
   vars_t[1,1] = var(boot.samples[,1], na.rm = TRUE);
   vars_t[1,2] = var(boot.samples[,2], na.rm = TRUE);
  }
 } else if(family.Y == "binomial" & family.X == "gaussian")
 {
  #truncation ...
  #
  predX0 = predict(model_X0);
  sdX = sqrt(var(model_X0$residuals)*n/model_X0$df.residual);
  delta0 = X0 - X;
  delta1 = X1 - X;
  H0 = dnorm(X-delta0, mean = predX0, sd = sdX);
  H1 = dnorm(X-delta1, mean = predX0, sd = sdX);
  HX = dnorm(X, mean = predX0, sd = sdX);

  Q1 = plogis(cbind(1, X1, U[,alpha_Y == 1])%*%coef(model_Y0));
  Q0 = plogis(cbind(1, X0, U[,alpha_Y == 1])%*%coef(model_Y0));
  QX = plogis(cbind(1, X, U[,alpha_Y == 1])%*%coef(model_Y0));

  # Risk difference
  betas_t[1,1] = mean((H1-H0)/HX*(Y - QX) + Q1 - Q0);

  # Relative risk
  D1 = H1/HX*(Y - QX) + Q1;
  D0 = H0/HX*(Y - QX) + Q0;
  Q1_1 = mean(D1);
  Q0_1 = mean(D0);
  betas_t[1,2] = Q1_1/Q0_1;
  if(var.comp == "asymptotic")
  {
   vars_t[1,1] = var((H1-H0)/HX*(Y - QX) + Q1 - Q0 - betas_t[1,1])/n;
   vars_t[1,2] = var((D1 - Q1_1)/Q0_1 - (D0 - Q0_1)*Q1_1/Q0_1**2)/n;
  } else #var == "bootstrap"
  {
   bootf = function(ds, i)
   {
    Y = ds[i,1];
    X = ds[i,2];
    U = ds[i, -c(1,2)];
    if(sum(alpha_Y) == 0)
    {
     model_Y0 = glm(Y~X, family = family.Y);
     model_X0 = glm(X~1, family = family.X);
    } else
    {
     model_Y0 = glm(Y~X + U[ , alpha_Y == 1], family = family.Y);
     model_X0 = glm(X~1 + U[ , alpha_Y == 1], family = family.X);
    }
    predX0 = predict(model_X0);
    sdX = sqrt(var(model_X0$residuals)*n/model_X0$df.residual);
    delta0 = X0 - X;
    delta1 = X1 - X;
    H0 = dnorm(X-delta0, mean = predX0, sd = sdX);
    H1 = dnorm(X-delta1, mean = predX0, sd = sdX);
    HX = dnorm(X, mean = predX0, sd = sdX);

    Q1 = plogis(cbind(1, X1, U[,alpha_Y == 1])%*%coef(model_Y0));
    Q0 = plogis(cbind(1, X0, U[,alpha_Y == 1])%*%coef(model_Y0));
    QX = plogis(cbind(1, X, U[,alpha_Y == 1])%*%coef(model_Y0));

    # Risk difference
    B1 = mean((H1-H0)/HX*(Y - QX) + Q1 - Q0);

    # Relative risk
    D1 = H1/HX*(Y - QX) + Q1;
    D0 = H0/HX*(Y - QX) + Q0;
    Q1_1 = mean(D1);
    Q0_1 = mean(D0);
    B2 = Q1_1/Q0_1;

    return(c(B1, B2));
   }
   boot.samples = boot(data = cbind(Y, X, U), statistic = bootf, R = B)$t;
   vars_t[1,1] = var(boot.samples[,1], na.rm = TRUE);
   vars_t[1,2] = var(boot.samples[,2], na.rm = TRUE);
  }
 }

 models_Y[1,] = alpha_Y;
 tested_models_Y[1] = paste0(as.numeric(alpha_Y), collapse = "");
 pY_tested[1] = bic_y;
 pY_tested3[1] = 1;

 ratio_init = 1; #posterior ratio vs initial model

 if(sum(alpha_Y != 0))
 {
  a = coef(model_Y0)[-(1:2)];
  b = diag(vcov(model_Y0))[-(1:2)]**0.5;
 } else 
 {
  a = NA;
  b = NA;
 }
 pY_tested1[1, 1:sum(alpha_Y)] = a;
 pY_tested2[1, 1:sum(alpha_Y)] = b;

 for(i in 2:niter)
 {
  alpha_Y_0 = alpha_Y; #Current alpha_Y
  if(sum(alpha_Y) < maxsize){
   alpha_Y_1 = alpha_Y; temp = sample(n_cov, 1); alpha_Y_1[temp] = (alpha_Y_1[temp] + 1)%%2; #Candidate alpha_Y
  } else {
   alpha_Y_1 = alpha_Y; temp = sample((1:n_cov)[alpha_Y == 1], 1); alpha_Y_1[temp] = (alpha_Y_1[temp] + 1)%%2; #Candidate alpha_Y
  }
  bic_y_0 = bic_y;
  a_0 = a;
  b_0 = b;
  char_alpha_Y_1 = paste0(alpha_Y_1, collapse = "");
  tested = which(tested_models_Y == char_alpha_Y_1);
  if(length(tested)) #The candidate model was already tested
  {
   a_1 = as.numeric(na.omit(pY_tested1[tested,]));
   b_1 = as.numeric(na.omit(pY_tested2[tested,]));
   bic_y_1 = pY_tested[tested];
   ratio = pY_tested3[tested]/ratio_init;
  } else #The candidate model was never tested
  {
   if(sum(alpha_Y_1) == 0)
   {
    model_Y1 = glm(Y~X, family = family.Y);
    model_X0 = glm(X~1, family = family.X);
    a_1 = NA;
    b_1 = NA;
   } else
   {
    model_Y1 = glm(Y~X + U[ , alpha_Y_1 == 1], family = family.Y);
    model_X0 = glm(X~1 + U[ , alpha_Y_1 == 1], family = family.X);
    a_1 = coef(model_Y1)[-(1:2)];
    b_1 = diag(vcov(model_Y1))[-(1:2)]**0.5;
   }
   pY_tested1[i, 1:sum(alpha_Y_1) ] = a_1;
   pY_tested2[i, 1:sum(alpha_Y_1) ] = b_1;
   bic_y_1 = BIC(model_Y1);
   tested_models_Y[i] = char_alpha_Y_1; #Recording info about candidate model
   pY_tested[i] = bic_y_1; #Recording info about candidate model
   if(family.Y == "gaussian" & family.X == "gaussian")
   {
    #TMLE
    B0 = coef(model_Y1)[2]; # Initial estimate of slope
    Q0 = predict(model_Y1); # Y hat
    r0 = Q0 - B0*X; 
    predX0 = predict(model_X0); # X hat
    #Truncation ...
    Hg = X - predX0; # Clever covariate
    epsilon = coef(lm(Y~-1+Hg, offset = Q0)); # Fluctuation of estimate
    B1 = betas_t[i] = B0 + epsilon; # Final TMLE estimate of slope

    if(var.comp == "asymptotic")
    {
     # Compute empirical influence curve for variance estimation
     r1 = r0 - epsilon*predX0; 
     Q1 = B1*X + r1; # Updated Y hat
     k = -mean(Hg*X); 
     D = Hg*(Y - Q1); # Score
     IC = 1/k*D; # Influence curve
     vars_t[i] = mean(IC**2)/n;
    } else #var == "bootstrap"
    {
     bootf = function(ds, i)
     {
      Y = ds[i,1];
      X = ds[i,2];
      U = ds[i, -c(1,2)];
      if(sum(alpha_Y_1) == 0)
      {
       model_Y1 = glm(Y~X, family = family.Y);
       model_X0 = glm(X~1, family = family.X);
      } else
      {
       model_Y1 = glm(Y~X + U[ , alpha_Y_1 == 1], family = family.Y);
       model_X0 = glm(X~1 + U[ , alpha_Y_1 == 1], family = family.X);
      }
      #TMLE
      B0 = coef(model_Y1)[2]; # Initial estimate of slope
      Q0 = predict(model_Y1); # Y hat
      r0 = Q0 - B0*X; 
      predX0 = predict(model_X0); # X hat
      #Truncation ...
      Hg = X - predX0; # Clever covariate
      epsilon = coef(lm(Y~-1+Hg, offset = Q0)); # Fluctuation of estimate
      B1 = B0 + epsilon; # Final TMLE estimate of slope
      return(B1);
     }
     vars_t[i] = var(boot(data = cbind(Y, X, U), statistic = bootf, R = B)$t, na.rm = TRUE);
    }
   } else if(family.Y == "gaussian" & family.X == "binomial")
   {
    pX1 = plogis(cbind(1, U[ , alpha_Y_1 == 1])%*%coef(model_X0));
    bounds = quantile(pX1, truncation); 
    pX1 = pmin(pX1, bounds[2]);
    pX1 = pmax(pX1, bounds[1]);

    H1 = X/pX1; H0 = (1 - X)/(1 - pX1);
    w = H1 + H0;
    Q1_0 = cbind(1, 1, U[,alpha_Y_1 == 1])%*%coef(model_Y1);
    Q0_0 = cbind(1, 0, U[,alpha_Y_1 == 1])%*%coef(model_Y1);
    QX_0 = cbind(1, X, U[,alpha_Y_1 == 1])%*%coef(model_Y1);

    #TMLE
    epsilon = coef(lm(Y~X+offset(QX_0), weights = w));
    Q1_1 = Q1_0 + sum(epsilon);
    Q0_1 = Q0_0 + epsilon[1];
    QX_1 = Q1_1*X + Q0_1*(1-X);
    betas_t[i] = mean(Q1_1 - Q0_1);
    if(var.comp == "asymptotic")
    {
     vars_t[i] = var((H1 - H0)*(Y - QX_1) + Q1_1 - Q0_1 - betas_t[i,1])/n;
    } else #var == "bootstrap"
    {
     bootf = function(ds, i)
     {
      Y = ds[i,1];
      X = ds[i,2];
      U = ds[i, -c(1,2)];
      if(sum(alpha_Y_1) == 0)
      {
       model_Y1 = glm(Y~X, family = family.Y);
       model_X0 = glm(X~1, family = family.X);
      } else
      {
       model_Y1 = glm(Y~X + U[ , alpha_Y_1 == 1], family = family.Y);
       model_X0 = glm(X~1 + U[ , alpha_Y_1 == 1], family = family.X);
      }
      pX1 = plogis(cbind(1, U[ , alpha_Y_1 == 1])%*%coef(model_X0));
      bounds = quantile(pX1, truncation); 
      pX1 = pmin(pX1, bounds[2]);
      pX1 = pmax(pX1, bounds[1]);

      H1 = X/pX1; H0 = (1 - X)/(1 - pX1);
      w = H1 + H0;
      Q1_0 = cbind(1, 1, U[,alpha_Y_1 == 1])%*%coef(model_Y1);
      Q0_0 = cbind(1, 0, U[,alpha_Y_1 == 1])%*%coef(model_Y1);
      QX_0 = cbind(1, X, U[,alpha_Y_1 == 1])%*%coef(model_Y1);

      #TMLE
      epsilon = coef(lm(Y~X+offset(QX_0), weights = w));
      Q1_1 = Q1_0 + sum(epsilon);
      Q0_1 = Q0_0 + epsilon[1];
      QX_1 = Q1_1*X + Q0_1*(1-X);
      B1 = mean(Q1_1 - Q0_1);
      return(B1);
     }
     vars_t[i] = var(boot(data = cbind(Y, X, U), statistic = bootf, R = B)$t, na.rm = TRUE);
    }
   }else if(family.Y == "binomial" & family.X == "binomial")
   {
    pX1 = plogis(cbind(1, U[ , alpha_Y_1 == 1])%*%coef(model_X0));
    bounds = quantile(pX1, truncation); 
    pX1 = pmin(pX1, bounds[2]);
    pX1 = pmax(pX1, bounds[1]);
    H1 = X/pX1; H0 = (1 - X)/(1 - pX1);
    w = H1 + H0;

    Q1_0 = plogis(cbind(1, 1, U[,alpha_Y_1 == 1])%*%coef(model_Y1));
    Q0_0 = plogis(cbind(1, 0, U[,alpha_Y_1 == 1])%*%coef(model_Y1));
    QX_0 = plogis(cbind(1, X, U[,alpha_Y_1 == 1])%*%coef(model_Y1));

    #TMLE
    epsilon = coef(glm(Y~X+offset(qlogis(QX_0)), family = "binomial", weights = w));
    Q1_1 = plogis(qlogis(Q1_0) + sum(epsilon));
    Q0_1 = plogis(qlogis(Q0_0) + epsilon[1]);
    QX_1 = Q1_1*X + Q0_1*(1-X);

    # Risk difference
    betas_t[i,1] = mean(Q1_1 - Q0_1);

    # Relative risk
    D1 = (H1*(Y - Q1_1) + Q1_1 - mean(Q1_1));
    D0 = (H0*(Y - Q0_1) + Q0_1 - mean(Q0_1));
    betas_t[i,2] = mean(Q1_1)/mean(Q0_1);

    if(var.comp == "asymptotic")
    {
     vars_t[i,1] = var((H1 - H0)*(Y - QX_1) + Q1_1 - Q0_1 - betas_t[i,1])/n;
     vars_t[i,2] = var(D1/mean(Q0_1) - D0*mean(Q1_1)/mean(Q0_1)**2)/n;
    } else #var == "bootstrap"
    {
     bootf = function(ds, i)
     {
      Y = ds[i,1];
      X = ds[i,2];
      U = ds[i, -c(1,2)];
      if(sum(alpha_Y_1) == 0)
      {
       model_Y1 = glm(Y~X, family = family.Y);
       model_X0 = glm(X~1, family = family.X);
      } else
      {
       model_Y1 = glm(Y~X + U[ , alpha_Y_1 == 1], family = family.Y);
       model_X0 = glm(X~1 + U[ , alpha_Y_1 == 1], family = family.X);
      }
      pX1 = plogis(cbind(1, U[ , alpha_Y_1 == 1])%*%coef(model_X0));
      bounds = quantile(pX1, truncation); 
      pX1 = pmin(pX1, bounds[2]);
      pX1 = pmax(pX1, bounds[1]);
      H1 = X/pX1; H0 = (1 - X)/(1 - pX1);
      w = H1 + H0;

      Q1_0 = plogis(cbind(1, 1, U[,alpha_Y_1 == 1])%*%coef(model_Y1));
      Q0_0 = plogis(cbind(1, 0, U[,alpha_Y_1 == 1])%*%coef(model_Y1));
      QX_0 = plogis(cbind(1, X, U[,alpha_Y_1 == 1])%*%coef(model_Y1));

      #TMLE
      epsilon = coef(glm(Y~X+offset(qlogis(QX_0)), family = "binomial", weights = w));
      Q1_1 = plogis(qlogis(Q1_0) + sum(epsilon));
      Q0_1 = plogis(qlogis(Q0_0) + epsilon[1]);
      QX_1 = Q1_1*X + Q0_1*(1-X);

      # Risk difference
      B1 = mean(Q1_1 - Q0_1);

      # Relative risk
      D1 = (H1*(Y - Q1_1) + Q1_1 - mean(Q1_1));
      D0 = (H0*(Y - Q0_1) + Q0_1 - mean(Q0_1));
      B2 = mean(Q1_1)/mean(Q0_1);
      return(c(B1, B2));
     }
     boot.samples = boot(data = cbind(Y, X, U), statistic = bootf, R = B)$t;
     vars_t[i,1] = var(boot.samples[,1], na.rm = TRUE);
     vars_t[i,2] = var(boot.samples[,2], na.rm = TRUE);
    }
   } else if(family.Y == "binomial" & family.X == "gaussian")
   {
    predX0 = predict(model_X0);
    #Truncation...
    sdX = sqrt(var(model_X0$residuals)*n/model_X0$df.residual);
    H0 = dnorm(X-delta0, mean = predX0, sd = sdX);
    H1 = dnorm(X-delta1, mean = predX0, sd = sdX);
    HX = dnorm(X, mean = predX0, sd = sdX);

    Q1 = plogis(cbind(1, X1, U[,alpha_Y_1 == 1])%*%coef(model_Y1));
    Q0 = plogis(cbind(1, X0, U[,alpha_Y_1 == 1])%*%coef(model_Y1));
    QX = plogis(cbind(1, X, U[,alpha_Y_1 == 1])%*%coef(model_Y1));

    # A-IPTW
    # Risk difference
    betas_t[i,1] = mean((H1-H0)/HX*(Y - QX) + Q1 - Q0);

    # Relative risk
    D1 = H1/HX*(Y - QX) + Q1;
    D0 = H0/HX*(Y - QX) + Q0;
    Q1_1 = mean(D1);
    Q0_1 = mean(D0);
    betas_t[i,2] = Q1_1/Q0_1;

    if(var.comp == "asymptotic")
    {
     vars_t[i,1] = var((H1-H0)/HX*(Y - QX) + Q1 - Q0 - betas_t[i,1])/n;
     vars_t[i,2] = var((D1 - Q1_1)/Q0_1 - (D0 - Q0_1)*Q1_1/Q0_1**2)/n;
    } else #var == "bootstrap"
    {
     bootf = function(ds, i)
     {
      Y = ds[i,1];
      X = ds[i,2];
      U = ds[i, -c(1,2)];
      if(sum(alpha_Y_1) == 0)
      {
       model_Y1 = glm(Y~X, family = family.Y);
       model_X0 = glm(X~1, family = family.X);
      } else
      {
       model_Y1 = glm(Y~X + U[ , alpha_Y_1 == 1], family = family.Y);
       model_X0 = glm(X~1 + U[ , alpha_Y_1 == 1], family = family.X);
      }
      predX0 = predict(model_X0);
      #Truncation...
      sdX = sqrt(var(model_X0$residuals)*n/model_X0$df.residual);
      H0 = dnorm(X-delta0, mean = predX0, sd = sdX);
      H1 = dnorm(X-delta1, mean = predX0, sd = sdX);
      HX = dnorm(X, mean = predX0, sd = sdX);

      Q1 = plogis(cbind(1, X1, U[,alpha_Y_1 == 1])%*%coef(model_Y1));
      Q0 = plogis(cbind(1, X0, U[,alpha_Y_1 == 1])%*%coef(model_Y1));
      QX = plogis(cbind(1, X, U[,alpha_Y_1 == 1])%*%coef(model_Y1));

      # A-IPTW
      # Risk difference
      B1 = mean((H1-H0)/HX*(Y - QX) + Q1 - Q0);

      # Relative risk
      D1 = H1/HX*(Y - QX) + Q1;
      D0 = H0/HX*(Y - QX) + Q0;
      Q1_1 = mean(D1);
      Q0_1 = mean(D0);
      B2 = Q1_1/Q0_1;
      return(c(B1, B2));
     }
     boot.samples = boot(data = cbind(Y, X, U), statistic = bootf, R = B)$t;
     vars_t[i,1] = var(boot.samples[,1], na.rm = TRUE);
     vars_t[i,2] = var(boot.samples[,2], na.rm = TRUE);
    }
   }	
   models_Y[i,] = alpha_Y_1;
   change_Y = alpha_Y_1 - alpha_Y_0;
   px = alpha_X[change_Y != 0];
   add_Y = sum(change_Y); #Is it proposed to add or remove a variable
   if(add_Y == 1)
   {
    delta = omega*(a_1[change_Y[alpha_Y_1 != 0] != 0]*su[change_Y != 0]/sy +
                   norm.sample*b_1[change_Y[alpha_Y_1 != 0] != 0]*su[change_Y != 0]/sy)**2;
    ratio_P_alpha = mean((px*(delta/(1 + delta)) + (1 - px)*0.5))/mean((px*(1/(1 + delta)) + (1 - px)*0.5));
   } else
   {
    delta = omega*(a_0[change_Y[alpha_Y_0 != 0] != 0]*su[change_Y != 0]/sy +
                   norm.sample*b_0[change_Y[alpha_Y_0 != 0] != 0]*su[change_Y != 0]/sy)**2;
    ratio_P_alpha = mean((px*(1/(1 + delta)) + (1 - px)*0.5))/mean((px*(delta/(1 + delta)) + (1 - px)*0.5));
   }
   ratio = exp(-bic_y_1/2 + bic_y_0/2)*ratio_P_alpha*prod((priorY/(1-priorY))**change_Y);
   pY_tested3[i] = ratio*ratio_init;
  }
		
  if(sample(c(1,0), size = 1, prob = c(min(1, ratio), 1 - min(1, ratio))))
  {
   alpha_Y = alpha_Y_1;
   bic_y = bic_y_1;
   a = a_1;
   b = b_1;
   ratio_init = ratio*ratio_init;
  }
  #if(!length(a)) break;
 }
 models_Y = na.omit(models_Y);
 pY_tested3 = as.numeric(na.omit(pY_tested3));
 pY_tested3[pY_tested3/max(pY_tested3) < 1/OR] = 0;
 pY_tested3 = pY_tested3/sum(pY_tested3);

 #TMLE
 betas_t = na.omit(betas_t);
 vars_t = na.omit(vars_t);

 models.Y = cbind(models_Y, betas_t, sqrt(vars_t), pY_tested3);

 if(family.Y == "gaussian")
 {
  colnames(models.Y) = c(1:n_cov, "Diff", "SE", "PostProb");
  ord = order(models.Y[,n_cov+3], decreasing = TRUE);
  beta_t = colSums(betas_t*pY_tested3);
  stderr_t = sqrt(colSums((vars_t + betas_t**2)*pY_tested3) - beta_t**2); 
 } else #family.Y == "binomial"
 {
  colnames(models.Y) = c(1:n_cov, "Diff", "RR", "SE.Diff", "SE.RR", "PostProb");
  ord = order(models.Y[,n_cov+5], decreasing = TRUE);
  beta_t = colSums(betas_t*pY_tested3);
  stderr_t = sqrt(colSums((vars_t + betas_t**2)*pY_tested3) - beta_t**2); 
  names(beta_t) = c("Diff", "RR");
  names(stderr_t) = c("Diff", "RR");
 }
 return(list(beta = beta_t, stderr = stderr_t, models.X = models.X, models.Y = models.Y[ord,]));
}







