#Packages
library(panelPomp)
library(adaptMCMC)

print("pt 1")
##### Imports current best fit parameters and jump covariances #####
# Gets the best fit parameters and the covariance matrix
get_best_fit = function(chain.list){
  L = length(chain.list)
  chain.scores = numeric()
  for(i in 1:L){
    chain.scores[i] = max(chain.list[[i]]$log.p)
  }
  list(chain.list[[which.max(chain.scores)]]$samples[which.max(chain.list[[which.max(chain.scores)]]$log.p),],
       chain.list[[which.max(chain.scores)]]$cov.jump)
  
}

print("pt 2")
# Read in the previous best-fit parameters and estimated var-covar matrix
samps = readRDS("Joint_GW_full_A.RDA")

print("pt 3")
pars = get_best_fit(samps)[[1]]
variances = get_best_fit(samps)[[2]]

###############################################################


##### Rebound experiment data #####
# Imports the data and experimental conditions from the Rebound experiment
Rebound_Time_Series <- read.csv("Data_Rebound_4-2-25.csv")
Rebound_exp_df <- read.csv("Rebound_experimental_conditions.csv")

# Defines the days that harvesting actually occurred
Harvest_days = c(8, 36, 64)

# This converts target of proportion mortality imposed for a single day into an instantaneous death rate for use in model
Rebound_exp_df[,"Harvest_rate"] = -log(1 - Rebound_exp_df[,"Harvest"]/100)*(Rebound_exp_df$Harvest_days %in% Harvest_days)
###############################################################
print("pt 4")
##### Rebound experiment pomp definition #####
# Initial conditions function for pomp
GW_R_rinit <- Csnippet("
  N = rpois(N0_R);
  J = rpois(J0_R);
  A = rpois(A0_R);
")

# r process model
GW_R_step_process <- Csnippet("
  double VOL = 15;
  
  double Pred_A = Harvest_rate;
  double Pred_J = Harvest_rate;
  double Pred_N = Harvest_rate * Sieve;
    if(Sieve == 0){ Pred_N = Pred_N * naup_catch;}
  
  double birth_rate = b_M*exp(-comp_b / VOL * (c_N * N + c_J * J + A));

  double d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A));
  double d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A));
  double d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A)) + cann*A/VOL;

  double m_N_c = m_N*exp(-comp_m / VOL * (c_N * N + c_J * J + A));
  double m_J_c = m_J*exp(-comp_m / VOL * (c_N * N + c_J * J + A));


  /* Specify adult rates */
  const double A_rates[2] = {d_A_c, Pred_A};
  double A_trans[2];
  reulermultinom(2, A, &A_rates, dt, &A_trans);
  
  /* Specify juvenile rates */
  const double J_rates[3] = {m_J_c, d_J_c, Pred_J};
  double J_trans[3];
  reulermultinom(3, J, &J_rates, dt, &J_trans);
  
  /* Specify nauplii rates */
  const double N_rates[3] = {m_N_c, d_N_c, Pred_N};
  double N_trans[3];
  reulermultinom(3, N, &N_rates, dt, &N_trans);
  
  /* births */
  int births = rpois(birth_rate * A / 2 * dt);
  
  /* New totals */
  int A_out = A - A_trans[0] - A_trans[1] + J_trans[0];
  int J_out = J - J_trans[0] - J_trans[1] - J_trans[2] + N_trans[0];
  int N_out = N - N_trans[0] - N_trans[1] - N_trans[2] + births;
  N = N_out >= 0 ? N_out : 0;
  J = J_out >= 0 ? J_out : 0;
  A = A_out >= 0 ? A_out : 0;
")

# r measure model
GW_R_rmeasure <- function(N, J, A, counted_volume, VOL = 15, SAMPLEVOL = 0.5, aliquot_volume = 10, k_R, ...) {
  c(
    N1 = rnbinom(n = 1, mu = N * counted_volume / aliquot_volume * SAMPLEVOL / VOL, size = k_R),
    JOA1 = rnbinom(n = 1, mu = (J + 0.5 * (1 + 1 / 3) * A) * counted_volume / aliquot_volume * SAMPLEVOL / VOL, size = k_R),
    AF1 = rnbinom(n = 1, mu = 0.5 * 2 / 3 * A * counted_volume / aliquot_volume * SAMPLEVOL / VOL, size = k_R)
  )
}
# d measure model
GW_R_dmeasure <- Csnippet("
  double VOL = 15;
  double SAMPLEVOL = 0.5;
  double aliquot_volume = 10;
  lik = dnbinom_mu(N1, k_R, fmax(N, 1) * counted_volume / aliquot_volume * SAMPLEVOL / VOL, give_log) +
        dnbinom_mu(JOA1, k_R, fmax(J + 0.5 * (1 + 1.0 / 3.0) * A, 1) * counted_volume / aliquot_volume * SAMPLEVOL / VOL, give_log) +
        dnbinom_mu(AF1, k_R, fmax(0.5 * 2.0 / 3.0 * A, 1) * counted_volume / aliquot_volume * SAMPLEVOL / VOL, give_log);
")
###############################################################
print("pt 5")
##### Constructs the panel pomp object for the Rebound experiment #####
# Template with updated parameters
template_R <- pomp(
  data = data.frame(
    Day = 0:17*7, #This needs to match the times in the actual dataset
    Tank = NA,
    N1 = NA, JOA1 = NA, AF1 = NA
  ),
  times = "Day", t0 = 0,
  rprocess = discrete_time(GW_R_step_process, delta.t = 1),
  rmeasure = GW_R_rmeasure,
  dmeasure = GW_R_dmeasure,
  rinit = GW_R_rinit,
  obsnames = c("N1", "JOA1", "AF1"),
  statenames = c("N", "J", "A"),
  covarnames = c("Sieve", "Harvest_rate", "counted_volume"),
  paramnames = c("N0_R", "J0_R", "A0_R", "b_M", "comp_b", "comp_d", "comp_m", "cann", "c_N", "c_J", "d_A", "m_J", "d_J", "m_N", "d_N", "naup_catch", "k_R"),
  params = pars[c("N0_R", "J0_R", "A0_R", "b_M", "comp_b", "comp_d", "comp_m", "cann", "c_N", "c_J", "d_A", "m_J", "d_J", "m_N", "d_N", "naup_catch", "k_R")]
)

print("pt 6")
# Checking how many tanks there are
U <- length(unique(Rebound_exp_df$Tank))
# Empty list for the pomp templates
poList <- setNames(vector(mode = "list", length = U),
                   nm = paste0("unit", 1:U))
# Populates each element of the list with the model template, covariates, and the data
for (u in seq_len(U)) {
  covariates <- subset(Rebound_exp_df, Tank == u, select = c("Harvest_days", "Sieve", "Harvest_rate", "counted_volume"))
  cov_table <- covariate_table(covariates, times = "Harvest_days")
  poList[[u]] <- pomp(template_R, covar = cov_table)
  data_u <- pomp(
    data = subset(Rebound_Time_Series, Tank == u, select = c("Day", "Tank", "N1", "JOA1", "AF1")),
    times="Day",
    t0=timezero(template_R)
  )
  poList[[u]]@data <- data_u@data
}

# Creates the panel Pomp object
GW_R_panel <- panelPomp(object = poList, shared = coef(template_R))
#plot(simulate(GW_R_panel, nsim = 1))
###############################################################

##### Brings in data for fish experiment #####
# Brings in data, note this has been processed in other scripts to avoid using tidyverse on cluster
Fish_Time_Series = read.csv("Fish_experiment_data.csv")
###############################################################

##### Pomp construction Fish experiment #####
# Initial conditions function for pomp
GW_F_rinit <- Csnippet("
  N = rpois(N0_F);
  J = rpois(J0_F);
  A = rpois(A0_F);
")

# r process model
GW_F_step_process <- Csnippet("
  double VOL = 40;
  
  /* These three equations specify predation rates on A,J,N respectiviely. f is predation rate on adults. f_j and f_n, are relative predation rates on juveniles and nauplii.*/
  /* Similarly, h is handeling time. h_J and h_N are relative handeling times. Lastly, i_P representats strength of interference.*/
  double Pred_A = f*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J + f_N*h_N*N)/VOL  + i_P*fmax(Preds/VOL - 1/VOL, 0)/VOL);
  double Pred_J = f*f_J*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J + f_N*h_N*N)/VOL  + i_P*fmax(Preds/VOL - 1/VOL, 0)/VOL);
  double Pred_N = f*f_N*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J + f_N*h_N*N)/VOL  + i_P*fmax(Preds/VOL - 1/VOL, 0)/VOL);
  
  double birth_rate = b_M*exp(-comp_b / VOL * (c_N * N + c_J * J + A));

  double d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A));
  double d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A));
  double d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A)) + cann*A/VOL;
  
  double m_N_c = m_N*exp(-comp_m / VOL * (c_N * N + c_J * J + A));
  double m_J_c = m_J*exp(-comp_m / VOL * (c_N * N + c_J * J + A));


  /* Specify adult rates */
  const double A_rates[2] = {d_A_c, Pred_A};
  double A_trans[2];
  reulermultinom(2, A, &A_rates, dt, &A_trans);
  
  /* Specify juvenile rates */
  const double J_rates[3] = {m_J_c, d_J_c, Pred_J};
  double J_trans[3];
  reulermultinom(3, J, &J_rates, dt, &J_trans);
  
  /* Specify nauplii rates */
  const double N_rates[3] = {m_N_c, d_N_c, Pred_N};
  double N_trans[3];
  reulermultinom(3, N, &N_rates, dt, &N_trans);
  
  /* births */
  int births = rpois(birth_rate * A / 2 * dt);
  
  /* New totals */
  int A_out = A - A_trans[0] - A_trans[1] + J_trans[0];
  int J_out = J - J_trans[0] - J_trans[1] - J_trans[2] + N_trans[0];
  int N_out = N - N_trans[0] - N_trans[1] - N_trans[2] + births;
  N = N_out >= 0 ? N_out : 0;
  J = J_out >= 0 ? J_out : 0;
  A = A_out >= 0 ? A_out : 0;
")

# r measure model - for some reason, Pomp will not compile when we make this a C-snippet (all others work just fine)
GW_F_rmeasure = function(N,J,A,VOL=40,SAMPLEVOL=0.5, counted_volume, aliquot_volume=10, k_F, ...){
  c(
    N1 = rnbinom(n=1, mu=N*counted_volume/aliquot_volume*SAMPLEVOL/VOL, size=k_F), 
    JOA1 = rnbinom(n=1, mu=(J + 0.5*(1 + 1/3)*A)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, size=k_F),
    AF1 = rnbinom(n=1, mu=0.5*2/3*A*counted_volume/aliquot_volume*SAMPLEVOL/VOL, size=k_F),
    N2 = rnbinom(n=1, mu=N*counted_volume/aliquot_volume*SAMPLEVOL/VOL, size=k_F), 
    JOA2 = rnbinom(n=1, mu=(J + 0.5*(1 + 1/3)*A)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, size=k_F),
    AF2 = rnbinom(n=1, mu=0.5*2/3*A*counted_volume/aliquot_volume*SAMPLEVOL/VOL, size=k_F))
}

# d measure model
GW_F_dmeasure = Csnippet("
  double VOL = 40;
  double SAMPLEVOL = 0.5;
  double aliquot_volume = 10;
  lik = dnbinom_mu(N1, k_F,fmax(N,1)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, give_log) +
        dnbinom_mu(JOA1, k_F,  fmax(J+ 0.5*(1 + 1.0/3.0)*A,1)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, give_log) +
        dnbinom_mu(AF1, k_F, fmax(0.5*2.0/3.0*A,1)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, give_log) +
        dnbinom_mu(N2, k_F, fmax(N,1)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, give_log) +
        dnbinom_mu(JOA2, k_F, fmax(J+ 0.5*(1 + 1.0/3.0)*A,1)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, give_log) +
        dnbinom_mu(AF2, k_F, fmax(0.5*2.0/3.0*A,1)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, give_log);"
)
###############################################################

##### Panel pomp construction for Fish experiment #####

# Template pomp to eventually construct the panelpomp
template_F = pomp(
  data =data.frame( #the template structure needs a dummy data set with the right column names
    Day = c(28, 42, 49, 56, 63, 70, 77, 84, 91, 98), #days are based on experimental sampling days
    Tank = NA, 
    total_observedN = NA, total_JA = NA), times="Day", t0=28,
  rprocess = discrete_time(GW_F_step_process, delta.t=1),
  rmeasure = GW_F_rmeasure,
  dmeasure = GW_F_dmeasure,
  rinit = GW_F_rinit,
  obsnames = c("N1", "JOA1", "AF1", "N2", "JOA2", "AF2"),
  statenames = c("N","J","A"),
  covarnames = c("Preds", "counted_volume"),
  paramnames = c("N0_F", "J0_F", "A0_F", "f", "f_J", "f_N", "h", "h_J", "h_N",  "b_M", "comp_b", "comp_d", "comp_m", "cann", "c_N", "c_J", "d_A", "m_J", "d_J", "m_N", "d_N", 
                 "k_F", "i_P"),
  params = pars[c("N0_F", "J0_F", "A0_F", "f", "f_J", "f_N", "h", "h_J", "h_N",  "b_M", "comp_b", "comp_d", "comp_m", "cann", "c_N", "c_J", "d_A", "m_J", "d_J", "m_N", "d_N", 
                    "k_F", "i_P")]
)

# Checking how many tanks there are
U <- length(unique(Fish_Time_Series$Tank))
# Empty list for the pomp templates
poList <- setNames(vector(mode = "list", length = U),
                   nm = paste0("unit", 1:U))

# This loop populates each element of the list with the model template, covariates, and the data (Fish_Time_Series)
for (u in seq_len(U)) {
  # For each tank, we need to grab the covariates, add the starting condition covariates, and get the covariate timing right
  #                for fish introduction
  covariates <- subset(Fish_Time_Series, Tank == u, select=c("Day", "Live_Fish", "counted_volume"))
  Pred_treat <- subset(Fish_Time_Series, Tank == u, select=c("Preds"))
  covariates <- rbind(covariates, c(50, max(Pred_treat, 0), 0)) # Get covariates for day of fish addition
  colnames(covariates)[1:2] <- c("FishDay", "Preds")  # time for covariates cannot have the same name as time for observed data
  # orders covariate table by time, which pomp() wants
  covariates = covariates[order(covariates$FishDay),]
  cov_table = covariate_table(covariates, times="FishDay")
  poList[[u]] <- pomp(template_F, covar = cov_table)
  data_u <- pomp(
    data = subset(Fish_Time_Series, Tank == u, select = c("Day", "Tank", "N1", "JOA1", "AF1", "N2", "JOA2", "AF2")),
    times="Day",
    t0=timezero(template_F)
  )
  poList[[u]]@data <- data_u@data
}

# Creates the panel Pomp object
GW_F_panel <- panelPomp(object = poList, shared = coef(template_F))
#plot(simulate(GW_F_panel, nsim = 1))
###############################################################

##### Define prior and model likelihoods to generate full likelihood #####

prior.likelihood = function(x){
  prior.lik = with(as.list(x),
                   # Initial conditions for both experiments
                   sum(dunif(c(N0_R, J0_R, A0_R, N0_F, J0_F, A0_F, cann), min=0, max=1000000, log=T)) + # maximum search rate of type 2 fxn response and N_0s are unknown
                     
                     # Informative priors from Ta'Nyia's small-scale experiments
                     dnorm(b_M, mean=5.67, sd=0.19, log=T) + # Estimated from the lab (egg data)                     
                     dnorm(m_N, mean = 0.457, sd = 0.11, log=T) + # maturation rate estimated in lab
                     dnorm(m_J, mean=0.067, sd = 0.035, log=T) + # maturation rate estimated in lab                     
                     dbeta(naup_catch, 416.25, 508.75, log=T) +
                     
                     # Relative competition parameters for nauplii, juveniles
                     sum(dbeta(c(c_J, c_N), 1, 1, log=T)) + # These parameters are unknown between 0-1
                     
                     # Adult competitive effects on births and deaths
                     sum(dunif(c(comp_b, comp_d, comp_m), min=0, max=0.01, log=T)) +
                     
                     # Background death rates
                     sum(dunif(c(d_N, d_J, d_A), min=0, max=1, log=T)) +
                     
                     # Predation parameters
                     sum(dbeta(c(f_J, f_N, h_J, h_N), 1, 1, log=T)) + # These parameters are unknown between 0-1
                     sum(dunif(c(f, i_P), min=0, max=1000000, log=T)) + # maximum search rate of type 2 fxn response 
                     dunif(h, min=0, max=0.02, log=T) + # handling time is less than 0.02 days
                     
                     # Observation model overdispersion parameter
                     sum(dunif(c(k_F, k_R), min=0, max=100, log=T))
                   
  )
  return(prior.lik)
}

prior.likelihood(pars)

full_LL = function(x){
  prior_LL = prior.likelihood(x)
  if(!is.finite(prior_LL)){return(prior_LL)}else{
    (prior_LL  + pfilter(GW_R_panel, shared = x, Np=50)@ploglik + pfilter(GW_F_panel, shared = x, Np=50)@ploglik)
  }
}

full_LL(pars)


var(replicate(full_LL(pars), n=100))



generate_pars_list = function(pars){
  LL = full_LL(pars)
  buffer = 1.1
  LL_check = (buffer+1)*LL
  attempts = 1
  while(LL_check < buffer*LL){
    #scale = sample(c(0.1, 0.2, 0.33, 0.5, rep(1, times=10), 2, 3, 5, 10), size=length(pars), replace=T)
    scale = rlnorm(n=length(pars), meanlog=0, sdlog=0.9)
    rand_pars = pars*scale
    rand_pars[c("c_N", "c_J", "naup_catch", "f_J", "f_N", "h_J", "h_N")] = pmin(rand_pars[c("c_N", "c_J", "naup_catch", "f_J", "f_N", "h_J", "h_N")], 0.98)
    #rand_pars[c("comp_b", "comp_d", "comp_m", "cann")] = 10^runif(4, min=-6, max=-2.5)
    LL_check = full_LL(rand_pars)
    attempts = attempts +1
  }
  print(paste0("parameter likelihood: ", LL_check, " in ", attempts, " attempts"))
  rand_pars
}

var_matrix = function(pars){
  diag(length(pars))*1e-4*pars^2
}

MCMC.parallel2 <- function (p, n, init, n.chain = 4, n.cpu, packages = NULL, dyn.libs = NULL, 
                            scale = rep(1, length(init)), adapt = !is.null(acc.rate), 
                            acc.rate = NULL, gamma = 2/3, list = TRUE, ...) 
{
  cl <- makeCluster(min(n.cpu, detectCores()))
  on.exit({
    stopCluster(cl)
    print("Cluster stopped.")
  })
  varlist <- unique(c(ls(), ls(envir = .GlobalEnv), ls(envir = parent.env(environment()))))
  clusterExport(cl, varlist = varlist, envir = environment())
  clusterSetRNGStream(cl)
  wd <- getwd()
  clusterExport(cl, varlist = c("packages", "dyn.libs", "wd"), 
                envir = environment())
  MCMC.wrap <- function(x, ...) {
    if (!is.null(packages)) 
      sapply(packages, function(x) require(x, character.only = TRUE))
    if (!is.null(dyn.libs)) {
      sapply(dyn.libs, function(x) dyn.load(paste(wd, x, 
                                                  sep = "/")))
      on.exit(sapply(dyn.libs, function(x) dyn.unload(paste(wd, 
                                                            x, sep = "/"))))
    }
    adaptMCMC::MCMC(...)
  }
  
  # list of dispersed starting parameters
  pars_list = c(list(init), replicate(n.chain-1, generate_pars_list(init), simplify = F))
  
  var_matrix = function(pars){
    variance = diag(length(pars))*1e-4*pars^2
  }
  variance_list = lapply(pars_list, FUN=var_matrix)

  result <- clusterMap(cl=cl, fun=MCMC.wrap, init = pars_list,  scale=variance_list, 
                       MoreArgs = list(p = p, n = n,adapt = adapt, acc.rate = acc.rate,
                                       gamma = gamma), ...)
  return(result)
}

################## Runs MCMC chains in parallel ###############
a = Sys.time()
print(paste0("starting MCMC at ", a))

p2 = MCMC.parallel2(p=full_LL, n=250000, init=pars, scale=variances, adapt=50000, n.chain=8, n.cpu = 8, 
                    packages=c("panelPomp"), acc.rate=0.3, list=T)
b = Sys.time()

b-a



saveRDS(p2, file = "Joint_GW_full_A2.RDA")
##############################################################
