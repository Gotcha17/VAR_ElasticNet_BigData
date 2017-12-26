#################################################################################################
#            Seminar in Econometrics: "Forecasting in data-rich environment"                    #
#          Quartely US GDP Forecast Experiment with elastic net shrinkage method                #
#                             Student: Mykhaylo Cherner                                         #
#                             Student-ID: 1010152                                               #
#################################################################################################


# Data preparation --------------------------------------------------------


### !!! Specify path of the xlsx data file here !!! ###
raw_data_file <- "D:/Dropbox/Uni/Kiel/Master/Semester5/Seminar_EconStats/SW_Updated.xlsx"


### Install packages
install.packages("readxl"); install.packages("zoo"); install.packages("tseries")
install.packages("doParallel"); install.packages("portes"); install.packages("caret"); 
install.packages("glmnet")

# Caret package has to be reinstalled if parallel option is desired
#devtools::install_github('topepo/caret/pkg/caret')


### Load packages
source("D:/Dropbox/Uni/Kiel/Master/Semester5/Seminar_EconStats/R_Script/enetVAR/enetVAR.R")
library(readxl); library(zoo); library(tseries); library(glmnet);
library(doParallel); library(caret); library(portes)


### Read column names from xlsx file
col_names <- names(read_excel(raw_data_file, range = "FRED!A1:EQ1"))                      # Column names


### Read data from xlsx file
raw_variables <- read_excel(raw_data_file, range = "'Monthly Data'!A12:EQ707", col_names = col_names)
raw_gdp <- read_excel(raw_data_file, range = "'US GDP'!A56:B287", col_names = c("date", "GDP"))

### Calculate first differences of the variables and aggregate monthly data to quartely
variables_ts <- zoo(raw_variables[,-1], raw_variables$`1date`)                            # Creates TS object

# Aggregation of monthly difference to quarters without Q1 1959 and replacing missing values with aggregated ones
variables_diff_quart <- aggregate(diff(variables_ts), as.yearqtr, frequency = 4)[-1,]


### Calculate first log differences of quartely GDP data
gdp_ts <- zoo(raw_gdp[,-1], as.yearqtr(strptime(raw_gdp$date,"%Y-%d-%m")))                # Date reformatting was needed
gdp_logdiff_quart <- diff(log(gdp_ts))                                                    # Logs to ensure non-increasing variance over time


### Most non-stationary variables are in currency units, thus
### log differences should be used instead of simple diffs
non_stat <- aug_dick_fuller(variables_diff_quart) 
non_stat


### Converting variables expressed in currency units in log diffs that are non-stationary
curr_var_names_1 <- names(read_excel(raw_data_file, range = "FRED!DK1:EB1"))              # Var names of variables in currency units
curr_var_names_2 <- names(read_excel(raw_data_file, range = "FRED!ED1:EP1"))              # Var names of variables in currency units
curr_var_names_3 <- names(read_excel(raw_data_file, range = "FRED!CW1:DA1"))              # Var names of variables in currency units
curr_var_names_all <- c(curr_var_names_1, curr_var_names_2, curr_var_names_3)             # Glueing together name vectors

# Only non-stationary variables and variables expressed in currency units with all values GTR zero are selected
curr_var_names_inf <- c()
non_stat_diff_inf <- c()
while (length(non_stat)>0){  
  curr_var_names <- c()
  non_stat_diff <- c()
  for (i in non_stat){
    if (any(i==curr_var_names_all[])==TRUE){
      curr_var_names_all <- curr_var_names_all[curr_var_names_all != i]
      if (all(na.omit(variables_ts[,i])>0)==TRUE){
        curr_var_names <- c(curr_var_names, i)
        # To keep track of which variables will be diff-loged
        curr_var_names_inf <- c(curr_var_names_inf, i)
      }else{
        non_stat_diff <- c(non_stat_diff, i)
      }
    }else{
      non_stat_diff <- c(non_stat_diff, i)
    }
    # To keep track of which variables will be differentiated
  }
  curr_var_names_inf <- c(curr_var_names_inf, curr_var_names)
  non_stat_diff_inf <- c(non_stat_diff_inf, non_stat_diff)
  # Selected non-stationary diffs of variables are replaced with log diff of these variables
  variables_diff_quart[,curr_var_names] <- aggregate(diff(log(variables_ts[,curr_var_names])), as.yearqtr, frequency = 4)[-1,]
  # Selected non-stationary diffs of variables are replaced diffs of diffs of these variables
  variables_diff_quart[,non_stat_diff] <- diff(variables_diff_quart[,non_stat_diff], na.pad = TRUE)
  # After adjustment of variables ADF Test is run
  non_stat <- aug_dick_fuller(variables_diff_quart)
}
# <-- --> "NONBORRES" should probably be exluded

### GDP time series and other variables are "glued" together
end_var <- merge(gdp_logdiff_quart, variables_diff_quart)


# Endogenous variables are split into a train and test sample
end_test <- end_var[(NROW(end_var)-67):NROW(end_var),]
end_train <- end_var[1:(NROW(end_var)-68),]

# Endogenous test and train variables are split into gdp and other variables
gdp_test <- window(gdp_logdiff_quart, start = time(end_test)[1], end = time(end_test)[length(time(end_test))])
var_test <- end_test[,2:NCOL(end_test)]
gdp_train <- window(gdp_logdiff_quart, start = time(end_train)[1], end = time(end_train)[length(time(end_train))])
var_train <- end_train[,2:NCOL(end_train)]


# Model selection ---------------------------------------------------------

### Selecting model to test

# Choosing econ variables
econ.vars.all <- c("GDP", "FEDFUNDS", "DPCERA3M086SBEA", "AWHI", "RPI", "CPIAUCSL", "TB3MS", "GS5",
                   "GS10", "M2SL", "SP500", "MCUMFN", "INDPRO", "UNRATE", "HOUST", "PPIACO", 
                   "PCEPI", "CES3000000008", "M1SL", "WTISPLC")
econ.vars.1 <- c("GDP", "DPCERA3M086SBEA")
econ.vars.2 <- c("GDP", "FEDFUNDS", "CPIAUCSL")
econ.vars.3 <- c("GDP", "DPCERA3M086SBEA", "CPIAUCSL", "TB3MS")
econ.vars.4 <- c("GDP", "DPCERA3M086SBEA", "FEDFUNDS", "AWHI", "RPI")
econ.vars.5 <- c("GDP", "FEDFUNDS", "DPCERA3M086SBEA", "AWHI", "RPI", "GS5",
                 "GS10", "M2SL", "SP500", "MCUMFN", "INDPRO", "UNRATE",
                 "HOUST", "PPIACO", "PCEPI", "CES3000000008", "M1SL", "WTISPLC")
econ.vars <- list(econ.vars.1, econ.vars.2, econ.vars.3, econ.vars.4, econ.vars.5)


# Variables choosen on the basis of autocorrelations
acf.selc.25 <- acf.var.selection2(y = end_train, lag = 20, maxnrvar = 25)
acf.selc.20 <- acf.var.selection2(y = end_train, lag = 20, maxnrvar = 20)
acf.selc.15 <- acf.var.selection2(y = end_train, lag = 20, maxnrvar = 15)
acf.selc.10 <- acf.var.selection2(y = end_train, lag = 20, maxnrvar = 10)
acf.selc.5 <- acf.var.selection2(y = end_train, lag = 20, maxnrvar = 5)

acf.test.vars <- list(acf.selc.5, acf.selc.10, acf.selc.15, acf.selc.20, acf.selc.25)


# Variables choosen on the basis of partial autocorrelations
pacf.selc.25 <- pacf.var.selection(y = end_train, lag = 8, maxnrvar = 25)
pacf.selc.20 <- pacf.var.selection(y = end_train, lag = 8, maxnrvar = 20)
pacf.selc.15 <- pacf.var.selection(y = end_train, lag = 8, maxnrvar = 15)
pacf.selc.10 <- pacf.var.selection(y = end_train, lag = 8, maxnrvar = 10)
pacf.selc.5 <- pacf.var.selection(y = end_train, lag = 8, maxnrvar = 5)

pacf.test.vars <- list(pacf.selc.5, pacf.selc.10, pacf.selc.15, pacf.selc.20, pacf.selc.25)

# !!!CAUTION!!! Takes very long time to run! #

# Variables choosen via the elastic net
#enet.selc.25 <- enetVARpreselection(as.matrix(na.omit(end_train)), maxnrvar = 25, lag = 3)
# Simulations have been run and results below have been obtained
enet.selc.25 <- c("GDP", "AHETPI", "CES0600000008", "CES3000000008", "CES2000000008",
                  "M2SL", "M1SL", "TCDSL", "CURRSL", "LOANINVNSA", "REALLN", "NONREVSL",
                  "MABMM301USM189S", "CUUR0000SAD", "M2REAL", "CUUR0000SEFV",
                  "DDURRG3M086SBEA", "CPIULFSL", "CUSR0000SAS", "INDPRO", "CUUR0000SA0L2",
                  "IPDMAT", "PCEPI", "DSERRG3M086SBEA", "M2MOWN")

#enet.selc.20 <- enetVARpreselection(as.matrix(na.omit(end_train)), maxnrvar = 20, lag = 3)
enet.selc.20 <- c("GDP", "AHETPI", "CES0600000008", "CES3000000008", "CES2000000008",
                  "M2SL", "M1SL", "TCDSL", "CURRSL", "LOANINVNSA", "REALLN", "NONREVSL",
                  "MABMM301USM189S", "CUUR0000SAD", "M2REAL", "CUUR0000SEFV",
                  "DDURRG3M086SBEA", "CPIULFSL", "CUSR0000SAS", "INDPRO")

#enet.selc.15 <- enetVARpreselection(as.matrix(na.omit(end_train)), maxnrvar = 15, lag = 3)
enet.selc.15 <- c("GDP", "AHETPI", "CES0600000008", "CES3000000008", "CES2000000008",
                  "M2SL", "M1SL", "TCDSL", "CURRSL", "LOANINVNSA", "REALLN", "NONREVSL",
                  "MABMM301USM189S", "CUUR0000SAD", "M2REAL")

#enet.selc.10 <- enetVARpreselection(as.matrix(na.omit(end_train)), maxnrvar = 10)
enet.selc.10 <- c("GDP", "AHETPI", "CES0600000008", "CES3000000008", "CES2000000008",
                  "M2SL", "M1SL", "TCDSL", "CURRSL", "LOANINVNSA")

#enet.selc.5 <- enetVARpreselection(as.matrix(na.omit(end_train)), maxnrvar = 5)
enet.selc.5 <- c("GDP", "AHETPI", "CES0600000008", "CES3000000008", "CES2000000008")

enet.test.vars <- list(enet.selc.5, enet.selc.10, enet.selc.15, enet.selc.20, enet.selc.25)

# All models combined, that have to be tested
test.vars <- c(econ.vars, acf.test.vars, pacf.test.vars, enet.test.vars)


# Out of sample expirement ------------------------------------------------

# !!!CAUTION!!! Takes very long time to run! #
require(caret)
require(glmnet)
require(zoo)
require(doParallel)
results <- list()

## Lag order is determined via information criteria and only models for which variables
#  have been preselected via the elastic net. 
model <- 0
for (i in enet.test.vars){
  # Train + Test data
  data <- na.omit(end_var[,i])
  # Train data only
  data.matrix <- as.matrix(na.omit(end_train[,i]))
  # Find optimal lag order based on information criteria
  lag.values <- enetVARselect(data.matrix, intercept = F)
  # FPE is not considered, only AIC, HQ, SC
  lag.min <- min(lag.values$IC_lag[2:4])
  lag.max <- max(lag.values$IC_lag[2:4])
  model <- model+1
  ### Elastic Net
  res1.sum <- list()
  iter <- 0
  for(j in c(lag.min, lag.max)){
    # Get best tune parameters for alpha & lambda based on train set
    iter <- iter + 1
    cores2 <- (detectCores()-1)
    cl2 <- makeCluster(cores2) # create a cluster with x=? cores
    registerDoParallel(cl2) # register the cluster
    tune <- enetVARtune(data.matrix, lag = j, init.window = 40,
                        horizon = 8, parallel = T)
    stopCluster(cl2)
    alpha <- tune[1,]
    lambda <- tune[2,]
    res1 <- modeltrain(data = data, lag = j, alpha = alpha,
                         lambda = lambda, const = F)
    res1.sum[[iter]] <- res1
  }
  names(res1.sum) <- c(as.character(lag.min), as.character(lag.max))
  results[[model]] <- res1.sum
}
names(results) <- c("enet.selc.5", "enet.selc.10", "enet.selc.15", "enet.selc.20",
                    "enet.selc.25")



# !!!CAUTION!!! Takes very long time to run! #
require(caret)
require(glmnet)
require(zoo)
require(doParallel)
results2 <- list()

# Lag order is predetermined based on "economic feel...", all models are tested, except for the
# models that have been preseleted via the elastic net
model <- 0
for (i in test.vars){
  # Train + Test data
  data <- na.omit(end_var[,i])
  # Train data only
  data.matrix <- as.matrix(na.omit(end_train[,i]))
  # FPE is not considered, only AIC, HQ, SC
  lag.min <- floor(24/NCOL(data)^(2/3))
  lag.max <- ceiling(24/NCOL(data)^(2/3)+1)
  model <- model+1
  ### Elastic Net
  res1.sum <- list()
  iter <- 0
  for(j in c(lag.min, lag.max)){
    # Get best tune parameters for alpha & lambda based on train set
    iter <- iter + 1
    cores2 <- (detectCores()-1)
    cl2 <- makeCluster(cores2) # create a cluster with x=? cores
    registerDoParallel(cl2) # register the cluster
    tune <- enetVARtune(data.matrix, lag = j, init.window = 40,
                        horizon = 8, parallel = T)
    stopCluster(cl2)
    alpha <- tune[1,]
    lambda <- tune[2,]
    res1 <- modeltrain(data = data, lag = j, alpha = alpha,
                       lambda = lambda, const = F)
    res1.sum[[iter]] <- res1
  }
  names(res1.sum) <- c(as.character(lag.min), as.character(lag.max))
  results2[[model]] <- res1.sum
}
names(results2) <- c("econ.vars.1", "econ.vars.2", "econ.vars.3", "econ.vars.4",
                    "econ.vars.5", "acf.selc.5", "acf.selc.10", "acf.selc.15",
                    "acf.selc.20", "acf.selc.25", "pacf.selc.5", "pacf.selc.10",
                    "pacf.selc.15", "pacf.selc.20", "pacf.selc.25")

#### Experiment with all variables without preselection. Variable selection is
#    only made with the 'natural' property of the elastic net.
pretune.allvars.enet.p1 <- enetVAR(as.matrix(na.omit(end_train)), p = 1, alpha = 0.75)
lambda1 <- sapply(pretune.allvars.enet.p1$model, function(x){return(x$lambda.min)})
all.vars.enet.p1 <- modeltrain(na.omit(end_var), alpha = 0.75, lambda = lambda1, lag = 1)

pretune.allvars.enet.p2 <- enetVAR(as.matrix(na.omit(end_train)), p = 2, alpha = 0.75)
lambda2 <- sapply(pretune.allvars.enet.p1$model, function(x){return(x$lambda.min)})
all.vars.enet.p2 <- modeltrain(na.omit(end_var), alpha = 0.75, lambda = lambda2, lag = 2)

pretune.allvars.enet.p3 <- enetVAR(as.matrix(na.omit(end_train)), p = 3, alpha = 0.75)
lambda3 <- sapply(pretune.allvars.enet.p3$model, function(x){return(x$lambda.min)})
all.vars.enet.p3 <- modeltrain(na.omit(end_var), alpha = 0.75, lambda = lambda3, lag = 3)

# Model testing -----------------------------------------------------------

## Example with ACF.Selc.15 model for the OOS-Experiment (does not take to long to run)
## Note, alpha and lambda has to be tuned beforehand on the training sample
#  otherwise some alpha has to be provided and lambda is tuned via cross-validation
#  that is adopted to time series data (folds are created in respect to time series)
test.model.acf15.p3 <- modeltrain(na.omit(end_var[,acf.selc.15]), lag = 3,
                                  alpha = 0.25, const = FALSE)

## Testing autocorrelation of residuals
# Calculating residuals for a given model
resids <- test.model.acf15.p3$residuals

# Portmanteau test for autocorrelation
Hosking(resids, order = 3)

# simulate AR(1) to get the for.err. and compute forecasts
ar1 <- ar1_train(na.omit(end_var[,1]))
# Example to get the forecast from  for.err.--> ar1$for.err$h1 + end_var[164:231,1]

## Testing if model forecasts differ systematically
# Clark-West Test for nested models
cw.test <- CW_test(e1 = test.model.acf15.p3$for.err$h2, e2 = ar1$for.err$h2,
                   yf1 = test.model.acf15.p3$forecast$h2, yf2 = (ar1$for.err$h2 + end_var[164:231,1]),
                   nwlag = 2)

# Diebold-Mariano test
dm.test <- DMtest(d = (test.model.acf15.p3$for.err$h2 - ar1$for.err$h2),
                  l = 2)
