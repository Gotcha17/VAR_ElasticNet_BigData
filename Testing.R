### !!! Specify path of the xlsx data file here !!! ###
raw_data_file <- "D:/Dropbox/Uni/Kiel/Master/Semester5/Seminar_EconStats/SW_Updated.xlsx"


### Install packages
#install.packages("readxl"); install.packages("zoo"); install.packages("tseries"); install.packages("zoo"); install.packages("gcdnet"); install.packages("sparsevar"); install.packages("caret"); install.packages("glmnet"); install.packages("mgm")
#install.packages("doParallel"); install.packages("klaR")
install.packages("MSBVAR")
library("MSBVAR")
install.packages("devtools")
library(devtools)
install_github("jeffwong/fastVAR")
library(fastVAR)
devtools::install_github('topepo/caret/pkg/caret')
install.packages("elasticnet")

### Load packages
library(readxl); library(zoo); library(tseries); library(vars); library(gcdnet); library(sparsevar) 
library(glmnet); library(mgm);  library(klaR); library(MASS)
library(elasticnet)

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


### Function for stationarity for p-th lag with Augmented Dickey-Fuller test
aug_dick_fuller <- function(x, p=7, crit=0.05){
  var_number <- NCOL(x)
  adf_result <- c()
  for (i in 1:var_number){
    adf_result[i] <- adf.test(na.omit(x[,i]), k=p)$p.value                                # From 23rd, because of NAs
  }
  non_stat <- names(variables_diff_quart[0, which(adf_result>=crit)])
  return(non_stat)                                                                        # List of non-stationary variables       
}


### Most non-stationary variables are in currency units, thus
### log differences should be used instead of simple diffs
non_stat <- aug_dick_fuller(variables_diff_quart) 
non_stat


### Converting variables expressed in currency units in log diffs that are non-stationary
curr_var_names_1 <- names(read_excel(raw_data_file, range = "FRED!DK1:EB1"))              # Var names of variables in currency units
curr_var_names_2 <- names(read_excel(raw_data_file, range = "FRED!ED1:EP1"))              # Var names of variables in currency units
curr_var_names_all <- c(curr_var_names_1, curr_var_names_2)                               # Glueing together name vectors

# Only non-stationary variables and variables expressed in currency units with all values GTR zero are selected
while (length(non_stat)>0){  
  curr_var_names <- c()
  non_stat_diff <- c()
  for (i in non_stat){
    if (any(i==curr_var_names_all[])==TRUE){
      curr_var_names <- curr_var_names[curr_var_names != i]
      if (all(variables_ts[,i]>0)){
        curr_var_names <- c(curr_var_names, i)
      }
    }else{
      non_stat_diff <- c(non_stat_diff, i)
    }
  }
  # Selected non-stationary diffs of variables are replaced with log diff of these variables
  variables_diff_quart[,curr_var_names] <- aggregate(diff(log(variables_ts[,curr_var_names])), as.yearqtr, frequency = 4)[-1,]
  # Selected non-stationary diffs of variables are replaced diffs of diffs of these variables
  variables_diff_quart[,non_stat_diff] <- diff(variables_diff_quart[,non_stat_diff], na.pad = TRUE)
  # After adjustment of variables ADF Test is run
  non_stat_tmp <- non_stat
  non_stat <- aug_dick_fuller(variables_diff_quart)
  # If no more variables can be adjusted, while loop is broken
  if (all(non_stat_tmp[]==non_stat[])){
    if (length(non_stat_diff)==0){
      break
    }
  }
}


### GDP time series and other variables are "glued" together
end_var <- merge(gdp_logdiff_quart, variables_diff_quart)

# All NAs are removed from final set of variables
end_var <- na.omit(end_var)

# Endogenous variables are split into a train and test sample
end_test <- end_var[(NROW(end_var)-67):NROW(end_var),]
end_train <- end_var[1:(NROW(end_var)-68),]

# Endogenous test and train variables are split into gdp and other variables
gdp_test <- end_test[,1]
var_test <- end_test[,2:NCOL(end_test)]
gdp_train <- end_train[,1:2]
var_train <- end_train[,2:NCOL(end_train)]


### selecting optimal variables for modeling
myTimeControl <- trainControl(method = "timeslice",
                              initialWindow = NROW(end_train),
                              horizon = 8,
                              fixedWindow = FALSE,
                              allowParallel = TRUE)

steplda.mod <- train(GDP ~ ., data = end_train,
                    method = 'lmStepAIC',
                    trControl = myTimeControl)

steplda <- stepclass(GDP ~ ., data = end_train, method = "lda", maxvar = 30, direction = "forward")
steplda.vars <- c("DDURRG3M086SBEA", "DNDGRG3M086SBEA", "DPCERA3M086SBEA", "DSERRG3M086SBEA", "PCEPI",
                  "RPI", "W875RX1", "INDPRO", "IPFINAL", "IPCONGD", "IPDCONGD", "IPNCONGD", "IPBUSEQ",
                  "IPMAT", "IPDMAT", "IPNMAT", "IPFPNSS", "IPFUELN", "TCU", "MCUMFN", "CLF16OV", "CE16OV",
                  "UNRATE", "UEMPMEAN", "UEMPLT5", "UEMP5TO14", "UEMP15OV", "UEMP15T26", "UEMP27OV", "PAYEMS")


### fitVar elastic net
fit.var.enet <- fitVAR(as.matrix(merge(gdp_train[,1], var_train[,25:47])), p=3, penalty = "ENET",
                       method = "timeSlice", scale = TRUE, alpha = 0.5,
                       parallel = TRUE, ncores = 4,
                       type.measure = "mse", horizon = 1, leaveOut = 10)

forc.var.enet <- computeForecasts(fit.var.enet, 8)
View(forc.var.enet)



### gcdnet elastic net
fit.gcd <- gcdnet(var_train, gdp_train, method = "ls", dfmax = 10)

### mgm mvar elastic net
fit.mvar <- mvar(as.matrix(end_train[,1:10]), type = rep("g", NCOL(end_train[,1:10])),
                 consec = 1:NROW(end_train[,1:10]), level = rep(1, NCOL(end_train[,1:10])),
                 lambdaSel = "EBIC", alphaSeq = 0.5, lags = 1:5, method = "glm", signInfo = FALSE,
                 scale = FALSE, saveData = TRUE)

pred.mvar <- predict(fit.mvar, data = as.matrix(end_test[,1:10]), consec = 1:NROW(end_test[,1:10]))


################################## caret
caretdata <- VAR.Z(as.matrix(na.omit(end_train[econ.vars.5])), p = 7, intercept = FALSE)
caret.y <- as.matrix(caretdata$y.p)
caret.x <- as.matrix(caretdata$Z)

end_var_cleaned_matrix <- as.matrix(na.omit(end_var))
gdp_cleaned_matrix <- as.matrix(end_var_cleaned_matrix[,1])
y = as.vector(gdp_cleaned_matrix[-1,])
x = end_var_cleaned_matrix[1:(nrow(end_var_cleaned_matrix)-1),-1]



library(caret);
library(doParallel);
clcaret <- makeCluster(7) # create a cluster with x=? cores
registerDoParallel(cores=clcaret)

myTimeControl <- trainControl(method = "timeSlice",
                              initialWindow = 159,
                              horizon = 1,
                              fixedWindow = FALSE,
                              allowParallel = TRUE)

glmFitTime <- train(x, y,
                    method = "glmnet",
                    family = "gaussian",
                    trControl = myTimeControl,
                    standardize = TRUE, maxit = 1000,
                    intercept = FALSE,
                    tuneGrid = expand.grid(.alpha = 0.5,
                                           .lambda = 10^seq(2,-2,length = 100)))

var.tune = apply(caret.y, 2, train, x=caret.x, method = "glmnet",
                 family = "gaussian", trControl = myTimeControl,
                 standardize = TRUE, maxit = 1000000,
                 tuneGrid = expand.grid(.alpha = seq(0, 1, by = 0.05),
                                        .lambda = 10^seq(2,-2,length = 100)))

var.elnet.tune = apply(caret.y, 2, train, x=caret.x, method = "enet",
                 family = "gaussian", trControl = myTimeControl,
                 preProc = c("center", "scale"))

stopCluster(clcaret)

best.tune = data.frame(lapply(var.tune, function(model) {
  t(model$bestTune)}))

source("D:/Dropbox/Uni/Kiel/Master/Semester5/Seminar_EconStats/R_Script/enetVAR/enetVAR.R")
y.matrix <- as.matrix(na.omit(end_var[,econ.vars.1]))

fit.gVAR <- enetVAR(y.matrix, p = 2, alpha = 0.25, intercept = FALSE)
pred.gVAR.h1 <- predict(fit.gVAR, n.ahead = 12, lambda = best.tune[2,])
View(pred.gVAR.h1)

fit.GroupElnetVAR <- GroupElnetVAR(y.matrix, p = 3, alpha = 0.4, intercept = FALSE)
pred.GroupElnetVAR.h1 <- predict(fit.GroupElnetVAR, n.ahead = 12)
View(pred.GroupElnetVAR.h1)


# Out of sample expirement ------------------------------------------------
ezlasso.all <- ezlasso(end_var, "GDP", alpha = 0.9, maxnrvar = 30)
ezlasso.all <- c("GDP", "DDURRG3M086SBEA", "DPCERA3M086SBEA", "W875RX1", "INDPRO",
                 "IPFINAL", "IPCONGD", "IPDCONGD", "IPNCONGD", "IPBUSEQ", "IPMAT",
                 "IPDMAT", "IPNMAT", "IPFPNSS", "IPFUELN", "TCU", "MCUMFN", "CLF16OV",
                 "CE16OV", "UNRATE", "UEMPMEAN", "UEMPLT5", "UEMP5TO14", "UEMP15OV",
                 "UEMP15T26", "UEMP27OV", "PAYEMS", "USPRIV", "CES1021000001", "FEDFUNDS")

y <- na.omit(end_var[,ezlasso.all])

# import tune for alpha and lambda
import_tune_file <- "D:/Dropbox/Uni/Kiel/Master/Semester5/Seminar_EconStats/R_Script/Best_tune_ezlasso_preselection_30_K_5_p.xlsx"
import_tune <- as.matrix(read_excel(import_tune_file, range = "Tabelle1!A1:AE4")[-1,-1])
import.tune <- matrix(as.numeric(import_tune), nrow = 2, ncol = ncol(import_tune))

result1 <- modeltrain(data = y, lag = 5, alpha = import.tune[1,], lambda = import.tune[2,], const = TRUE)
#result1 mse sum: 0.015914531  <-- without prior tuning of alpha and lambda
#result1 mse sum: 0.009757152  <-- with prior tuning of alpha and lambda
result2 <- modeltrain(data = y, lag = 5, const = TRUE)
#result2 mse sum: 0.00907175092168175

result3 <- modeltrain1(data = y, lag = 1, const = TRUE)
#result3 mse sum: 0.0107635485836803

testing.sample <- na.omit(end_var[,econ.vars.3])
result4 <- modeltrain(data = testing.sample, lag = 3, const = TRUE)
testing1.sample <- na.omit(end_var[,econ.vars.2])
result4.1 <- modeltrain(data = testing1.sample, alpha = 0.95, lag = 4, const = TRUE)
#result4 mse sum: 0.01482406

result5 <- modeltrain(data = y, lag = 5, alpha = 0.5, const = TRUE)
#result5 mse sum: 0.009076151


# Large Testing experiment ------------------------------------------------
### Selecting model to test

# Choosing econ variables
econ.vars.all <- c("GDP", "FEDFUNDS", "DPCERA3M086SBEA", "AWHI", "RPI", "CPIAUCSL", "TB3MS", "GS5",
                   "GS10", "M2SL", "SP500", "MCUMFN", "INDPRO", "UNRATE", "HOUST", "PPIACO", 
                   "PCEPI", "CES3000000008", "M1SL", "WTISPLC")
econ.vars.1 <- c("GDP", "DPCERA3M086SBEA")
econ.vars.2 <- c("GDP", "FEDFUNDS", "CPIAUCSL")
econ.vars.3 <- c("GDP", "DPCERA3M086SBEA", "CPIAUCSL", "TB3MS")
econ.vars.4 <- c("GDP", "DPCERA3M086SBEA", "FEDFUNDS", "AWHI", "RPI")
econ.vars.5 <- c("GDP", "FEDFUNDS", "DPCERA3M086SBEA", "AWHI", "RPI", "GS5", "GS10", "M2SL",
                 "SP500", "MCUMFN", "INDPRO", "UNRATE", "HOUST", "PPIACO", "PCEPI",
                 "CES3000000008", "M1SL", "WTISPLC")
econ.vars <- list(econ.vars.1, econ.vars.2, econ.vars.3, econ.vars.4, econ.vars.5)

# Variables choosen on the basis of autocorrelations
acf.selc.25 <- acf.var.selection2(y = end_train, lag = 6, maxnrvar = 25)
acf.selc.20 <- acf.var.selection2(y = end_train, lag = 7, maxnrvar = 20)
acf.selc.15 <- acf.var.selection2(y = end_train, lag = 8, maxnrvar = 15)
acf.selc.10 <- acf.var.selection2(y = end_train, lag = 9, maxnrvar = 10)
acf.selc.5 <- acf.var.selection2(y = end_train, lag = 10, maxnrvar = 5)

acf.test.vars <- list(acf.selc.5, acf.selc.10, acf.selc.15, acf.selc.20, acf.selc.25)


# Variables choosen on the basis of partial autocorrelations
pacf.selc.25 <- pacf.var.selection(y = end_train, lag = 6, maxnrvar = 25)
pacf.selc.20 <- pacf.var.selection(y = end_train, lag = 7, maxnrvar = 20)
pacf.selc.15 <- pacf.var.selection(y = end_train, lag = 8, maxnrvar = 15)
pacf.selc.10 <- pacf.var.selection(y = end_train, lag = 8, maxnrvar = 10)
pacf.selc.5 <- pacf.var.selection(y = end_train, lag = 8, maxnrvar = 5)

pacf.test.vars <- list(pacf.selc.5, pacf.selc.10, pacf.selc.15, pacf.selc.20, pacf.selc.25)

# All models combined, that have to be tested
test.vars <- c(econ.vars, acf.test.vars, pacf.test.vars)
test.vars <- list(pacf.selc.5, pacf.selc.10)


require(foreach)
cl1 <- makeCluster(6) # create a cluster with x=? cores
registerDoParallel(cl1) # register the cluster
for (m in c("FALSE", "TRUE")){
  model <- 0
  data <- na.omit(end_var)
  lag.min <- 20
  lag.max <- 24
  model <- model+1
  ### Elastic Net
  res1 <- foreach(i = seq(0.25, 0.75, by = 0.25),
                  .packages = c("glmnet",
                                "zoo")) %do% {
                                  resb <- foreach(j = seq(lag.min, lag.max, by = 1),
                                                  .packages = c("glmnet",
                                                                "zoo")) %dopar% {
                                                                  res2 <- mean(unlist((modeltrain.slim(data = data, lag = j, alpha = i,
                                                                                       const = as.logical(m)))$var.msfe))
                                                    
                                                                }
                                  
                                }
  res10 <- matrix(unlist(res1), ncol = 3, nrow = 5)
  colnames(res10) <- c(paste(paste("alpha = ", as.character(seq(0.25, 0.75, by = 0.25))), sep = ","))
  rownames(res10) <- c(paste(paste("p = ", as.character(seq(lag.min, lag.max))), sep = ","))
  result.file.path <- gsub(" ", "", (paste(getwd(), "/result_all_variables_",
                                           as.character(ncol(data)), "_K_", as.character(lag.max),
                                           "_p_", "model_",as.character(model),"_", 
                                           m, ".csv")), fixed = TRUE)
  write.csv(res10, file = result.file.path)
}

stopCluster(cl1)

############## Large testing with list elements
for (m in c("TRUE", "FALSE")){
  model <- 0
  for (i in test.vars){
    data <- na.omit(end_var[,i])
    data.matrix <- as.matrix(na.omit(end_train[,i]))
    lag.values <- enetVARselect(data.matrix, intercept = as.logical(m))
    lag.min <- min(lag.values$IC_lag[2:4])
    lag.max <- max(lag.values$IC_lag[2:4])
    model <- model+1
    ### Elastic Net
    cores <- 3
    cl1 <- makeCluster(cores) # create a cluster with x=? cores
    registerDoParallel(cl1) # register the cluster
    res1 <- foreach(i = seq(0.25, 0.75, by = 0.25),
                    .export= c("theils_u_ar1", "ar1_train"),
                    .packages = c("glmnet",
                                  "zoo")) %dopar% {
                                    resb <- foreach(j = c(lag.min, lag.max),
                                                    .packages = c("glmnet",
                                                                  "zoo")) %do% {
                                                                    result <- modeltrain(data = data, lag = j, alpha = i,
                                                                                         const = as.logical(m))
                                                                  }
                                    names(resb) <- c(as.character(lag.min), as.character(lag.min))
                                  }
    names(res1) <- c("alpha = 0.25", "alpha = 0.50", "alpha = 0.75")
    stopCluster(cl1)
    results[[model]] <- res1
  }
  names(results) <- c("econ.vars.1", "econ.vars.2", "econ.vars.3", "econ.vars.4",
                      "econ.vars.5", "acf.selc.5", "acf.selc.10", "acf.selc.15",
                      "acf.selc.20", "acf.selc.25", "pacf.selc.5", "pacf.selc.10",
                      "pacf.selc.15", "pacf.selc.20", "pacf.selc.25")
  result.f[[(2-1*as.logical(m))]] <- results
}

a<-lapply(result, lapply, lapply, sapply, function(x){mean(unlist(x$var.msfe))})


library(forecast)
accuracy(seq(1,10,by=0.5), seq(0,8,by=0.5), f)

########### fitting arma
ar1.fit <- arima(na.omit(gdp_logdiff_quart), order = c(1,0,0), transform.pars = FALSE)
ar1.pred <- predict(ar1.fit, n.ahead = 8)
ar1.pred.values <- ar1.pred[pred.ind]
ar1.err <- cbind(ar1.err, (ar1.pred-true.values))


############# DM test
DM.test.result <- DMtest(d = (result4.1$for.err$h8^2-result4$for.err$h8^2), l = 8)

############# CW test
y1.pred.h1 <- result4$for.err$h8+as.vector(gdp_test)
y2.pred.h1 <- result4.1$for.err$h8+as.vector(gdp_test)
CW.test.result <- CW_test(result4$for.err$h8, result4.1$for.err$h8, y1.pred.h1, y2.pred.h1, 8)


######## residuals of glmnet
y <- na.omit(end_var[,econ.vars.1])
y.matrix <- as.matrix(y)
### enetVAR Elastic net
fit.gVAR <- enetVAR(y.matrix, p = 2, alpha = 0.25, intercept = FALSE)
residuals <- residuals.enetVAR(fit.gVAR)
Box.test(residuals[,1],lag=3, fitdf=0, type="Lj")

install.packages("portes")
library(portes)
# H0: no autocorrelation
LjungBox(residuals, lags = seq(6,18,3), order = 3)



# CSV File Result Reading -------------------------------------------------

### csv file read and format
res3_1 <- read.csv(gsub(" ", "", (paste(getwd(), "/acf_variables_experiment/result_Enet_ezlasso_variables_",
                                        as.character(23), "_K_", as.character(7),
                                        "_p_", "model_",as.character(5),"_", 
                                        "enet", ".csv")), fixed = TRUE))
rownames(res3_1) <- res3_1[,1]
res3_1 <- res3_1[,-1]

pacf.var.selection(end_var, lag = 8)
acf.var.selection2(end_var, lag = 8)

############# diff log graphix
variables_ts_quart <- aggregate(variables_ts, as.yearqtr, frequency = 4)
DNDGRG3M086SBEA_diff10 <- diff(variables_ts[1:300,2], differences = 1)
DNDGRG3M086SBEA_log_diff10 <- diff(log(variables_ts[1:300,2]), differences = 1)
log_diff_example <- merge(DNDGRG3M086SBEA_diff10, DNDGRG3M086SBEA_log_diff10)
plot(log_diff_example, screens = c(1,2), plot.type ="multiple",
     col = c("2", "4"), ylim = list(c(-0.5,1),c(-0.05,0.1)),
     main = "DNDGRG3M086SBEA Series", xlab = "Time",
     ylab = c("Diffs", "Diffs of log"))



################# information criteria
### Calculation of information criteria (FPE,AIC,HQ,SC) for Elnet Penalty with adopted dof
infCrit2 <- function(enetVAR){
  T <- NROW(enetVAR$var.z$y.p)
  K <- NCOL(enetVAR$var.z$y.orig)
  m <- enetVAR$var.z$p
  alpha <- enetVAR$alpha
  lambda <- mean(sapply(enetVAR$model, function(x){x$lambda.min}))
  U_hat <- t(residuals.enetVAR(enetVAR))
  Sigma_hat <- (U_hat %*% t(U_hat))/T
  det_Sigma_hat <- max(det(Sigma_hat), 0)
  # y = (Z' kron I_k)*beta+u
  B <- t(coef.enetVAR.enetVAR(enetVAR))                             # B=(nu, A1,...,Ap) coeffs
  beta <- t(as.vector(B))                                         # VAR(p) for LS estimation form
  A <- (beta>0 | beta<0)                                          # Support of beta
  # Degrees of freedom for Elastic Net in accordance to TIBSHIRANI AND TAYLOR (2011)
  dof <- sum(A)
  FPE <- (1 + dof/T)/(1 - dof/T) * det_Sigma_hat
  AIC <- log(det_Sigma_hat) + 2/T * dof
  HQ <- log(det_Sigma_hat) + 2*log(log(T))/T * dof
  SC <- log(det_Sigma_hat) + (log(T))/T * dof
  return(list(FPE = FPE, AIC = AIC, HQ = HQ, SC = SC))
}

enetVARselect2 <- function(data, max.lag = 30, alpha = 0.25, intercept = F){
  iteration <- 0
  result <- c()
  for (p in (1:max.lag)){
    fit.enetVAR <- enetVAR(data, p = p, alpha = alpha, intercept = intercept)
    ic <- infCrit2(fit.enetVAR)
    result <- cbind(result, ic)
    iteration <- iteration + 1
    tic <- matrix(unlist(result), ncol = iteration, nrow = 4)
    if (iteration > 3){
      inf.sum <- sum(tic[,iteration]==(-Inf))
      if (inf.sum>2) break
      ch <- 0
      for (i in 1:4){
        min <- min(tic[i,])
        ch <- ch + (sum(result[i,(iteration-3):iteration]>min)>3)
      }
      if (ch==4) break
    }
  }
  FPE_p <- which(tic[1,]==min(tic[1,]))
  AIC_p <- which(tic[2,]==min(tic[2,]))
  HQ_p <- which(tic[3,]==min(tic[3,]))
  SC_p <- which(tic[4,]==min(tic[4,]))
  colnames(result) <- seq(1,iteration)
  ic_matrix <- matrix(c(FPE_p, AIC_p, HQ_p, SC_p), ncol = 4, nrow = 1)
  colnames(ic_matrix) <- c("FPE", "AIC", "HQ", "SC")
  return(list(IC_lag = ic_matrix, IC_value = result))
}







##########################################

l1 <- function(xs){
  y <- sapply(xs, function(x){(sqrt((1 - sqrt(x^2))^2))})
  return(y)
}

l2 <- function(xs){
  y <- sapply(xs, function(x){sqrt(1 - x^2)})
  return(y)
}

enet <- function(xs, z){
  y <- sapply(xs, function(x){((2 - 2 * x - 2 * z + 4 * x * z -
                    (4 * z^2
                     - 8 * x * z^2
                     + 8 * x^2 * z^2
                     - 16 * x^2 * z^3
                     + 8 * x * z^3 + 4 * x^2 * z^4)^(1 / 2)
                  - 2 * x * z^2) / (2 - 4 * z))})
  return(y)
}


# Contour plot of L1,L2 and Elastic Net penalties -------------------------


####plot cross
cross <- function(ext){
  plot(c(-ext, ext), c(0, 0), type = "l", ylim = c(-3,3), xlim = c(-3,3),
       xlab = expression(beta[1]), ylab = expression(beta[2]),
       main = "LASSO, Ridge and Elastic Net penalties")
  
  lines(c(0, 0), c(-ext, ext), type = "l")
}

xs <- seq(0, 1, length.out = 1000)

alpha <- 0.501

cross(5)

l1_col <- 1
l2_col <- 1
enet_col <- 1
lw <- 2

lines(xs, l1(xs), col = l1_col, lwd = lw, lty = 4)
lines(xs, (-l1(xs)), col = l1_col, lwd = lw, lty = 4)
lines(-xs, l1(xs), col = l1_col, lwd = lw, lty = 4)
lines(-xs, (-l1(xs)), col = l1_col, lwd = lw, lty = 4)

lines(xs, l2(xs), col = l2_col, lwd = lw, lty = 3)
lines(xs, (-l2(xs)), col = l2_col, lwd = lw, lty = 3)
lines(-xs, l2(xs), col = l2_col, lwd = lw, lty = 3)
lines(-xs, (-l2(xs)), col = l2_col, lwd = lw, lty = 3)

lines(xs, enet(xs, alpha), col = enet_col, lwd = lw, lty = 1)
lines(xs, -enet(xs, alpha), col = enet_col, lwd = lw, lty = 1)
lines(-xs, enet(xs, alpha), col = enet_col, lwd = lw, lty = 1)
lines(-xs, -enet(xs, alpha), col = enet_col, lwd = lw, lty = 1)

legend(x = "topright", legend = c("LASSO", "Ridge", "Elastic Net"), lty = c(4,3,1))


# Elipse centered at (x0, y0) with axes _a_ and _b_:
a <- 1; b <- 0.5; alpha = 1
theta <- seq(0, 2 * pi, length=(1000))
x <- 0.7 + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
y <- 1.8 + a * cos(theta) * sin(alpha) + b * sin(theta) * cos(alpha)
lines(x, y, type = "l")



################ Unlisting values



### Unlisting values
forc.h1 <- cbind(sapply(results, sapply, function(x){
  return(c(x$var.msfe$h1, x$model$var.z$p)) 
}))

forc.h1.U <- cbind(sapply(results, sapply, function(x){
  return(c(x$theils_u_ar1$h1, x$model$var.z$p)) 
}))

forc.h2 <- cbind(sapply(results, sapply, function(x){
  return(c(x$var.msfe$h2, x$model$var.z$p)) 
}))

forc.h2.U <- cbind(sapply(results, sapply, function(x){
  return(c(x$theils_u_ar1$h2, x$model$var.z$p)) 
}))

forc.h4 <- cbind(sapply(results, sapply, function(x){
  return(c(x$var.msfe$h4, x$model$var.z$p)) 
}))

forc.h4.U <- cbind(sapply(results, sapply, function(x){
  return(c(x$theils_u_ar1$h4, x$model$var.z$p)) 
}))

forc.h8 <- cbind(sapply(results, sapply, function(x){
  return(c(x$var.msfe$h8, x$model$var.z$p)) 
}))

forc.h8.U <- cbind(sapply(results, sapply, function(x){
  return(c(x$theils_u_ar1$h8, x$model$var.z$p)) 
}))

View(rbind(forc.h1, forc.h2, forc.h4, forc.h8))
View(rbind(forc.h1.U, forc.h2.U, forc.h4.U, forc.h8.U))


# Further model testing ---------------------------------------------------
y <- na.omit(end_train[,acf.selc.5])
y.matrix <- as.matrix(y)
### enetVAR Elastic net
fit.gVAR <- enetVAR(y.matrix, p = 2, alpha = 0.25, intercept = FALSE)
pred.gVAR.h1 <- predict(fit.gVAR, n.ahead = 8)
View(pred.gVAR.h1)

### enetVAR Group Elastic net
fit.GroupElnetVAR <- GroupEnetVAR(y.matrix, p = 2, alpha = 0.25, intercept = FALSE, parallel = FALSE)
pred.GroupElnetVAR.h1 <- predict(fit.GroupElnetVAR, n.ahead = 8)
View(pred.GroupElnetVAR.h1)


### Elastic Net 
require(foreach)
cl1 <- makeCluster(5) # create a cluster with x=? cores
registerDoParallel(cl1) # register the cluster
res1 <- foreach(i = seq(0.05, 0.95, length.out = 5),
                .combine = "rbind",
                .packages = c("glmnet",
                              "zoo")) %dopar% {
                                resb <- foreach(j = seq(8, 10, by = 1),
                                                .combine = "cbind",
                                                .packages = c("glmnet",
                                                              "zoo")) %do% {
                                                                result <- modeltrain(data = y, lag = j, alpha = i,
                                                                                     const = TRUE)
                                                              }
                              }
stopCluster(cl1)
rownames(res1) <- c(paste(paste("alpha = ", as.character(seq(0.05, 0.95, length.out = 5))), sep = ","))
colnames(res1) <- c(paste(paste("p = ", as.character(seq(4, 5, by = 1))), sep = ","))
result.file.path <- gsub(" ", "", (paste(getwd(), "/result_Enet_econ_variables_",
                                         as.character(ncol(data)), "_K_", as.character(lag.max),
                                         "_p_", "model_",as.character(model),"_", 
                                         m, ".csv")), fixed = TRUE)
write.csv(res1, file = result.file.path)


####### Plotting forecasts
# Plotting Forecasted Series in absolute GDP values
plot(zoo(diff_log2norm(results$acf.selc.15$`1`$forecast$h1, 12323.336),
         index(gdp_test)), type = "l", ylab = "US GDP in Billions of Chained 2009 Dollars",
         xlab = "Time", main = "One step ahead forecast comparison")
# Adding True TS in absolute GDP values
lines(zoo(raw_gdp[(NROW(raw_gdp)-67):NROW(raw_gdp),2], index(gdp_test)), col = 2)
# Adding AR1 TS in absolute GDP values
ar1_ts <- ar1_train(na.omit(gdp_logdiff_quart))
lines(zoo(diff_log2norm(ar1_ts$forecast$h1, 12323.336), index(gdp_test)), col = 3)
# Adding a legend
legend(x = "topleft", legend = c("ACF.Selc.15, lag=1", "True", "AR1"), lty = c(1,1,1), col = c(1,2,3))

# Plotting Forecasted Series errors for h=1
plot(zoo(results2$acf.selc.10$`7`$for.err$h1, index(gdp_test)),
     type = "h", lwd = 8, col = 1,
     ylab = "Forecast errors in log differences",
     xlab = "Time", main = "One step ahead forecast error comparison for US GDP")
# Adding another forecast
lines(zoo(results$acf.selc.15$`2`$for.err$h1, index(gdp_test)), col = 6, lwd = 2, type = "h")
# Adding a legend
legend(x = "topleft", legend = c("ACF.Selc.10, lag=7", "ACF.Selc.15, lag=2"),
       lty = c(1,1), col = c(1,6), lwd = c(8,2))

# Plotting Forecasted Series errors for h=8
plot(zoo(results2$acf.selc.10$`7`$for.err$h8, index(gdp_test)),
     type = "h", lwd = 8, col = 1,
     ylab = "Forecast errors in log differences",
     xlab = "Time", main = "Eight steps ahead forecast error comparison for US GDP")
# Adding another forecast
lines(zoo(results$acf.selc.15$`2`$for.err$h8, index(gdp_test)), col = 6, lwd = 2, type = "h")
# Adding a legend
legend(x = "topleft", legend = c("ACF.Selc.10, lag=7", "ACF.Selc.15, lag=2"),
       lty = c(1,1), col = c(1,6), lwd = c(8,2))

# Plotting Forecasted Series in absolute GDP values h=1
plot(zoo(diff_log2norm(results2$acf.selc.10$`7`$forecast$h1, 12323.336),
         index(gdp_test)), type = "l", lwd = 2, ylim = c(12000, 17000),
     ylab = "US GDP in Billions of Chained 2009 Dollars",
     xlab = "Time", main = "One step ahead forecast comparison")
# Adding AR1 TS in absolute GDP values
lines(zoo(diff_log2norm(results$acf.selc.15$`2`$forecast$h1, 12323.336), index(gdp_test)), col = 3, lwd = 2)
# Adding AR1 TS in absolute GDP values
lines(zoo(diff_log2norm(results3$enet.selc.25$`3`$forecast$h1, 12323.336), index(gdp_test)), col = 4, lwd = 2)
# Adding True TS in absolute GDP values
lines(zoo(raw_gdp[(NROW(raw_gdp)-67):NROW(raw_gdp),2], index(gdp_test)), col = 2, lwd = 2)
# Adding a legend
legend(x = "topleft", legend = c("ACF.Selc.10, lag=7", "ACF.Selc.15, lag=2", "Econ.Selc.25, lag=3", "True"),
       lty = c(1,1,1), col = c(1,3,4,2), lwd = c(2,2,2))

# Plotting Forecasted Series in absolute GDP values h=8
plot(zoo(diff_log2norm(results2$acf.selc.10$`7`$forecast$h8, 12323.336),
         index(gdp_test)), type = "l", ylim = c(12000, 17000),
     lwd = 2, ylab = "US GDP in Billions of Chained 2009 Dollars",
     xlab = "Time", main = "Eight steps ahead forecast comparison")
# Adding AR1 TS in absolute GDP values
lines(zoo(diff_log2norm(results$acf.selc.15$`2`$forecast$h8, 12323.336), index(gdp_test)), col = 3, lwd = 2)
# Adding True TS in absolute GDP values
lines(zoo(raw_gdp[(NROW(raw_gdp)-67):NROW(raw_gdp),2], index(gdp_test)), col = 2, lwd = 2)
# Adding a legend
legend(x = "topleft", legend = c("ACF.Selc.10, lag=7", "ACF.Selc.15, lag=2", "True"),
       lty = c(1,1,1), col = c(1,3,2), lwd = c(2,2,2))
