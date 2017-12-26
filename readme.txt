Files that have to be available: SW_Updated.xlsx, Main.R, enetVAR.R

SW_Updated.xlsx contains the data; Main.R contains scripts that have to be run; enetVAR.R contains used functions


How to run the scripts:

1. Please open the file: Main.R

2. In the 13-th line change the path to the location of the SW_Updated.xlsx file

3. Lines: 17 - 19 install packages, that are needed to execute the scripts

4. In line 26 you have to specify the path of the enetVAR.R file

5. Execute until line 150, in these lines data is handled. Aggregation, differenciation, stationarity, etc.
   Furthermore, lines 129 - 133 perform the ACF preselection procedure.
   Lines 139 - 143 perform the PACF preselection procedure.

6. Since this procedure takes a rather long time to run, it was executed and results were already
   put in the according variables. You may execute this line to varify the results. This procedure
   performs variables preselection using the VAR model with elastic net and evaluates the performance
   of the models via the SC information criterion

7. You may run the script until line 185. Until this line, all the models have been constructed that have
   to be tested.

8. Lines 185 to 227 contain procedures to get the optimal lag order based on information criteria,
   to tune alpha and lambda parameters on the training sample and then the out-of-sample experiment
   is performed with each model. This is an example of what has been done with all models.
   (On my computer it runs about 6 hours...)
   However, if you wish to run this, you need to go to line 22 and to install another version of the
   {caret}-package, otherwise parallel execution is not supported and the parallel option has to be
   set to FALSE in line 215!

9. Lines 232 to 274 contain basically the same procedures, one main difference is, that the lag order
   is not selection with the information criteria, but is preselected. It was done to test the models
   with higher lags, than the ones obtained with the information criteria.
   (On my computer it runs about 12 hours...)
   However, if you wish to run this, you need to go to line 22 and to install another version of the
   {caret}-package, otherwise parallel execution is not supported and the parallel option has to be
   set to FALSE in line 260!
	

10. Lines 278 to 288 contain procedures to tune lambda for a given alpha and lag order for a model
    containing all variables. Afterwards model is tested with the out-of-sample experiment. This
    is also a rather time consuming procedure

11. To demonstrate a working example with one model, you can run lines 296 - 319.

12. Line 296 starts the modeltrain() function, which performs the out-of-sample experiment for
    predefined alpha and a predefined lag order. If no lambda is given, then lambda values
    are obtained with an internal cross-validation procedure, that is adjusted to time series data.
    Returned result is a list of different objects. There the MSFE (MSPE) and the Theil's U against
    a random walk and agains the AR(1) process are included. 

13. Line 301 extracts the residuals of the aformentioned model

14. Line 304 performs the "Hosking" Portmanteau test for autocorrelation

15. Line 307 trains a AR(1) process with the OOS experiment.

16. Line 312 performs the Clark and West test 

17. Line 317 performs the Diebold-Mariano test


Some further notes on the scripts. Of course more scripts had to be run to obtain all the results, that are
presented in the paper. However, the included scripts have to be modified just slightly (just change the needed
data inputs) to obtain all the results from the paper. Further details can be given.



Some notes on the function included in the enetVAR.R file:

Since it was not possible for me to find a R package, that performs VAR elastic net regression properly,
I had to take an old R {fastVAR}-package (was not maintained for 5 years) and to rewrite some functions.
To some degree these functions are based on the functions from the {fastVAR}-package:

	.enetVAR; enetVAR; coef.enetVAR.enetVAR; predict.enetVAR.enetVAR; VAR.Z; GroupEnetVAR;
	coef.enetVAR.GroupEnetVAR; predict.enetVAR.GroupEnetVAR; 

And since the functions for the Clark and West test and the Diebold-Mariano test were given for MatLab,
I have rewritten them into R code. Including the Newey-West function nw().

All other functions have been written by myself.