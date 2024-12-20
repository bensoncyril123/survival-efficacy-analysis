/* 
** Import the dataset from the specified file path and handle missing values.
** The 'missover' statement ensures that missing values are treated as missing, 
** and 'dsd' handles delimiters correctly.
*/

data Pbc3;
    infile '/home/u63669348/Applied Survival/Final Project/Pbc3.txt' delimiter='09'x firstobs=2 missover dsd;
    
    /* Define the input format for the dataset */
    input 
        SN : $4.
        ptno : best32.
        unit : best32.
        tment : best32.
        sex : best32.
        age : best32.
        stage : $2.
        gibleed : best32.
        crea : best32.    /* Creatinine - treat NA as missing */
        alb : best32.     /* Albumin - treat NA as missing */
        bili : best32.    /* Bilirubin */
        alkph : best32.   /* Alkaline Phosphatase */
        asptr : best32.   /* Aspartate Aminotransferase */
        weight : best32.  /* Weight - treat NA as missing */
        days : best32.    /* Time-to-event variable (Days) */
        status : best32.  /* Status variable indicating if the patient was censored */
    ;

    /* Handle missing values explicitly for specific variables */
    if alb = 'NA' then alb = .;
    if bili = 'NA' then bili = .;
    if weight = 'NA' then weight = .;
run;

/* 
** Sort the data by the 'days' variable (time-to-event) in ascending order. 
** This helps in better analysis of survival times.
*/
proc sort data=Pbc3 out=Pbc3;
    by days;
run;


/* Step 1: Descriptive Statistics for Continuous Variables */

proc means data=Pbc3 n mean median std min max;
   var age crea alb bili alkph asptr weight days;
run;


/* Step 2: Histogram for Time to Failure (Days) */
proc univariate data=Pbc3;
   var days;
   histogram / normal;  /* Add normal curve to histogram */
   inset mean std median min max / position=ne;  /* Display summary statistics */
run;

/* Step 3: Scatter Plots for Relationships Between Continuous Variables and Time to Failure (Days) */
title "Scatter Plots: Relationships Between Age and Time to Failure";
proc sgscatter data=Pbc3;
   plot days*age / grid;
   plot days*weight / grid;
   plot days*crea / grid;
run;

title "Scatter Plots: Relationships Between Weight and Time to Failure";
proc sgscatter data=Pbc3;
   plot days*weight / grid;
run;


/* Step 4: Boxplot for Time to Failure by Treatment Type */
proc sgplot data=Pbc3;
   vbox days / category=tment;
   title "Time to Failure by Treatment Type";
run;


/* Step 5: Pearson Correlation for Continuous Variables */
proc corr data=Pbc3;
   var days age crea alb bili alkph asptr weight;
run;


/* Step 6: Log-Rank Test for Treatment Type (Categorical Variable) */
proc lifetest data=Pbc3;
   time days*status(0);  /* status=0 indicates censored data */
   strata tment;        /* Stratify by treatment type */
   test tment;          /* Perform log-rank test */
run;


/* Step 7: Multicollinearity Check Using VIF (Variance Inflation Factor) */
proc reg data=Pbc3;
   model days = age sex tment stage crea alb bili alkph asptr weight;
   collin;  /* Request multicollinearity diagnostics */
run;


/* Step 8: Survival Analysis by Age Group */

/* Create Age Groups */
data Pbc3;
   set Pbc3;
   if age < 40 then age_group = "Under 40";
   else if 40 <= age < 50 then age_group = "40-49";
   else if 50 <= age < 60 then age_group = "50-59";
   else if 60 <= age < 70 then age_group = "60-69";
   else age_group = "70 and above";
run;

/* Kaplan-Meier Survival Analysis by Age Group */
proc lifetest data=Pbc3 plots=survival(cl);
   time days*status(0); /* Define survival time and censoring status */
   strata age_group;    /* Compare survival across age groups */
   title "Survival Analysis by Age Group using Log-Rank Test";
run;


/* Step 9: Survival Analysis by Weight Group */

/* Create Weight Groups */
data Pbc3;
   set Pbc3;
   if weight < 50 then weight_group = "Under 50";
   else if 50 <= weight < 60 then weight_group = "50-59";
   else if 60 <= weight < 70 then weight_group = "60-69";
   else weight_group = "70 and above";
run;

/* Kaplan-Meier Survival Analysis by Weight Group */
proc lifetest data=Pbc3 plots=survival(cl);
   time days*status(0); /* Define survival time and censoring status */
   strata weight_group; /* Compare survival across weight groups */
   title "Survival Analysis by Weight Group using Log-Rank Test";
run;


/* Step 10: Survival Analysis by Bleeding Status */

/* Kaplan-Meier Survival Curves for Bleeding Status */
proc lifetest data=Pbc3 plots=(survival(cl));
   time days*status(0);  /* Define survival time and censoring status */
   strata gibleed;       /* Stratify by bleeding status (0 = No, 1 = Yes) */
   title "Kaplan-Meier Survival Curves with 95% CI for Bleeding Status";
run;

/* Log-Rank Test for Bleeding Status */
proc lifetest data=Pbc3;
   time days*status(0);  /* Define survival time and censoring status */
   strata gibleed;       /* Compare survival between bleeding vs. no bleeding */
   title "Log-Rank Test for Bleeding Status";
run;


/* Step 11: Weibull Regression Model */

proc lifereg data=Pbc3;
    class sex tment stage status;  /* Define categorical variables */
    model days*status(0) = age sex tment stage crea alb bili alkph asptr weight / dist=weibull;
    output out=weibull_output Cresidual=csr sresidual=sr;  /* Output residuals */
run;

/* Plot the Standardized Residuals from Weibull Regression */
proc gplot data=weibull_output;
    title "Standardized Residuals from Weibull Regression";
    symbol v=circle c=blue i=none h=1 w=1;
    plot sr*days;  /* Plot standardized residuals against time (days) */
run;

/* Plot the Cox-Snell Residuals from Weibull Regression */
proc gplot data=weibull_output;
    title "Cox-Snell Residuals from Weibull Regression";
    symbol v=circle c=red i=none h=1 w=1;
    plot csr*days;  /* Plot Cox-Snell residuals against time (days) */
run;

/* Step 12: Assessing the Distribution of Cox-Snell Residuals */
proc lifetest data=weibull_output plots=(ls) notable graphics;
    time csr*status(0);  /* Use Cox-Snell residuals as survival time */
run;


 

 
/* Step 1: Fit Weibull regression model and compute residuals */
proc lifereg data=Pbc3;
    class sex tment stage status;  /* Define categorical variables */
    model days*status(0) = age sex tment stage crea alb bili alkph asptr weight / dist=weibull;
    output out=weibull_output Cresidual=csr sresidual=sr;  /* Output Cox-Snell residuals and standardized residuals */
run;

/* Step 3: Plot the standardized residuals */
proc gplot data=weibull_output;
    title "Standardized Residuals from Weibull Regression";
    symbol v=circle c=blue i=none h=1 w=1;
    plot sr*days;  /* Plot standardized residuals against time (days) */
run;

/* Step 4: Plot the Cox-Snell residuals */
proc gplot data=weibull_output;
    title "Cox-Snell Residuals from Weibull Regression";
    symbol v=circle c=red i=none h=1 w=1;
    plot csr*days;  /* Plot Cox-Snell residuals against time (days) */
run;

/* Step 5: Assess whether Cox-Snell residuals follow an Exponential(1) distribution */
proc lifetest data=weibull_output plots=(ls) notable graphics;
    time csr*status(0);  /* Use Cox-Snell residuals as survival time */
run;































































