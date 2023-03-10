Libname respir "C:\Users\ferruccilu\Box\PAD_respirometry\";
Options Linesize=132;
* 
* Imports the data on Respirometry
* from C:\Users\ferruccilu\Box\PAD_respirometry\simple.csv
*;
PROC IMPORT OUT= RESPIR.SIMPLE
            DATAFILE= "C:\Users\ferruccilu\Box\PAD_respirometry\simple.csv"
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2;
RUN;
option ls=120; run;
proc print data=respir.simple; run;

Data respir.analytical;
        set respir.simple;
		keep Code PAD_non_PAD Complex_II Complex_IV
	     logComplex2 logComplex4 missComplex2 missComplex4  
         n_missing anymissing;
		 n_missing=numero_missing;
		 if (Complex_II eq .) then missComplex2=1;
		 else                      missComplex2=1;
		 if (Complex_IV eq .) then missComplex4=1;
		 else                      missComplex4=1;
		 logComplex2=log2;
		 logComplex4=log4;
        if (code eq .) then code=1;
run;
*
* Now look at the missing data (Supplementary figure 8.A)
*;
Proc freq data=respir.analytical;
        tables PAD_non_PAD*missComplex2 PAD_non_PAD*missComplex4 
        PAD_non_PAD*n_missing PAD_non_PAD*anymissing/ chisq;
run;
*
* Exports the analytical database into analytical./CSV
*;
PROC EXPORT DATA= RESPIR.ANALYTICAL
            OUTFILE= "C:\Users\ferruccilu\Box\PAD_respirometry\analytical.csv"
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
*
* Export data for the boxplot in supplemental Figure 8 
* (right) that will be plotted in Origins Pro
*;
data boxplot1temp (keep=PAD log2 log4)
     boxplot2temp (keep=PAD log2 log4);
        set respir.simple;
        if (pad eq 0) then output boxplot1temp;
        else if (PAD eq 1) then output boxplot2temp;
run;
data boxplot1;
        set boxplot1temp;
        Nopad=pad;
        NoPADlog2=log2;
        NoPADlog4=log4;
        drop pad nopad log2 log4;
run;
data boxplot2;
        set boxplot2temp;
        Nopad=pad;
        PADlog2=log2;
        PADlog4=log4;
        drop pad nopad log2 log4;
run;
data boxplot;
        merge boxplot1 boxplot2;
run;
*
* Exports the data 
* for the bpxplots in Figure S8 (right)
*;
PROC EXPORT DATA= boxplot
            OUTFILE= "C:\Users\ferruccilu\Box\PAD_respirometry\boxplot.csv"
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
*
* Now imports the data on covariates
* from "Complex2_4_Scores_Reduced.csv"
*;
PROC IMPORT OUT= WORK.respiratory_data
            DATAFILE= "C:\Users\ferruccilu\Box\PAD_respirometry\Complex2_4_Scores_Reduced.csv"
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2;
RUN;
proc contents data=work.respiratory_data position; run;
Data proteomics;
	set work.respiratory_data;
	keep code ABI Six_minute_walk__meters__x Group GroupReg Age Sex
         BlocksWalkedLog2 BlocksWalked_x ABI_Score MitoScore
         Cplx2Score Cplx4Score;
run;
proc sort data=respir.simple; by code; run;
proc sort data=work.proteomics; by code; run;
data mixed;
        merge respir.simple (in=a) proteomics (in=b);
        by code;
        AA=A;
        BB=B;
        if (AA and BB);
run;
*
* Now recreates this mixed file, 
* for the 39 participants who have some data on 
* respiratory complexes;
*;
data mixed_all;
        merge respir.simple (in=a) proteomics (in=b);
        by code;
        AA=A;
        BB=B;
*
* Creates a classification for the table;
*;
        if (BB eq 1) then             class1="All participants    ";
        if (AA eq 1 and BB eq 1) then class2="Specimen available  ";
        if (Complex_II ne .) then     class3="With Complex II data";
        if (Complex_IV ne .) then     class4="With Complex IV data";
*
* selects case with data on respirometry
*;
		if (Class3 ne "" or class4 ne "");
		sixmin_walk=Six_minute_walk__meters__x;
		drop Six_minute_walk__meters__x;
run;
proc contents data=mixed_all position; run;



*********************************************************************************
* Spearman correlations and Robust regression for Figure 4 for the PAD paper
*;
proc sort data=mixed_all;
        by pad_non_pad;
run;
Proc corr data=mixed_all spearman fisher;
        var log2 log4;
        by pad_non_pad;
run;
proc robustreg data=mixed_all;
        model log2=log4;
        by pad_non_pad;
run;
*
* Figures on the right
*;
options nocenter;
Proc corr data=mixed_all spearman fisher;
        var log2 mitoscore;
        by pad_non_pad;
run;
proc robustreg data=mixed_all;
        model log2=mitoscore;
        by pad_non_pad;
run;
Proc corr data=mixed_all spearman fisher;
        var log4 mitoscore;
        by pad_non_pad;
run;
proc robustreg data=mixed_all;
        model log4=mitoscore;
        by pad_non_pad;
run;

Proc corr data=mixed_all spearman fisher;
        var log2 cplx2score;
        by pad_non_pad;
run;
proc robustreg data=mixed_all;
        model log2=cplx2score;
        by pad_non_pad;
run;

Proc corr data=mixed_all spearman fisher;
        var log4 cplx4score;
        by pad_non_pad;
run;
proc robustreg data=mixed_all;
        model log4=cplx4score;
        by pad_non_pad;
run;

*********************************************************************************
*********************************************************************************
* Robust Regression for supplementary figure 9
*;
proc sort data=mixed_all; by pad_non_pad; run;
Proc corr data=mixed_all;
        var sixmin_walk log2;
        by pad_non_pad;
run;
proc robustreg data=mixed_all;
        model sixmin_walk=log2;
        by pad_non_pad;
run;
Proc corr data=mixed_all spearman fisher;
        var sixmin_walk log4;
        by pad_non_pad;
run;
proc robustreg data=mixed_all;
        model sixmin_walk=log4;
        by pad_non_pad;
run;
Proc corr data=mixed_all spearman fisher;
        var log2 abi;
        by pad_non_pad;
run;
proc robustreg data=mixed_all;
        model log2=abi;
        by pad_non_pad;
run;
Proc corr data=mixed_all spearman fisher;
        var log4 abi;
        by pad_non_pad;
run;
proc robustreg data=mixed_all;
        model log4=abi;
        by pad_non_pad;
run;
*********************************************************************************;
*********************************************************************************;
*  Comparison of complexII and complexIV between PAD and non PAD.
*  Using the Wilcoxon test
*;  
proc NPAR1WAY data=mixed_all wilcoxon;
	title "Nonparametric test to compare complex II activity between PAD and controls";
	class pad_non_pad;
	var log2;
	exact wilcoxon;
run;
proc NPAR1WAY data=mixed_all wilcoxon;
	title "Nonparametric test to compare complex IV activity between PAD and controls";
	class pad_non_pad;
	var log4;
	exact wilcoxon;
run;
*
* Attached the following CVS files
* SIMPLE.CVS
* COMPLEX2_4_SCORES_REDUCED.CVS
* BOXPLOT.CVS
* Thia code was created in SAS 9.4 TS and saved as PAD_muscle_paper.sas;
* reviewed on Feb 18th 2023
*;

 

