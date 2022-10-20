#setwd()
#install.packages("survival")
#install.packages("survminer")
library(survival)
library(survminer)

data = read.table('CGGA.mRNAseq_693_clinical.20200506.txt',header = T, sep = '\t',row.names = 1)
print(data)

Fit = survfit(Surv(OS,Censor..alive.0..dead.1.) ~ Gender , data = data)
print(Fit)

ggsurvplot(Fit,data = data)

ggsurvplot(Fit,data = data, surv.median.line = 'hv',pval = T,risk.table = T)

Fit = survfit(Surv(OS,Censor..alive.0..dead.1.) ~ IDH_mutation_status , data = data)
print(Fit)

Fit = survfit(Surv(OS,Censor..alive.0..dead.1.) ~ IDH_mutation_status + Gender , data = data)
print(Fit)

