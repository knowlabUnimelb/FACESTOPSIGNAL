load wiener
load glm
load dic
model in "C:\Users\littled\Dropbox\Work\FACESTOPSIGNAL\analysis\Bayesian_Group_Analysis\stopSignalModel.txt"
data in jagsdata.R
compile, nchains(1)
parameters in jagsinit2.R
initialize
update 5000
monitor set mu_m, thin(10)
monitor set prec_m, thin(10)
monitor set mu_T0, thin(10)
monitor set prec_T0, thin(10)
monitor set mu_tau, thin(10)
monitor set prec_tau, thin(10)
monitor set m, thin(10)
monitor set T0, thin(10)
monitor set tau, thin(10)
monitor set d, thin(10)
monitor set c, thin(10)
monitor set dprime, thin(10)
monitor set msamp, thin(10)
monitor set T0samp, thin(10)
monitor set tausamp, thin(10)
monitor deviance
update 50000
coda *, stem('CODA2')
