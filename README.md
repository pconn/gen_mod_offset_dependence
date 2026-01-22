Welcome to a repository hosting code for several simulation studies and an example analysis that examines 
the potential for bias in generalized models for count data when there is correlation between the offset (e.g., 
effort/exposure) and the underlying process intensity.  We include a simple Poisson GLM simulation, as well
as a more sophisticated simulation study examining the performance of abundance estimators when density surface models (DSMs; sensu Miller et al. 2013)
are fitted to animal transect data.  In particular, we examine the case where detection and density processes 
are correlated in space.  We also analyze a beluga data set.

This code is associated with the paper, "Dependence between effort offsets and process intensity 
in generalized models for count data: spatial models of animal abundance"
by P.B. Conn and M.C. Ferguson (in revision at Environmetrics).

The Poisson glm simulation study can be replicated vi athe R script, /inst/sim_pois_exposure.R
The DSM simulation study can be replicated via the R script, /inst/sim_study.R
The beluga analysis can be replicated via the R script /inst/xy_tw_DSM_2022.R
