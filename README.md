# SOAP-opera

An attempt to combine SOAP stellar activity models with ABC

Should I find a backronym for *opera*?

## Basic idea

SOAP is a bit slow to create a simulated data set of flux, RVs and CCF indicators from an active star.
But we would like to fit observed data with a physically motivated activity model. Maybe ABC can help, if we manage to find good descriptive statistics of the observables.

The idea is to build an ABC model with `DNest4` where the data is generated by a call to SOAP. The parameters being inferred would be the active region properties (maybe the number of active regions as well?). 
This can be used to fit observed or simulated data of active stars.

Why `DNest4`? Because it provides an efficient implementation of diffusive nested sampling, an algorithm which is still efficient in sampling multimodal distributions with big degeneracies.


Further down the line, the model can be extended by considering two types of objects whose number needs to be inferred: active regions and planets. `N_ar` and `N_p` would be parameters in the MCMC, included in a run of `DNest4`'s `RJObject`. 
This huge model can then produce transit data, RVs and CCF indicators, to be compared with observations.
I'm guessing it will be slow... But let's try it!

## What is needed

- [x] Create basic `DNest4` or `RJObject` template
- [ ] Create and set priors for active region parameters
- [ ] Implement a call to SOAP
- [ ] *enhancement:* could also call the `macula` code to create photometry
- [ ] Implement ABC "likelihood"
- [ ] Find descriptive statistics
- [ ] *enhancement:* is the CCF itself a better descriptive statistic then RVs, FWHM, BIS, etc?
- [ ] Create simulated data for tests
- [ ] Extend `RJObject` to allow for two classes of objects, active regions and planets
- [ ] Sample........
