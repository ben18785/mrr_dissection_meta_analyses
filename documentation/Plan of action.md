## Essential

### Parity

* [ ] Write up method and data description for parity analysis
  * [x] Include Stan code for model fitting
  * [x] Include priors
* [x] For parity series, compare fit of species, group, geography (Africa, Asia, Americas), and overall models. Redo using K fold CV
* [ ] Redo CV for continent vs complex-level because - at the moment - I am not using all data points (because of being constrained by species)
* [x] Redo parity analysis where I ensure I don't get duplicates between species and complex estimates (e.g. *A. moucheti*) and where I remove species or groups with fewer than 5 obs
* [ ] Rewrite SOM section on data criteria for inclusion (and edits made to data) to incorporate
  * [ ] 5 or more data points, for species and complex estimates and cross validation (also continental CV)
  * [ ] all data for overall and continental comparison
  * [ ] Include information on the CV method used
* [ ] Make estimation notebook which does lifespan estimation (and creates folds), then produces PPCs
  * [ ] species (check numbers of species match with lifespan diagrams)
  * [ ] complex
  * [ ] all 
* [x] Tidy up visualisations of lifespan for species and complex
* [x] Make visualisations for continental and overall
* [x] Posterior predictive checks for parity analysis. Redo with eCDFs and include in SOM
* [ ] Individual lifespan estimates

### Overall

* Gonotrophic cycle estimates by genus (perhaps even revisit this approach as not convinced it's ideal)
* Scrutinise species classifications:
  * [ ] MRR
  * [ ] Dissection
  * [x] Parity
* Dissection power analysis SOM (I think I just wanted to look at senescence here?)
* Look at Github issues
* Table of references for included studies (or does this exist in the Excel files?)
* Quantify difference in lifespan from dissection vs parity
  * Can do a power analysis too

## Nice to have

* Make overall estimates where we combine all three data sources
* Group-level estimates of lifespan? Based on scrutinising the species classifications
* Identify which species are where:
  * MRR
  * Dissection
* Plot individual estimates of lifespan on a map:
  * MRR
  * Dissection
  * Parity
* Do model fitting for MRR:
  * Male / female vs no differences (have I done this before?)
  * Pre-feeding
  * Temperature