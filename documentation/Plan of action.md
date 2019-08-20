### Parity

* [ ] Write up method and data description for parity analysis
  * [x] Include Stan code for model fitting
  * [x] Include priors
* [x] For parity series, compare fit of species, group, geography (Africa, Asia, Americas), and overall models. Redo using K fold CV
* [x] Redo CV for continent vs complex-level because - at the moment - I am not using all data points (because of being constrained by species). Think this is fine -- leave as is.
* [x] Redo parity analysis where I ensure I don't get duplicates between species and complex estimates (e.g. *A. moucheti*) and where I remove species or groups with fewer than 5 obs
* [ ] Rewrite SOM section on data criteria for inclusion (and edits made to data) to incorporate
  * [x] 5 or more data points, for species and complex estimates and cross validation (also continental CV)
  * [x] all data for overall and continental comparison
  * [x] Include information on the CV method used
* [ ] Make estimation notebook which does lifespan estimation (and creates folds), then produces PPCs
  * [ ] species (check numbers of species match with lifespan diagrams)
  * [ ] complex
  * [ ] all 
* [x] Tidy up visualisations of lifespan for species and complex
* [x] Make visualisations for continental and overall
* [x] Posterior predictive checks for parity analysis. Redo with eCDFs and include in SOM
* [ ] Individual lifespan estimates
* [x] Write up results in 'parity_results.tex' file.
* [ ] Gambiae:
  * [ ] Estimate gambiae s.l. lifespan by country (for countries with at least a minimum number of studies) and compare with grouped estimates via k-fold CV
  * [ ] Estimate lifespan by month in a particular location/country
* [ ] Weather vs lifespan investigations:
  * [x] Night temperature range (00:00)
  * [ ] Day-night temperature range
  * [ ] Monthly std
  * [x] Mean temperature
  * [ ] Monthly min
  * [ ] Monthly max
  * [ ] Precipitation
  * [ ] Wet vs dry season

### Physiological age

- [ ] Temperature investigation:
  - [ ] Extract temperature measures for studies
  - [ ] Cross-validation with mean, daily range, standard deviation in mean, min, max
- [ ] Group results by complex
- [ ] Group according to Africa, Asia, the Americas (and potentially estimate at this level)

### MRR

- [ ] Redo temperature extraction using 0.1 degrees and also extract mins and maxes
  - [ ] Cross-validation with mean, daily range, standard deviation in mean, min, max
- [ ] Group results by complex
- [ ] Group according to Africa, Asia, the Americas (and potentially estimate at this level)

### Overall

* Gonotrophic cycle estimates by genus (perhaps even revisit this approach as not convinced it's ideal)
* Dissection power analysis SOM (I think I just wanted to look at senescence here?)
* Look at Github issues
* Table of references for included studies (or does this exist in the Excel files?)
* Quantify difference in lifespan from dissection vs parity
  * Can do a power analysis too

## Nice to have

* Make overall estimates where we combine all three data sources
* Do model fitting for MRR:
  * Male / female vs no differences (have I done this before?)
  * Pre-feeding