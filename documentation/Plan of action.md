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
* [x] Gambiae:
  * [x] Estimate gambiae s.l. lifespan by country (for countries with at least a minimum number of studies) and compare with grouped estimates via k-fold CV
  * [ ] Estimate lifespan by month in a particular location/country
* [ ] Weather vs lifespan investigations:
  * [x] Night temperature range (00:00)
  * [ ] Day-night temperature range
  * [ ] Monthly std
  * [x] Mean temperature
  * [ ] Monthly min
  * [ ] Monthly max
  * [ ] Precipitation
  * [ ] Wet vs dry season from downloaded data
  * [ ] Wet vs dry season from database
* [ ] Insecticide use (recorded in database)
* [x] Gonotrophic cycle estimation using citations
  * [x] Recheck for duplicates
  * [x] Review methods for each gonotrophic cycle estimate
  * [x] Need to determine what measures these were and what method was used (29 references in total); no need, I think
  * [x] Combine with my literature search for physiological age

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

## Gonotrophic cycle estimates

- [ ] Replace old gonotrophic cycle figure in SOM with new ones with updated data
- [ ] ANOVA on medians:
  - [ ] Genus differences for both 1st and 2nd durations
  - [ ] No differences between 1st and 2nd durations when controlling for Genus
- [ ] Update SOM text:
  - [ ] on method used to analyse
  - [ ] literature search (explain inclusion criteria and how the Massey database was included); describe how - if no mention was made of 1st or subsequent cycle, we assumed it was a measure of both -- meaning the observations often overlapped
- [ ] Include Mathematica plots of the estimated densities for 1st and subsequent gonotrophic cycle durations

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

## Needed for write up

* Lifespan estimates (use posterior medians!):
  * gambiae s.l.:
    * [x] MRR: 4.4 days
    * [x] Detinova: 8.7 days
    * [x] Polovodova: 6.4 days
  * funestus s.l.:
    * MRR: 4.2 days
    * Detinova: 13 days
  * albopictus:
    * MRR: 11.6
    * Polovodova: NA
  * aegypti: 
    * MRR: 6.2
    * Polo: 4.7
  * Male vs female: shouldn't have changed as based only on MRRs
* Add notes on:
  * [x] Mansonia: assumed overall gonotrophic estimates parameters
  * [x] Add notes on conversion from parity to chronological age: assumed a single gonotrophic cycle length
  * [x] Add notes on physiological age conversion
* Replace figures:
  * [x] Dissection lifespans
  * [x] Gonotrophic cycle plots and estimates
* [x] Check MCMC running characteristics (number of chains, iterations):
  * [x] Parity (think this is wrong in SOM)
  * [x] Polovodova
* [x] Add notes on truncation for Polovodova dissection analysis: used 14 age classes if there was no thresholding, so should be a reasonable approximation in those cases; where thresholding occurs, there is no difference
* [x] Create a file that runs Polovodova analysis: no need.
* [x] Add glossary of terms: MRR, Polovodova-type dissection, Detinova-type
* [x] Do k fold with new kappa prior for gambiae sl
* [ ] Make all axis titles capitalised
* [x] Look why moucheti doesn't seem to be being fit with species / complex fits. Think this is a mistake because Moucheti has only a single species. Need to rerun! Has been rerun. No need to rerun k-fold for species vs grouped since Moucheti not in species. Nor for group vs continent as seems I was using 778 data points before.
* [x] Added estimates for Detinova without insecticide
* [x] check that quadrimaculatus isn't being missed off of detinova estimates: it gets removed since it has too few data

# November picking things up

* [x] Add details of insecticide analysis to supplementary
* [x] Check out Polovodova estimates for whether or not we do actually record insecticide use (we say we do in the supplementary methods): we do but going to leave out as it seems unreliable.
* [x] Rerun Polovodova with more iterations: actually looks fine in terms of Rhat with 16 chains
* [x] Find latest Detinova estimates: they are generated in "s_fitting_insecticide" but the Stan fits no longer exist! I must have deleted them. There are no summary lifespans, so I need to regenerate these.
* [x] Rerun Detinova fits with 4000 iterations (some hadn't converged before)
  * [x] Check Rhat for each fit (all below 1.1)
  * [x] If Rhat ok, update SOM with new numbers of iterations
* [x] Regenerate Detinova lifespans: it's in combined_lifespans_no_insecticide
  * [x] gonotrophic cycles
  * [x] calendar
* [x] Detinova: Update numbers of observations in cross-validation of SOM to reflect no insecticide case
* [x] Update main Detinova lifespan figure and move to write up repo
* [x] Update Detinova lifespan comparison figure and include in figure file
* [x] Include Detinova in EIP figure
* [x] Check and probably change gonotrophic cycle length quoted in main text
* [x] Go through SOM gonotrophic cycle duration bit and update numbers of studies included. Looks fine.
* [x] Detinova: sort k fold
  * [x] Rerun k-fold cross validation with more iterations (and update SOM if needs be). Some issue with Africa-species -- sorted
  * [x] Update main text with numbers from this
* [x] Check and include Polovodova senescence estimates
* [x] Update titles for supplementary figures
* [ ] Go through all numbers in text. To indicate that they have been checked this time, I am using a green colour.
* [x] Update text about pairwise figure
* [x] Update text about EIP figure
* [x] Fill in "navigating_analyses" markdown text
  * [x] Also add to it the locations of the three main folders
* [x] Make all axes labels match (e.g. detinova boxplot and polovodova boxplot)
* [ ] Make all axis titles capitalised:
  * [x] Main figures
  * [ ] SOM materials figures
* [x] Check all tables have updated lifespan estimates:
  * [x] Detinova
  * [x] Polovodova
  * [x] MRR
* [x] Weather analysis section update
* [x] Clear up comments
* [ ] Add references to text as comments
* [x] Go through manuscript and add references to relevant SOM figures
* [ ] Update PPCs for parity data
* [x] Go through all numbers in text. To indicate that they have been checked this time, I am using a green colour