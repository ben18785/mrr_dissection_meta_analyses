A Meta-analysis of Longevity Estimates of Mosquito Vectors of Disease

Ben Lambert^1,2^, Ace North^1^ & H. Charles J. Godfray^1^

^1^ Department of Zoology, University of Oxford, South Parks Road,
Oxford OX1 3PS, United Kingdom

Corresponding author:

Phone: 01865 271176

^2^ Present address: MRC Centre for Outbreak Analysis and Modelling,
Infectious Disease Epidemiology, Imperial College London, London W2 1PG,
UK.

\
=

Abstract
========

Mosquitoes are responsible for more human deaths than any other
organism, yet we still know relatively little about their ecology.
Mosquito lifespan is a key determinant of transmission strength for the
diseases they vector, but the field experiments used to determine this
quantity – mark-release-recapture (MRR) studies and wild-caught
dissection of female mosquitoes – produce estimates with high
uncertainty. In this paper, we use Bayesian hierarchical models to
analyse a previously-published database of 232 MRR experiments and
compile then analyse a database of 131 dissection studies to produce the
first ever species- and genus-level estimates of mosquito lifespan. Due
to the assumptions required to analyse the field data, we term our
estimates lower bounds on lifespan (LBL). Notably, for the major African
malaria vector *Anopheles gambiae s.l.*, we estimate LBLs of 4.4 days
(median estimate; 25%-75% CI: 3.8-5.1 days for unfed female mosquitoes
from the MRR analysis) and 7.6 days (median estimate; 25%-75% CI:
5.2-11.0 days from the dissection analysis); and an LBL of 4.2 days
([]{#_Hlk530592196 .anchor}median estimate; 25%-75% CI: 3.6-4.8 days,
only present in the MRR database) for the predominantly East-African
vector *A. funestus s.l*. We estimate LBLs of 6.2 days (median estimate;
25%-75% CI: 4.5-8.5 days from the MRR analysis) and 4.3 days (median
estimate; 25%-75% CI: 3.5-5.1 days from the dissection analysis) for
*Aedes aegypti*; and 11.6 days (median estimate; 25%-75% CI: 10.0-13.7
days from the MRR analysis) for *Ae.* *albopictus* – the predominant
vectors of dengue fever, chikungunya and Zika. Our estimates indicate
that there is significant variation in lifespan across species, with
most variation explained by differences between genera. In
correspondence with laboratory studies, we estimate that female
mosquitoes outlive males by 0.9 days on average (median estimate;
25%-75% CI: 0.3-1.6 days). We fit models incorporating mosquito
senescence to the data, which allows us to assess evidence for
age-dependent mortality in mosquitoes across different species. We
determine that 8 of 33 species included in the MRR database indicated
evidence for senescence, versus only 2 of 25 species from the dissection
database. Our analysis applies a common framework to the analysis of
databases of MRR and dissection-based experiments, allowing us to
produce robust estimates of lower bounds on lifespan. It also enables us
to critically appraise each field method, which highlights a need for
alternative field methods for measuring this important mosquito
characteristic.

Keywords
========

mosquitoes, mortality, senescence, vector-borne disease, Bayesian,
hierarchical model.

Introduction
============

Some of the most important infectious diseases afflicting humans are
transmitted by mosquitoes, including pathogens such as the causative
agent of malaria that have been associated with humans throughout our
evolutionary history as well recently emergent infections such as the
Zika virus. Most mosquito species have a “gonotrophic cycle” involving
successive episodes of vertebrate blood feeding, egg maturation and
oviposition. In order for a mosquito to transmit a pathogen it must feed
on an infectious person and live long enough to complete at least one
gonotrophic cycle and feed on an uninfected and susceptible individual.
Adult lifespan is thus a critical determinant of the ability of a
mosquito population to allow the persistence of an indirectly
transmitted infection. Lifespan can of course be straightforwardly
assessed in the laboratory, but it is generally accepted that
measurements under relatively benign laboratory conditions are likely to
have limited relevance in the field, and much effort has been directed
at estimating this parameter in the vector’s natural environment. Most
work has focused on assessing average daily mortality rates, and the
simplest assumption is that these do not vary with mosquito age – in
this case longevity is simply the reciprocal of mortality. Testing this
assumption and discovering whether mosquitoes senesce or show other
types of age-dependent mortality has also been studied in the field.

There are two main strategies to estimate mosquito mortality rates and
longevity. The first is through mark-release-recapture (MRR)
experiments, a technique that is widely applied to estimate these
parameters in many types of animal. As applied to mosquitoes, insects
are caught in the field or reared in the laboratory and then marked,
typically with fluorescent dust. The mosquitoes are then released into
the field and then efforts are made to recapture them, for example using
human baits or light traps, usually over an extended period of time.
Mortality rates can be statistically estimated from the numbers of
recaptures given certain assumptions. The main challenges with MRR is
ensuring the marking technique does not affect recapture probability,
and distinguishing mortality from mosquitoes dispersing out of range of
being recaptured. Also, releasing insects that can transmit disease
(especially if this increases ambient population levels) raises
important ethical issues.

The second technique is specific to mosquitoes and makes use of their
gonotrophic cycle. The mosquito ovary is made up of ovarioles, each of
which typically produces one egg every gonotrophic cycle. After the egg
passes into the oviduct the distended ovariole does not completely
recover its previous form but a discrete dilation remains which can be
detected by dissecting the female reproductive organs. Data on the
fraction of females that have oviposited provides some information about
mortality rates. However, a skilled dissector can distinguish the number
of dilations from multiple gonotrophic cycles so providing much richer
data on longevity. The challenges of this method include the amount of
time and expertise it takes to collect the data, establishing the
relationships between physiological and chronological time (though the
distribution of the number of gonotrophic cycles wild-caught mosquitoes
have gone through is of direct epidemiological relevance) and the fact
that it only applies to females.

An issue with both methods is that they require logistically difficult
and expensive field campaigns. There is thus value in conducting a
meta-analysis of existing data to explore consistency across studies,
identify correlates of lifespan and to learn lessons for further
studies. Here we analyse data from 232 MRR and 131 dissection studies
using a common statistical methodology. For MRR we make use of a very
valuable database of 394 mosquito studies assembled by Guerra et al.
(2014) while the dissection studies we extracted from the literature
ourselves. We concentrated on the three major genera of mosquito
vectors, *Anopheles, Aedes* (in its traditional sense) and *Culex*,
which constitute the majority of the data.

Methods
=======

In recent years many important vectors of disease have been shown to be
complexes of very closely related species, biotypes or forms that cannot
be distinguished morphologically (for example the morphospecies
*Anopheles gambiae sensu lato* is now separated into the widespread
*gambiae, coluzzii, arabiensis* and a number of more local species). As
the majority of studies analysed here took place before molecular
techniques allowed these taxa to be separated we work here chiefly with
morphospecies.

Mark-release-recapture
----------------------

Data from MRR experiments in the Guerra et al. (2014) database were
examined and those with fewer than six recaptures and species with only
a single MRR study were excluded for the hierarchical analysis. Of the
232 data sets, 177 involved only females, 35 males, and 18 both sex
releases. For 102 data sets the age of the released mosquitoes was known
(the average age of released mosquitoes was 4.0 days) while in the other
cases it was unknown or unrecorded; in these cases we assumed the
mosquitoes were newly emerged at the time of release and return to this
assumption later.

We analysed all MRR experiments within the same statistical framework
(for full details see the Supplementary Online Material). In the
simplest case *N~R~* mosquitoes are released at day zero and the
probability that they remain in the recapture area until day *t* is
*S*(*t*) when they are recaptured with probability *ψ*. We model the
number of mosquitoes recaptured on day *t* using a negative binomial
sampling model with mean
$\left( N_{R} - Y\left( t - 1 \right) \right)S\left( t \right)\psi$,
where $Y\left( t - 1 \right)$ is cumulative captures before day *t*, and
shape parameter *κ*. The negative binomial has been used previously in
MRR analyses (ref) because of its ability to represent temporal
over-dispersion in recaptures most likely caused by variable weather. A
slight modification was required for studies with multiple releases (see
SOM).

The simplest model for *S*(*t*) assumes there is a constant probability
(*λ*) that a mosquito dies or leaves the recapture area so that the
numbers remaining after time *t* are given by the exponential
distribution, $exp\left\lbrack - \lambda t \right\rbrack$. We utilised
this form extensively but in testing for senescence used five other
models where *λ*(*t*) varies with time so that

  ------------------------------------------------------------------------------- -----
  $$S\left( t \right) = e^{- \int_{0}^{t}{\lambda\left( t \right)\text{dτ}}}.$$   (1)
  ------------------------------------------------------------------------------- -----

Details of the five models (Gompertz, Weibull, Gompertz-Makeham,
Logistic and Logistic-Makeham), which vary in their ability to detect
different forms of age-dependent mortality, are given in the SOM. Using
multiple different types of models increased our chances of detecting
senescence though, as discussed below, also increases the likelihood of
false positives.

Parameters were estimated using Bayesian techniques with relative
uninformative priors for *κ* and the parameters of *λ*(*t*), but
assuming a prior for $\psi$ indicating a low recapture probability
(bounded in part by knowledge of the maximum daily recapture rates). We
used a Bayesian hierarchical model to estimate distributions of lifespan
at the species and the genus levels, and across the complete data set.
This procedure assumes that there is a distribution of lifespan
parameters for each species from which those governing individual MRR
time series are sampled, and similarly a distribution at the genus level
from which those for individual species are derived (rather akin to
random effects in classical statistics). Within this framework we can
also allow the parameters for individual time series to be influenced by
co-variates such as differences in experimental methodology. As in the
estimation of the parameters of the individual experiments, relative
uninformative priors were set for the parameters of the hierarchical
models except for $\psi$ where again a distribution representing low
recapture probabilities was assumed. Posterior distributions were
derived using Markov Chain Monte Carlo (MCMC) methods with convergence
assessed using the $\hat{R}$ statistic (Gelmin & Rubin 1992). The
predictive power of the model was assessed using *K*-fold cross
validation which tests the ability of the model fitted to part of the
data to predict the rest using multiple different partitions. Further
details of the prior specification, fitting and validation through
posterior predictive checks are given in the SOM.

Two studies of *Anopheles balabacensis* reported capture rates
increasing with time, presumably reflecting a violation of our
assumption of constant recapture probabilities. We omitted this species
from the analysis.

The Guerra et al. database included the latitude and longitude of each
study, along

with the date when the study began. We used this information to find
estimates of the

air temperature for each study using the European Centre for Medium
Range Weather Forecasts' ERA Interim Daily historical database. For each
study we calculated the mean monthly temperature across a spatial area
of (latitude +- 1 degree, longitude +- 1 degree), for the month at which
each study was carried out. The records for this database begin in 1979,
which pre-dates the study date for 65 of our 232 MRR time-series. For
these time-series, we chose to estimate the air temperature by an
average of the corresponding monthly temperatures over the years
1979-89.

Dissection
----------

Studies using dissection to estimate mosquito longevity were located in
literature databases using relevant keyword, citation and author
searches, and by checking previous studies cited by the papers located
(see SOM). The list of studies located with associated metadata is
available as a Supplementary Online File.

Most dissection studies recorded the distribution of the number of
gonotrophic cycles in mosquito samples collected over a specific period
of time. Overall, we found 568 time series in 72 published articles.
Because seldom were sufficient ancillary data available data to analyse
trends over time, we aggregated the data from each time series. We
further omitted time series with fewer than 100 mosquitoes and for
species with only one data set leaving 131 studies of mosquitoes in the
genera *Anopheles, Aedes*, *Culex* and *Mansonia*.

To compare lifespan estimates from dissection and MRR studies we need to
convert physiological age (the number of gonotrophic cycles) into
chronological age. Using a literature search and a review by Silver
(2007) we found 79 estimates in 42 published articles. Most estimates
were obtained by dissecting females recaptured in MRR studies or by
observations in the laboratory, the latter tending to give longer
durations. Studies differed greatly in how (if at all) they represented
uncertainty in their estimate of the duration of the gonotrophic cycle.
Where confidence limits were given we treated these as the relevant
quantiles of a normal distribution, where a range was stated (e.g. “4-6
days”) we interpreted the bounds as the 2.5% and 97.5% quantiles of a
normal distribution, and where a single figure was quoted we assumed
this was the mean this distribution. Using the quantiles of the normal
distribution, we estimated its mean and standard deviation by regression
(see SOM). Initially we calculated distributions of gonotrophic cycle
lengths at the species and then genus levels, but because of the paucity
of data for many species and the lack of significant differences we
aggregated the data into a single distribution. We converted
physiological age to chronological age by sampling from this
distribution to obtain a particular gonotrophic cycle length for each
mosquito (we also explored sampling from this distribution to obtain the
duration of *each* gonotrophic cycle which increased the uncertainty in
lifespan estimate but did not affect any of the conclusions).

We modelled the number of mosquitoes found by dissection to be of age
*a* using the negative binomial distribution with mean
$\text{ΨS}\left( a \right)$ and shape parameter *κ*, where *Ψ* is the
product of the recruitment rate of adult mosquitoes, which we assume is
constant over time, and the probability of being captured for
dissection, and *S*(*a*) is the probability of surviving until age *a*.
We used the number of females that have yet to lay eggs (nulliparous) to
estimate the recruitment rate as described further in the SOM. Initial
examination revealed that in some data sets the number of nulliparous
females was anomalously low, something that has been noticed before
(Gillies & Wilkes 1965). As some studies have suggested that the first
gonotrophic cycle tends to be longer than the subsequent ones, this is
probably due to differences in capture probability. In data sets where
the fraction of nulliparous females was less than 90% the uniparous
(completed on gonotrophic cycle) we excluded the nulliparous
observation.

Data was analysed using a Bayesian framework as with the MRR data with
minor differences in the specification of the priors (see SOM). Some
published studies do not distinguish the number of gonotrophic studies
beyond a threshold (the more ovariole dilations there are the harder it
is to count them) which is akin to censoring the data. Because this
censoring involves only a relatively small number of mosquitoes in each
time series (median = 2%), and because of the technical difficulties of
allowing for censoring in our Bayesian estimation procedures, we assume
the females died at the threshold age.

Results
=======

Lifespan estimates from MRR
---------------------------

MRR estimates the length of time a mosquito remains alive and is still
in the area available for recapture. It thus provides a lower bound to
lifespan which we shall refer to as LBL. In 187 of the 230 MRR time
series the estimated LBL was less than 10 days (Fig. 1). The smallest
estimate was &lt;1 day for the Asian malaria vector *Anopheles subpictus
s.l.* which is unfeasibly short and almost certainly reflects dispersal
out of the recapture zone or a violation of the assumptions of our
analyses. The longest estimate was 18.3 days for the temperate species
*Aedes simpsoni s.l.* which is a vector of yellow fever in Africa. There
are multiple data sets for the most important vector species such as
*Anopheles gambiae, Aedes aegypti* and *albopictus* and *Culex tarsalis*
all of which show considerable variation. For example, there are 54
estimates of LBL for *A. aegypti* which range from 2.2 days to 38.3 days
with a mean of 8.3 days and coefficient of variation of 0.7. There are
significant differences in LBL amongst species (ANOVA on median LBL
controlling for sex and pre-release feeding: *F*~37,194~ = 2.5, *p*
&lt;0.01; the non-parametric Kruskal Wallace: ${\chi^{2}}_{38}$,
*p*&lt;0.01).

The estimated median LBL for female mosquitoes not pre-fed with blood or
sugar before release for *Culex, Anopheles* and *Aedes* were 2.5, 5.0
and 6.9 days respectively with an overall estimate of 4.6 days (Fig. 2).
Differences between genera were significant (ANOVA on median LBL
controlling for sex and pre-release feeding: *F*~2,229~ = 12.4,
*p*&lt;0.01; Kruskal Wallace: ${\chi^{2}}_{2} = 30.8,p < 0.01$).
*K*-fold cross validation suggests that after the effect of genus is
accounted for the incorporation of a species term provides little
predictive power (Fig. S\#; in part explained by the latter model
over-fitting the data where there are few time series per species).

We reasoned that if dispersal out of the recapture area was reducing the
LBL below the true lifespan then there should be a positive correlation
between the spatial extent of the recapture zone and LBL. We found no
such pattern (Fig. S\#), although there was a positive correlation
between LBL and trap density (Fig. S\#).

The MRR experiments included a mixture of male-only and female-only
releases, and releases of both sexes. We estimated average male and
female LBL at the genus level (Fig. 3; there were too few studies to
make comparisons at the species level). There was a consistent trend for
females to live longer than males for each of the genera, with the
difference largest for *Aedes* (2.5 days; fraction of pairwise posterior
samples of females versus males where difference was less than zero,
p&lt;0.01), followed by *Anopheles* (2.0 days; p=0.17) and *Culex* (0.3
days; p=0.34). Overall, female mosquitoes were estimated to live 0.9
days longer than males (p=0.10).

The MRR experiments included information on whether mosquitoes were
pre-fed with sugar (41 time series), blood (71), both (4) or
alternatively unfed (116). We estimate that female mosquitoes that were
fed on sugar pre-release lived on average for 0.7 days longer than those
that were not fed (Fig. S\#; a pattern that was consistent across the
genera). There were insufficient males that were either fed or unfed
with sugar prior to release to make a meaningful comparison. Females
that were blood-fed prior to release on average lived 1.5 days longer
than those who were not fed for *Aedes* but this trend was reversed for
*Anopheles* meaning that there was little difference overall.

To access whether temperature is associated with LBL we used weather
records to calculate average temperatures at the MRR sites (see
Methods). Using both linear and quadratic temperature terms in
regressions, we found no significant relationship between study-site
temperature and LBL (overall or within genus) for the 238 datasets we
analysed (Fig. 4). This result held if, instead of pooling results from
all time series, we considered the four species with the most data (*Ae.
aegypti*, *Cx. tarsalis*, *A. gambiae s.l.* and *A. culicifacies s.l.*)
individually (Fig. S\#).

Number of gonotrophic cycles estimates from dissection
------------------------------------------------------

Dissection allows the number of completed gonotrophic cycles to be
counted and from this the mean number of cycles before death was
estimated. Overall, the mean number of cycles completed in a lifetime
was 1.2 (Fig. 5; posterior median) and across the 131 studies, 95% of
the individual time series estimates were less than 3 (Fig. 6). The
estimated greatest number of cycles was for *Anopheles sergentii* (2.5
cycles; posterior median) which is adapted to desert conditions (it is
known as the “oasis vector” of malaria) and may have evolved greater
longevity. The important African malaria vector *An. gambiae s.l.* was
estimated to be the second longest living (1.9 cycles; posterior
median). The smallest estimated mean number of gonotrophic cycles was
for *Anopheles* *bellator* (0.5 cycles; posterior median) which
transmits malaria in Brazil’s Atlantic Forest. There were significant
differences in estimated lifetime gonotrophic cycles amongst species
(ANOVA: *F*~24,106~ =2.2, *p* &lt;0.01; the non-parametric Kruskal
Wallace: ${\chi^{2}}_{24}$, *p*&lt;0.01).

The estimated lifetime gonotrophic cycles for the different genera were
*Anopheles,* 1.4; *Mansonia*, 1.1; *Culex,* 1.0; and *Aedes* 0.8 (Fig.
5) and the differences between the genera were significant (ANOVA:
*F*~3,127~ =3.4, *p* =0.02; the non-parametric Kruskal Wallace:
$\chi_{3}^{2} = 21.7$, *p*&lt;0.01).

Comparison of longevity estimates from two methods
--------------------------------------------------

To compare the two methods, we converted numbers of gonotrophic cycle
(physiological age) into lifespan (chronological age) as described in
the Methods (Fig. S\#). For ten species, we had enough data from both
species to make a comparison, and there was a positive correlation (not
statistically significant; Pearson correlation ρ=0.42, n=10, p=0.23)
between the two measures (Fig. 7), and in only one case – for *A.
darlingi* - (Table S\#) there was a significant difference in the
time-series level LBLs.

Evidence for age-dependent mortality
------------------------------------

The survival model upon which the above analyses are based is the
single-parameter exponential model which assumes an age-invariant
mortality hazard. We also fitted five multi-parameter models that allow,
in different ways, mortality to vary with age. We did this to maximise
our chance of detecting age-varying mortality (though aware of the risks
of false positives with multiple estimations).

In Fig. 8 we compare the performance of the six models for describing
lifespan in MRR studies of 33 species using K-fold cross-validation. We
categorised the evidence for age-dependent mortality in each species
according to the performance of the five age-dependent models versus the
exponential: ‘+’ indicated that all age-dependent models outperformed
the exponential; ‘?’ indicated that the exponential outperformed one or
more age-dependent models; and ‘-’ indicated that the exponential
performed at least as well as all other models. Overall, we estimated
that there were 8 ‘+’ species, where age-dependent mortality fit the
data better; 11 ‘?’ species where the evidence was mixed; and 14 species
where constant mortality models performed at least as well. The species
where age-dependent mortality best fit the data included the major
vector of dengue fever, Zika and chikungunya, *Ae. Aegypti*. These
studies also tended to include multiple release MRR studies which, on
average, were conducted over a longer period of time than the others,
which may be why we failed to detect age-dependence in the latter (Fig
S\#).

In Fig. 9, we compare the performance of the six models for describing
lifespan in dissection studies of 25 species using K-fold
cross-validation, and categorise the evidence in the same way as for the
MRR analysis. By our metric, we determined that there were only two
species with evidence for age-dependent mortality (the major African
malaria vector *A. gambiae s.l.* and *A. minimus*, a malaria vector in
Asia).

Overall, we conclude that there is mixed evidence for age-dependent
mortality from studies of mosquitoes in the field. It is possible that
some of the sampled mosquito species did not live long enough in the
wild to experience physiological decline. A Spearman’s rank correlation
test indicated that there was a correlation between the ranked estimated
LBLs of the species and the ranked mean predictive accuracy of
age-dependent models for the MRR analysis (ρ=0.19, p=0.01), however was
not significant for the dissection analysis (ρ=0.07, p=0.43). Similarly,
a recent study determined that the degree of senescence varies according
to season for semi-wild populations of *Ae. aegypti* (Hugo et al.,
2014), and it is possible that by pooling data from different
geographies and seasons that we failed to detect age-dependent mortality
in some cases.

Estimates of the fraction mosquitoes capable of transmitting disease
--------------------------------------------------------------------

We can use the posterior parameter estimates from our Bayesian analysis
to estimate the fraction of mosquitoes that live beyond a certain age.
In order to transmit a disease, a mosquito must live longer than the
length of the intrinsic incubation period (the time taken for a pathogen
ingested in one blood meal to be ready to be transmitted during a future
feeding event). This is a lower bound as it does not include the waiting
time to find a host after feeding or egg maturation. In Fig. 10 we plot
the fraction of the mosquito population that pass this threshold using
estimates from both MRR and dissection studies for vector species (see
SOM for references used to identify species as vectors) and their most
significant diseases.

For malaria, estimates of the minimum fraction of the population that
can transmit the disease vary from &lt;0.1% for *A. subpictus* (from the
MRR analysis, as noted above likely to be due to the LBL substantially
underestimating lifespan) to 52% for the drought-adapted and long-lived
*A sergentii.* The proportions surviving long enough to become
infectious for the two major African malaria vectors were: for *A.
gambiae s.l.*: 10% (MRR) and 27% (dissection); and for *A. funestus
s.l.*: 9% (MRR). Using the individual time series estimates, however,
there was inconclusive evidence for a difference in EIP between the
species (MRR ANOVA: *F*~14,58~ = 1.8, *p*=0.07; MRR Kruskal-Wallis:
$\chi_{14}^{2} = 30.2$, *p* &lt;0.01; dissection ANOVA: *F*~11,64~ =
1.8, *p*&lt;0.01; dissection Kruskal-Wallis: $\chi_{11}^{2} = 38.9$,
p&lt;0.01).

*Ae. aegypti* and *Ae. albopictus* are the main vectors of dengue,
chikungunya and Zika viruses. Because of their short intrinsic
incubation periods a greater fraction of mosquito potentially live long
enough to transmit diseases (Fig. 10), rising to a maximum of 84% for
*Ae. albopictus* transmitting chikungunya.

Discussion
==========

In this study, we applied a Bayesian hierarchical framework to the
analysis of a database of mark-release-recapture experiments and a
database of mosquito dissection studies to estimate mosquito lifespan.
By applying a single framework, this allows us to effectively synthesise
information from the disparate experiments which, individually, estimate
lifespan with considerable uncertainty. Due to the assumptions required
to analyse the field data, our estimates represent lower bounds on
lifespan (LBL). Across both meta-analyses, the estimated LBLs were
mostly less than 10 days, hinting that only a small proportion of
mosquitoes may live long enough to transmit disease. We determined that
LBL varies across species and genera, although most variance is
explained by genus. The MRR analysis includes experiments conducted on
each sex individually, and we estimate that, on average, males live
shorter lives than females. Pre-release feeding with sugar also
lengthens lifespan across all three genera, although this effect is less
marked than the sex differences. In contrast to a number of lab-based
experiments, temperature was not determined to significantly impact
lifespan. By fitting a range of survival models to the data in both
meta-analyses, we could assess evidence for age-dependent mortality.
Overall, we conclude that the evidence is mixed: in the MRR experiments,
in 8 of 33 species we found evidence for mosquito senescence, whereas in
only 2 of 25 species included in the dissection analysis were better fit
by a model incorporating an increasing risk of mortality with age.

MRR experiments are known to produce downwardly-biased estimates of
lifespan. Lab experiments have demonstrated that marking can negatively
impact survival (Verhurst et al., 2013; Dickens and Brant, 2014),
resulting in artificially depressed survival. MRR studies typically
cannot differentiate between a mosquito dying and dispersal from the
study area meaning that lifespan will be underestimated. In this study,
we found a positive correlation between lifespan estimates and the
density of traps, indicating that better trapping coverage likely raises
estimates towards their real value. We conducted an *in silico* Monte
Carlo study to determine how accurately we could estimate mosquito
lifespan given study parameters in an ideal MRR experiment, where the
assumptions of no emigration and harmless marking are fully satisfied
(see SOM for full details). This work indicated that for many of the
experiments, the short study lengths or typical numbers of mosquitoes
released, results in considerable uncertainty in lifespan estimates (Fig
S\#). This indicates that statistical power can be substantially
increased by pooling data across experiments as we did using a Bayesian
hierarchical model.

The key assumptions of dissection based methods to determine
chronological age are: (i) physiological age can be accurately
determined by dissection of female specimens (unlike MRR, this method
can only be applied to one sex), (ii) the relationship between
physiological and chronological age is known, (iii) the population being
sampled is in equilibrium (recruitment matches mortality) and (iv)
individual mosquitoes can be randomly sampled from the population. The
reliability and accuracy of dissection has been questioned. The
objections include the impracticality of dissecting more than a small
proportion of ovarioles (Hoc and Wilkes, 1995), particularly in African
vector species (Gillies and Wilkes, 1965), the related issue of locating
ovarioles whose count of dilations represents true physiological age
(Fox and Brust, 1994), and the variation in numbers of ovariolar
dilations for mosquitoes of the same, known, physiological age (Kay,
1979; Russell, 1986; Hugo et al., 2008). Indeed there is considerable
uncertainty concerning the fundamental question of how dilations in
ovarioles form in the first place. Whilst the \`Old School' of thought
(a term coined by Fox and Brust, 1994) headed by Polovodana (Polovodova,
1949) and Detinova (Detinova, 1962) considers dilations to result from
normal oogenesis, a \`New School' headed by Lange and Hoc (Lange and
Hoc, 1981) has challenged this assertion. The New School believe that
only abortive oogenesis results in follicular dilations because normal
oogenesis destroys the sack-like structures (Fox and Brust, 1994). This
means that Polovodana's method requires dissecting large numbers of
ovarioles to uncover those with the most dilations, where abortive
oogenesis has occurred in each gonotrophic cycle. They deem these
ovarioles \`diagnostic' since only in these cases the number of
dilations equals the number of gonotrophic cycles that have occurred. As
a mosquito ages, the number of diagnostic ovarioles diminishes, since
the random occurrence of normal oogenesis in a particular ovariole means
its dilation count does not equal the number of gonotrophic cycles
undertaken. This increased difficulty of finding diagnostic ovarioles as
a mosquito ages would elevate the chance of age \`hypodiagnosis' for
older specimens (Fox and Brust, 1994), and likely biases lifespan
estimates downwards. The difficulty of locating diagnostic ovarioles has
been investigated using lab populations of *Culex* and *Aedes*
mosquitoes by Hugo et al., 2008, who conclude that only a small
percentage of ovarioles are diagnostic. The exchange rate between
physiological age and chronological age is the duration of gonotrophic
cycles. Two methods are commonly used to estimate the duration of
gonotrophic cycles: MRR studies (see, for example, Gillies and Wilkes,
1965), where marked mosquitoes are recaptured and dissected to determine
the number of gonotrophic cycles occurring since release; and
laboratory-based observations of colonies of (typically) wild-caught
females, or their progeny (see, for example, Afrane et al., 2005).
Whilst it is unclear how each method could bias estimated gonotrophic
cycle duration, in our analysis, laboratory-based studies indicated a
longer gonotrophic cycle (Fig. S\#). The distributions we used to
convert physiological age into calendar age were calculated by pooling
data across both approaches, to incorporate uncertainty from both
experimental procedures. It is possible, however, that this aggregate
approach may induce biases in estimates and an approach more entrenched
in experimental knowledge would fare better. If a population of
mosquitoes is shrinking, this leads to a relative under-abundance of
young mosquitoes, and a flattening of the survival curve, resulting in
over-estimates of lifespan. For stable populations, periods when
shrinking occurs must result in equal changes in the population size
compared to those when it expands. If mosquito collections occur with
equal frequency in each of these two modes, then aggregating the data
across all sampling times and estimating a single model, as we do here,
should yield an approximately unbiased estimate of lifespan. The
additional uncertainty of a fluctuating population size, however, could
lead us to understate the uncertainty in estimates. Field entomologists
have challenged the assumption of random sampling the mosquito
population, although there are conflicting opinions as to whether this
results in a relative paucity (Gillies and Wilkes, 1965) or abundance
(Clements and Paterson, 1981) of nulliparous individuals. In our
database, there are cases where there was an obvious deficit of
nulliparous individuals, which has previously been ascribed to the
differing distribution of resting females between indoor and outdoor
traps (Detinova, 1962; Gillies and Wilkes, 1965). We chose to not
include those counts of nulliparous individuals in our analysis where
their number was less than 90% of the uniparous. Whilst we see no
obvious differences in lifespan according to collection method (data not
shown) or location, it is possible that the assumption of random
sampling is violated, although the directionality of the bias induced by
this is unclear. Overall, the assumptions underpinning estimates from
dissection studies indicate that our estimates represent lower bounds on
lifespan. The alternative dissection-based approach of Detinova
(Detinova, 1962), based on dichotomous categorisation of female mosquito
specimens as \`parous’ or \`unparous’ relies on fewer assumptions, and
is widely used. Further work examining parity rates in field specimens
may be fruitful although, in principle, it offers less information on
the age structure of a population than Polovodova’s approach.

By applying a common method to analysing all studies in our databases,
it is possible that we may have missed patterns of mortality that would
have been evident from using a more bespoke approach. As our *in silico*
analysis of MRR experiments indicates, however, the overdispersed data
from single experiments results in high measurement error (Fig. S\#). By
applying different methods to each study, this could lead us to falsely
detect patterns when none are present, and we prefer a pooled approach.

The different nature of the assumptions of each of the two methods means
they offer complimentary information on mosquito survival. We also note
that Polovodova’s dissection-based studies require specialised expertise
which will often be unavailable, whereas MRR methods can more readily be
used. Furthermore, most if not all dissection methods that have been
used previously are only applicable to female mosquitoes, whereas MRR
can be applied to either sex and can additionally be used to determine
other ecological parameters (for example, population size and
dispersal). Although dissection data gives detailed of age-structure, we
thus foresee a continued reliance on MRR experiments in field
entomological experiments. Efforts to use both approaches concurrently
will be particularly useful and will allow quantification of the biases
induced by the assumptions of each. Similarly, MRR experiments releasing
large numbers of marked mosquitoes and recording
spatiotemporally-disaggregated captures of wild and re-caught marked
mosquitoes will be especially useful to estimate lifespan and dispersal.

To compare estimates of lifespan derived from MRR with those from
dissection-based methods, we display the estimates of lifespan from
those ten species occurring in both databases in a single plot (Fig. 7).
In is reassuring that there is correlation between estimates from both
approaches, although the small sample size likely hindered our ability
to determine statistical significance. In both cases, we estimate that
*A. sergentii* was amongst the longest lived of the anopheline species
with a posterior median LBL of 8.7 days (median estimate; 25%-75% CI:
5.9-13.8 days) from the MRR analysis and 10.0 days (median estimate;
25%-50% CI: 7.6-14.0 days) from the analysis of dissection studies. This
species is a major vector of malaria in the Sahara (Sinka et al., 2010),
where to act as a disease vector it must persevere through these hard
conditions. It is reasonable to hypothesise that this species should
live longer than those in environments where the potential for
blood-feeding and oviposition is greater. The species with the greatest
discrepancy in the estimates was the major African malaria vector *A.
gambiae s.l.*, where we estimated LBLs of 4.4 days (median estimate;
25%-75% CI: 3.8-5.1 days for unfed female) from the MRR analysis and 7.6
days (median estimate; 25%-75% CI: 5.2-11.0) from the dissection
analysis. Across genera, the greatest discrepancy in estimates was for
*Culex* (a posterior median of 2.5 days from the MRR versus 4.4 days
from the dissection analysis) with the smallest discrepancy for
*Anopheles* (4.8 versus 5.0 days). By contrast for *Aedes* the estimates
from the MRR studies (6.9 days) exceed those of dissection-based studies
(3.5 days). Across all studies we estimate from the MRR analysis that
mean mosquito lifespan is 4.6 days versus 5.1 days from the
dissection-based studies. Some of the differences in these group-level
estimates between the two approaches is likely due to environmental and
genetic differences between mosquitoes in the experiments that were
analysed in each meta-analysis. However, we believe that part of the
discrepancy can be explained by the methodological differences in
approaches. We speculate that differences in dispersal rate can explain
some of the discrepancy. Both *Anopheles* and *Culex* mosquitoes are
generally thought to fly farther during their lifetimes than *Aedes*
\[Charles, do you have a reference here?\], meaning that the estimates
from MRR-based approaches will be most downwardly-biased for these
genera. This is supported by our results since the dissection-based
estimates (themselves not reliant on assumptions about dispersal) exceed
the MRR estimates for *Anopheles* and *Culex* mosquitoes, but not for
*Aedes*.

It is widely believed mosquitoes live artificially long under the benign
conditions of the laboratory. We find it informative to consider
estimates of lifespan derived from observations of such populations as
they constitute an upper bound on the lifespan of wild populations.
Also, since the numbers of mosquitoes involved in large cage experiments
often numbers in the thousands, these estimates have lower uncertainty
than those from field experiments although are typically conducted on
highly inbred mosquito strains. Styer et al. (2007), using colonies of
45,054 female and 55,997 male *Ae. aegypti*, determined that females
lived nearly twice as long as males; the median lifespan was estimated
as 31.69 days +- 0.06 days for females and 16.39 +- 0.03 days for males.
A similar study by Dawes et al., 2009 with a lab colony of over 1000
female *A. stephensi* found similar estimates for median lifespan (31-42
days). These estimates are many multiples of the average estimates that
result from our analysis of field data which, as discussed, represent
lower bound estimates. Without an unbiased method to measure mosquito
lifespan, however, it is difficult to quantify and explain the gap that
exists between field and laboratory lifespans. The development of
additional methods to estimate mosquito age, such as \`Near-Infrared
Spectroscopy’ (Mayagaya et al., 2009; Sikulu et al., 2011; Lambert et
al., 2018), if they are proven to work in the field, may be of
considerable worth here.

We conducted a power analysis of MRR experiments (Fig. S\#) to determine
whether typical experimental characteristics could detect senescence.
Here we calculated the power of a maximum likelihood estimator of the
\`senescence parameter' β of the Gompertz survival function (see Table
S\#) for case study populations with three different levels of
senescence (Fig. S\#). This analysis indicated that power to detect
senescence strongly depends on study length (Fig. S\#) but is
insensitive to release size (Fig. S\#). Clements and Patterson (1981)
conducted a meta-analysis of MRR and dissection-based field experiments
and found evidence of an increasing risk of mortality hazard with age
that is similar in magnitude to that of the \`mild' case considered
above. For this case, detecting senescence with a power of 80% requires
a study length of at least 18 days. Since the median study length for
experiments included in our analysis was 10 days (Table S\#) this could
partly explain our failure to detect senescence at the species level. A
number of experiments have found evidence of age-dependence in
laboratory populations (Styer et al., 2007; Dawes et al., 2009).
However, the artificially benign environment of the laboratory means
mosquitoes live considerably longer than in the wild, where they may die
because of exogenous factors, before the effects of physiological
decline have had time to manifest. Field experiments have also found
evidence for age-dependent mortality. Harrington et al. (2008) conducted
a field experiment where mosquitoes reared under laboratory conditions
were marked and released at different ages. Analysis of the resultant
MRR time-series indicated that mosquito mortality increases with age at
release. It is possible, however, that this field experiment suffers
from the same biases as laboratory-based approaches, because the
released mosquitoes were often of ages considerably higher (up to 20
days) than typical estimates of wild mosquito lifespan.

As ethical concerns of contributing to disease burden are more often
considered, it is now less common for MRR experiments to release female
mosquitoes versus males than historically (Fig. S\#). Our analysis
indicates that females outlive male mosquitoes by approximately 0.9 days
(Fig. \#), meaning that differences between the sexes may exist for
other ecological parameters determinable by MRR. This suggests that
continued field entomological work on contained releases of mosquitoes
in semi-field sites or large microcosms may be a valuable source of
information on female mosquitoes.

Our estimates of LBL indicate that mosquitoes that were sugar-fed prior
to release lived on average 0.7 days longer than those that were unfed
(, suggesting the potential value of this underappreciated aspect of the
mosquito ecology to the insects. It may also partly explain the recent
successes in the use of Attractive Toxic Sugar Baits (ATSB) as a vector
control intervention (Schlein and Muller references). More research is
needed, however, to identify the sugar-feeding frequency and food
sources for wild populations.

There is considerable evidence from the laboratory that temperature
modulates mosquito ecology and behaviour (Matt Thomas references). The
locations and times of year over which the MRR studies were conducted
encompassed a large range of average air temperatures, from
approximately 10 ^o^C to 35 ^o^C and, within this, we determined no
relationship between lifespan and temperature across all time series
(Fig. 4) or, for any of the species with the most data (Fig. S\#). It is
possible that by considering a raw average of air temperature across the
month, this ignored other temperature-dependent factors of importance
for longevity, such as temperature range. It is also possible that by
ignoring the effects of rainfall (the historical data on rainfall is
less likely to be reliable for a given location), that this masked a
more complex interaction between longevity and temperature. The observed
laboratory relationship between lifespan and temperature, however, may
not be as robust in the field if mosquitoes adjust their behaviours
(such as, by seeking shade) in reaction to changes in temperature. More
work exploring the relationship between mosquito ecology and temperature
in semi-field experiments may be useful in probing these interactions
further.

In this work, we have used modern statistical methods to synthesise
precious field data conducted by entomologists past and present, to
produce lower bound estimates of mosquito lifespan. The importance of
vector mortality for disease transmission has long been recognised,
however, since even before 1957, when George Macdonald formulated the
now famous Ross-Macdonald equation of R~0~ for malaria. Indeed, the
recent declines in malaria prevalence in Sub-Saharan Africa were likely
due to the upscaling of interventions (insecticide-treated bednets and
indoor residual spraying) that aim to reduce mosquito lifespan (Bhatt et
al., 2015). Worryingly, resistance to pyrethroids, the only class of
insecticide used in current insecticide-treated bednets and likely the
only product to come to market in the near future, has been determined
to be widespread and increasing in intensity across Sub-Saharan Africa
(WHO reference 2018). This alarming trend highlights the need for
continued MRR and dissection-based studies to monitor the effectiveness
of bednets and determine whether more expensive alternatives, such as
nets incorporating piperonyl butoxide be deployed. It also emphasises
the need for investment in new tools for real time monitoring of
mosquito populations. In recent years, there has been considerable
funding allocated to molecular and genomic research into mosquitoes that
strengthens existing interventions and suggest novel control strategies.
Without commensurate funding allocated to applied vector ecology, our
lack of knowledge in this area threatens our opportunity to capitalise
on molecular advances and potentially hinders our ability to control of
mosquito-borne disease.

Acknowledgements
================

The authors would like to thank the following for useful conversations
throughout the course of this work: Austin Burt, Thomas Churcher, Steve
Lindsay, Mike Bonsall and Ellie Sherrard-Smith.

References
==========

ADDIN EN.REFLISTX

Figure Legends
==============

Figure 1:
