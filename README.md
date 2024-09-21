# Estimation of Covid19 reproduction number via nonsmooth convex optimization

**Contributors:** P. Abry (1), J. Du (1), C.-G. Lucas (1), B. Pascal (2), and N. Pustelnik (1).

(1) CNRS, ENS de Lyon Laboratoire de physique F-69007 Lyon, France.  
(2) Nantes Université, École Centrale Nantes, CNRS, LS2N, UMR 6004, F-44000 Nantes, France.

Fundings: ANR-23-CE48-0009 [OptiMoCSI](https://optimocsi.cnrs.fr/) and CNRS 80PRIME-2021 «CoMoDécartes» project.

## Regularized estimates of Covid19 pandemic reproduction number

This project contains the `Matlab` codes associated to the journal papers:

> Abry, P., Pustelnik, N., Roux, S., Jensen, P., Flandrin, P., Gribonval, R., Lucas, C.-G., Guichard, É., Borgnat, P., & Garnier, N. (2020). Spatial and temporal regularization to estimate COVID-19 reproduction number R(t): Promoting piecewise smoothness via convex optimization. *PlosOne*, 15(8), e0237901.
> [url](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0237901)
>
> Pascal, B., Abry, P., Pustelnik, N., Roux, S., Gribonval, R., & Flandrin, P. (2022). Nonsmooth convex optimization to estimate the Covid-19 reproduction number space-time evolution with robustness against low quality data. *IEEE Transactions on Signal Processing*, 70, 2859–2868.
> [arxiv:2109.09595](https://arxiv.org/abs/2109.09595)


`Python` versions of these codes developed from the present toolbox by J. Du in collaboration with P. Abry and B. Pascal is available in [Covid-R-Estim](https://github.com/juliana-du/Covid-R-estim).

Three demonstration scripts are provided:
- [`demo_R_World`](https://github.com/bpascal-fr/Covid-Estim-R/blob/master/demo_R_World.m) for the daily estimation of the reproduction number in the 200+ countries monitored by [Johns Hopkins University](https://coronavirus.jhu.edu/) independently;
- [`demo_R_France`](https://github.com/bpascal-fr/Covid-Estim-R/blob/master/demo_R_France.m) for the daily estimation of the reproduction number in the 104 French departments for which daily infection counts are reported by [Santé Publique France](https://www.data.gouv.fr/fr/datasets/donnees-de-laboratoires-pour-le-depistage-a-compter-du-18-05-2022-si-dep/), including possible multivariate estimation leveraging *spatial* regularization;
- [`demo_R_Advanced`](https://github.com/bpascal-fr/Covid-Estim-R/blob/master/demo_R_Advanced.m) for the daily estimation of the reproduction number in the 200+ countries monitored by [Johns Hopkins University](https://coronavirus.jhu.edu/) independently under a more advanced model in which the reproduction number is fixed for two days, proposed in:
> Abry, P., Chevallier, J., Fort, G., & Pascal, B. (2023, December).  Pandemic intensity estimation from Stochastic Approximation-based algorithms. *In 2023 IEEE 9th International Workshop on Computational Advances in Multi-Sensor Adaptive Processing (CAMSAP)* (pp. 356-360). IEEE. [hal-04174245](https://hal.science/hal-04174245);
- [`demo_Daily_to_Weekly`](https://github.com/bpascal-fr/Covid-Estim-R/blob/master/demo_Daily_to_Weekly.m) for aggregating daily infection counts into weekly infection counts and computing the associated global infectiouness with an adapted *weekly* discretized serial interval function.

Animated maps presenting the obtained estimates have been developed by the consortium and are available at [covidatlas.eu](https://www.covidatlas.eu/World/).

## Project description

In an epidemic outbreak context, such as for example during the Covid19 pandemic, it is of utmost importance for Health Authorities to achieve close daily surveillance of the pathogen transmission.
The intensity of an epidemic can be quantified at day $t$ through its *effective reproduction number* $\mathsf{R}_t$, defined as the *average number of secondary infections caused by a standard contagious individual* [1]. The major advantage of the effective reproduction number is that it has a straightforward interpretation:

- if $\mathsf{R}_t > 1$, then the epidemic is spreading exponentially fast, and an epidemic wave is to be expected in the absence of countermeasures; 
  
- while if $\mathsf{R}_t < 1$, the epidemic is regressing, indicating for example the effectiveness of some social distancing rules.

The estimation of the reproduction number on a daily basis crucially relies on collected new infection counts.
In an emergency situation, these counts are of medium to low quality, due to errors and imprecisions occurring during the measurement and reporting processes, such as irrelevant
or missing counts, pseudo-seasonalities, week-end, etc. The low quality of collected data severely incurs accurate estimation of $\mathsf{R}_t$ from state-of-the-art epidemiological tools [1].
To provide robust estimates of the daily reproduction number, capable to manage low quality Covid19 data, [2,3] elaborated on the propagation model proposed in [1] and designed regularizing convex functionals balancing the fidelity to the epidemiological mechanisms and spatiotemporal regularity constraints, leading to four estimators:

- Univariate ($\mathsf{U}$): enforcing temporal piecewise linearity of the estimated $\widehat{\mathsf{R}}_t$;

- Univariate Corrected ($\mathsf{U}$-$\mathsf{C}$): enforcing temporal piecewise linearity of the estimated $\widehat{\mathsf{R}}_t$ and including the estimation of a sparse corrective term accounting for misreported infection counts;
  
- Multivariate ($\mathsf{M}$): enforcing temporal piecewise linearity and spatial piecewise constancy of the estimated $\widehat{\mathsf{R}}_{c,t}$ in county $c$ at day $t$, based on a prior connectivity structure between counties;

 - Multivariate Corrected ($\mathsf{M}$-$\mathsf{C}$): enforcing temporal piecewise linearity and spatial piecewise constancy of the estimated $\widehat{\mathsf{R}}_{c,t}$ based on a prior connectivity structure between counties, and including the estimation of a sparse corrective term accounting for misreported infection counts.

Fast and memory efficient implementation of the nonsmooth convex optimization problems defining the aforementioned estimators permits daily estimation of the reproduction number in hundreds, possibly connected, territories simultaneously.

[1] Cori, A., Ferguson, N. M., Fraser, C., & Cauchemez, S. (2013). A new framework and software to estimate time-varying reproduction numbers during epidemics. *American Journal of Epidemiology*, 178(9), 1505-1512.

[2] Abry, P., Pustelnik, N., Roux, S., Jensen, P., Flandrin, P., Gribonval, R., Lucas, C.-G., Guichard, É., Borgnat, P., & Garnier, N. (2020). Spatial and temporal regularization to estimate COVID-19 reproduction number R(t): Promoting piecewise smoothness via convex optimization. *PlosOne*, 15(8), e0237901.

[3] Pascal, B., Abry, P., Pustelnik, N., Roux, S., Gribonval, R., & Flandrin, P. (2022). Nonsmooth convex optimization to estimate the Covid-19 reproduction number space-time evolution with robustness against low quality data. *IEEE Transactions on Signal Processing*, 70, 2859–2868.

## Installation and dependencies

To download and install the toolbox:  
> - open a Terminal in your working folder,
> - execute `git clone https://github.com/bpascal-fr/Covid-Estim-R.git`.
>

The Matlab toolbox [stein-piecewise-filtering](https://github.com/bpascal-fr/stein-piecewise-filtering) should be downloaded and placed in the folder containing the toolbox:
> - `cd Covid-Estim-R`,
> - `git clone https://github.com/bpascal-fr/stein-piecewise-filtering.git`.