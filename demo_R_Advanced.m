clear all
close all
clc

addpath(genpath(pwd))

% All estimators process infection count times series of size C x T
% - C: number of territories monitored
% - T: number of days in the time period.
%
% Univariate estimators are computed in parallel for each territory.


%% LOAD AND DISPLAY NEW INFECTION COUNTS

% List all countries referenced in JHU repository
list             = 1;
get_AllCountries(list)
% list = 0: do not display countries
% list = 1: display all countries names

% List of countries to monitor
User_Countries             = ["France","United Kingdom","Spain"]; % If [] all 201 countries referenced in the JHU repository are selected.
% Countries are stored in the order indicated in User_Countries or by alphabetic order if User_Countries = [].

% Time period
opts_load.LastDay           = '2023-02-05'; % Last day of the time period in format 'YYYY-MM-DD'. If not provided set to March 9, 2023.
opts_load.T                 = 70 + 2;       % + 2 due to two first days of fixed R. Length of the time period. If -1 or not provided entire available time period. 

% Load data
[Z, Zphi, Dates, Countries] = load_JHU_World(User_Countries,opts_load);

% Store countries, dates, infection counts and infectiousness
results.Countries = Countries;
results.Dates     = Dates;
results.Z         = Z;
results.Zphi      = Zphi;

% Countries to be displayed
DisplayCountries    = ["France","United Kingdom"];

% Display infection counts time series for monitored countries
results.FontSize  = 22.5; % Adjust font size in plots
display_Counts_World(DisplayCountries,results) % If DisplayCountries = [] all monitored countries are displayed.

% Prepare storage of estimates
results.Estimates = [];

%% MAXIMUM LIKELIHOOD ESTIMATOR AND BAYESIAN ESTIMATOR

% Following the Poisson epidemiological model of 
% Cori et. al, 2013, Am. J. Epidemiol. the maximum likehood
% estimator is R_MLE(t) = Z(t)/Zphi(t)

R_MLE             = R_MaxLikelihood(Z,Zphi);

results.Estimates = [results.Estimates, "MLE"];
results.MLE       = R_MLE;


% Mean A Posteriori with a Gamma prior on R(t) and assuming that
% R(t) is constant across a window of tau days proposed by
% Cori et. al, 2013, Am. J. Epidemiol.

tau               = 7; %(default choice in Cori et al.) 
R_Gamma           = R_EpiEstim(Z,Zphi,tau);

results.Estimates = [results.Estimates, "Gamma"];
results.Gamma     = R_Gamma;

%% UNIVARIATE ESTIMATION OF THE REPRODUCTION NUMBER 

% Regularized log-likelihood estimator enforcing smooth temporal behavior.
% From Abry et al., 2020, PlosOne
% 
% Temporal regularization parameter
lambda_T = 50;

custom   = 0;
% 0: all parameters of the estimator are set to default values
% 1: some parameters are set manually (more parameters listed in R_Univariate)

if custom == 0
    
    % Perform estimation via functional minimization
    [R_U_12,obj_U_12,incr_U_12] = R_Univariate_R1R2(Z,Zphi,lambda_T);
    
else

    % Customized settings of the variational estimator
    opts_U_12.R1                = R_MLE(:,1);      % initialization
    opts_U_12.R2                = R_MLE(:,2);      % initialization
    opts_U_12.Ri                = R_MLE(:,3:end);  % initialization
    opts_U_12.prec              = 1e-7;            % required precision on increments for convergence
    opts_U_12.iter              = 1e5;             % maximal number of iterations
    
    % Perform estimation via functional minimization
    [R_U_12,obj_U_12,incr_U_12] = R_Univariate_R1R2(Z,Zphi,lambda_T,opts_U_12);

end

results.Estimates      = [results.Estimates, "U-12"];
results.U_12           = R_U_12;

%% UNIVARIATE SIMULTANEOUS ESTIMATION OF THE REPRODUCTION NUMBER AND CORRECTION OF MISREPORTED COUNTS

% Regularized log-likelihood estimator enforcing smooth temporal behavior
% and sparsity of the corrective term. From
% Pascal et al., 2022, IEEE Trans. Sig. Process and
% Abry et al., 2023, IEEE CAMSAP.
%
% Extended version keeping fixed the two first values R(1) and R(2).

% Temporal regularization parameter
lambda_T = 3.5;

% Sparsity of the corrective term
lambda_O = 0.05;

custom   = 0;
% 0: all parameters of the estimator are set to default values
% 1: some parameters are set mcanually (more parameters listed in R_Univariate_Correct)

if custom == 0
    
    % Perform estimation via functional minimization
    [R_U_C_12,O_U_C_12,obj_U_C_12,incr_U_C_12] = R_Univariate_Correct_R1R2(Z,Zphi,lambda_T,lambda_O);
    
else

    % Customized settings of the variational estimator
    opts_U_12.R1   = R_MLE(:,1);      % initialization
    opts_U_12.R2   = R_MLE(:,2);      % initialization
    opts_U_12.Ri   = R_MLE(:,3:end);  % initialization
    opts_U_C_12.prec = 1e-7;          % required precision on increments for convergence
    opts_U_C_12.iter = 1e7;           % maximal number of iterations

    % Perform estimation via functional minimization
    [R_U_C_12,O_U_C_12,obj_U_C_12,incr_U_C_12] = R_Univariate_Correct_R1R2(Z,Zphi,lambda_T,lambda_O,opts_U_C_12);

end

results.Estimates    = [results.Estimates, "U-C-12"];
results.U_C_12       = R_U_C_12;


%% DISPLAY INFECTION COUNTS AND ESTIMATED REPRODUCTION NUMBER

% Estimators to be compared
DisplayEstimates    = ["MLE","Gamma","U-12","U-C-12"];
% - MLE:   maximum likelihood estimate
% - Gamma: maximum a posteriori under Gamma prior from Cori et al. (2013)
% - U:     univariate estimate with piecewise linear temporal behavior from Abry et al. (2020)
% - U-C:   univariate estimate with piecewise linear temporal behavior and estimation of corrective terms from Pascal et al. (2022)
% If DisplayEstimates = [] all estimates available will be displayed.

% Display infection counts and estimates
display_Estim_World_R1R2(DisplayCountries,DisplayEstimates,results) % DisplayCountries = [] displays all analyzed countries.


