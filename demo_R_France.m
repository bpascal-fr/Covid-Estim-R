clear all
close all
clc

addpath(genpath(pwd))

% All estimators process infection count times series of size D x T
% - D: number of French departments 
% - T: number of days in the time period.
%
% Univariate estimators are computed in parallel for each department.

%% LOAD AND DISPLAY NEW INFECTION COUNTS IN FRANCE

% List all French departments, including overseas territories
list                = 1; % 1 display the list, 0 do not show it
AllDepartments      = get_FrenchDepartments(list); 

% List of departments to monitor
User_Departments    = ["44","59","69","75","971","975"]; % Official French numerotation (correspondence provided by get_FrenchDepartments, second column: Number)
% If [] all 104 French departments are selected.
% Departments are stored in the order indicated in User_Departments or in the official INSEE order with Corsica between Finistère (29) and Gard (30) if User_Departments = [].

% Time period
opts_load.LastDay   = '2023-02-05'; % Last day of the selected time period in format 'YYYY-MM-DD'. If not provided set to June 27, 2023.
opts_load.T         = 70;           % Length of the time period. If -1 or not provided entire available time period.

% Load data
[Z, Zphi, Dates, Departments]    = load_SIDEP_France(User_Departments,opts_load);

% Store departments, dates, infection counts and infectiousness
results.Departments = Departments;
results.Dates       = Dates;
results.Z           = Z;
results.Zphi        = Zphi;

% Departments to be displayed
DisplayDepartments    = ["59","69","75","975"]; % Official French numerotation (correspondence provided by get_FrenchDepartments, second column: Number) 

% Display infection counts time series for monitored countries
display_Counts_France(DisplayDepartments,results) % If DisplayDepartments = [] all monitored departments are displayed.

% Prepare storage of estimates
results.Estimates = [];

%% MAXIMUM LIKELIHOOD ESTIMATOR AND BAYESIAN ESTIMATOR

% Following the Poisson epidemiological model of 
% Cori et. al, 2013, Am. J. Epidemiol. the maximum likehood
% estimator is R_MLE(t) = Z(t)/Zphi(t)

R_MLE             = Z./Zphi;

results.Estimates = [results.Estimates, "MLE"];
results.MLE       = R_MLE;


% Maximum A Posteriori with a Gamma prior on R(t) and assuming that
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

custom = 0;
% 0: all parameters of the estimator are set to default values
% 1: some parameters are set manually (more parameters listed in R_Univariate)

if custom == 0
    
    % Perform estimation via functional minimization
    scale              = std(Z,[],2);   % scale of infection counts
    [R_U,obj_U,incr_U] = R_Univariate(Z./scale,Zphi./scale,lambda_T);
    
else

    % Customized settings of the variational estimator
    opts_U.Ri   = R_MLE;  % initialization
    opts_U.prec = 1e-7;   % required precision on increments for convergence
    opts_U.iter = 1e5;    % maximal number of iterations
    
    % Perform estimation via functional minimization
    scale              = std(Z,[],2);   % scale of infection counts
    [R_U,obj_U,incr_U] = R_Univariate(Z./scale,Zphi./scale,lambda_T,opts_U);

end

results.Estimates      = [results.Estimates, "U"];
results.U              = R_U;


%% UNIVARIATE SIMULTANEOUS ESTIMATION OF THE REPRODUCTION NUMBER AND CORRECTION OF MISREPORTED COUNTS

% Regularized log-likelihood estimator enforcing smooth temporal behavior
% and sparsity of the corrective term. From
% Pascal et al., 2022, IEEE Trans. Sig. Process.

% Temporal regularization parameter
lambda_T = 50;

% Sparsity of the corrective term
lambda_O = 0.2;

custom = 0;
% 0: all parameters of the estimator are set to default values
% 1: some parameters are set mcanually (more parameters listed in R_Univariate_Correct)

if custom == 0
    
    % Perform estimation via functional minimization
    scale                          = std(Z,[],2);   % scale of infection counts
    [R_U_C,O_U_C,obj_U_C,incr_U_C] = R_Univariate_Correct(Z./scale,Zphi./scale,lambda_T,lambda_O);
    
else

    % Customized settings of the variational estimator
    opts_U_C.xi   = R_MLE;  % initialization
    opts_U_C.prec = 1e-7;     % required precision on increments for convergence
    opts_U_C.iter = 1e7;  % maximal number of iterations

    % Perform estimation via functional minimization
    scale                          = std(Z,[],2);   % scale of infection counts
    [R_U_C,O_U_C,obj_U_C,incr_U_C] = R_Univariate_Correct(Z./scale,Zphi./scale,lambda_T,lambda_O,opts_U_C);

end

results.Estimates = [results.Estimates, "U-C"];
results.U_C       = R_U_C;


%% MULTIVARIATE ESTIMATION OF THE REPRODUCTION NUMBER 

% Regularized log-likelihood estimator enforcing smooth temporal behavior
% and spatial regularity across connected territories. From 
% Abry et al., 2020, PlosOne
% 

% Temporal regularization parameter
lambda_T = 50;

% Graph discrete gradient operator on connected departments
G = Laplacian_French_Dept(Departments);

% Spatial regularization parameter
lambda_S = 0.005;

custom = 1;
% 0: all parameters of the estimator are set to default values
% 1: some parameters are set manually (more parameters listed in R_Multivariate)

if custom == 0
    
    % Perform estimation via functional minimization
    scale              = std(Z,[],2);   % scale of infection counts
    [R_M,obj_M,incr_M] = R_Multivariate(Z./scale,Zphi./scale,lambda_T,G,lambda_S);
    
else

    % Customized settings of the variational estimator
    opts_M.Ri   = R_MLE;  % initialization
    opts_M.prec = 1e-7;   % required precision on increments for convergence
    opts_M.iter = 1e5;    % maximal number of iterations

    % Perform estimation via functional minimization
    scale              = std(Z,[],2);   % scale of infection counts
    [R_M,obj_M,incr_M] = R_Multivariate(Z./scale,Zphi./scale,lambda_T,G,lambda_S,opts_M);

end

results.Estimates      = [results.Estimates, "M"];
results.M              = R_M;


%% MULTIVARIATE SIMULTANEOUS ESTIMATION OF THE REPRODUCTION NUMBER AND CORRECTION OF MISREPORTED COUNTS

% Regularized log-likelihood estimator enforcing smooth temporal behavior
% and spatial regularity across connected territories and sparsity of the 
% corrective term. From 
% Pascal et al., 2022, IEEE Trans. Sig. Process.
% 

% Temporal regularization parameter
lambda_T = 50;

% Graph discrete gradient operator on connected departments
G = Laplacian_French_Dept(Departments);

% Spatial regularization parameter
lambda_S = 0.005;

% Sparsity of the corrective term
lambda_O = 0.2;

custom = 0;
% 0: all parameters of the estimator are set to default values
% 1: some parameters are set manually (more parameters listed in R_Multivariate_Correct)

if custom == 0
    
    % Perform estimation via functional minimization
    scale                           = std(Z,[],2);   % scale of infection counts
    [R_M_C,O_M_C,obj_M_C,incr_M_C] = R_Multivariate_Correct(Z./scale,Zphi./scale,lambda_T,G,lambda_S, lambda_O);
    
else

    % Customized settings of the variational estimator
    opts_M_C.Ri   = R_MLE;  % initialization
    opts_M_C.prec = 1e-7;   % required precision on increments for convergence
    opts_M_C.iter = 1e5;    % maximal number of iterations
    
    % Perform estimation via functional minimization
    scale                           = std(Z,[],2);   % scale of infection counts
    [R_M_C,O_M_C,obj_M_C,incr_M_C] = R_Multivariate_Correct(Z./scale,Zphi./scale,lambda_T,G,lambda_S,lambda_O,opts_M_C);

end


results.Estimates      = [results.Estimates, "M-C"];
results.M_C            = R_M_C;



%% DISPLAY INFECTION COUNTS AND ESTIMATED REPRODUCTION NUMBER

% Estimators to be compared
DisplayEstimates    = ["MLE","Gamma","U","U-C","M","M-C"];
% - MLE:   maximum likelihood estimate
% - Gamma: maximum a posteriori under Gamma prior from Cori et al. (2013)
% - U:     univariate estimate with piecewise linear temporal behavior from Abry et al. (2020)
% - U-C:   univariate estimate with piecewise linear temporal behavior and estimation of corrective terms from Pascal et al. (2022)
% - M:     multiivariate estimate with piecewise linear temporal behavior and piecewise constant spatial behavior from Abry et al. (2020)
% - M-C:   univariate estimate with piecewise linear temporal behavior, piecewise constant spatial behavior and estimation of corrective terms from Pascal et al. (2022)
% If DisplayEstimates = [] all estimates available will be displayed.

% Display infection counts and estimates per department
display_Estim_France(DisplayDepartments,DisplayEstimates,results) % DisplayDepartments = [] displays all analyzed departments.

% Departments to be compared
CompareDepartments   = ["44","59","69","75","971"]; 

% Compare reproduction number estimates in different departments
compare_Estim_France(CompareDepartments,DisplayEstimates,results)
