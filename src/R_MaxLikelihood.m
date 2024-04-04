% Estimation of the daily reproduction number R(t) of an epidemic under a
% Poisson epidemiological model
%
%         Z(t) | Z(1), ..., Z(t-1) = Poisson(R(t) * Zphi(t))
%
% with Zphi(t) = sum_s phi(s) * Z(t-s) where phi is the serial interval
% function.
%
% Maximum Likelihood Estimator writes
%
%         R(t) = Z(t)/Zphi(t),
%
% with by convention R(t) = 0 if Zphi(t) = 0.
%
% Epidemiological model from
% - Cori, A., Ferguson, N. M., Fraser, C., & Cauchemez, S. (2013).
% A new framework and software to estimate time-varying reproduction
% numbers during epidemics. American journal of epidemiology, 178(9),
% 1505-1512.



function R_MLE = R_MaxLikelihood(Z,Zphi)

    % Inputs:  - Z: new infection counts
    %          - Zphi: global infectiousness defined as a weighted sum of past counts
    %
    % Output:  - R: estimated maximum likelihood reproduction number

    R_MLE = Z./Zphi;

    % Handle trivial estimates
    R_MLE(Zphi == 0) = 0;


end