% Estimation of the daily reproduction number R(t) of an epidemic under a
% Poisson epidemiological model
%
%         Z(t) | Z(1), ..., Z(t-1) = Poisson(R(t) * Zphi(t))
%
% with Zphi(t) = sum_s phi(s) * Z(t-s) where phi is the serial interval
% function.
%
% Mean A Posteriori with a Gamma prior on R(t) and assuming that R(t) is
% constant across the past tau days:
%    
%         R_EpiEstim(t) = ( a + <Z>_tau(t)) / (1/b + <Zphi>_tau(t))
%
% accompanied with alpha credibility intervals.
%
% from
% - Cori, A., Ferguson, N. M., Fraser, C., & Cauchemez, S. (2013).
% A new framework and software to estimate time-varying reproduction
% numbers during epidemics. American journal of epidemiology, 178(9),
% 1505-1512.



function [R_Gamma,CI] = R_EpiEstim(Z,Zphi,tau,alpha)

    % Inputs:  - Z: new infection counts
    %          - Zphi: global infectiousness defined as a weighted sum of past counts
    %          - tau: number of days during which R is assumed constant
    %          - alpha: coverage of the credibility interval (default: 0.95)
    %
    % Output:  - R_Gamma: estimated maximum a posteriori reproduction number
    %          - CI_Gamma: alpha credibility interval
    %                    CI.low: lower bound of the interval
    %                    CI.upp: upper bound of the interval

    % Resize input
    [d1,d2]     = size(Z);
    if min(d1,d2) == 1
        Z       = reshape(Z,1,max(d1,d2));
        Zphi    = reshape(Zphi,1,max(d1,d2));
    end
    
    % default credibility interval
    alpha   = 0.95;

    % default parameters of the prior Gamma(a,b)
    a       = 1; 
    b       = 5; 

    % define the window extracting past tau days
    win     = zeros(1,2*tau-1); win(end-tau+1:end) = 1;

    % compute the cumulative <Z>_tau and <Zphi>_tau
    sZ      = convn(Z,win,'same');
    sZphi   = convn(Zphi,win,'same');

    % explicit expression of the Maximum A Posteriori
    shape   = (a+sZ);
    scale   = 1./(1/b + sZphi);
    R_Gamma = scale.*shape;

    % compute the alpha credibility interval
    delta   = (1-alpha)/2;
    CI.low  = gaminv(delta,shape,scale);
    CI.upp  = gaminv(1-delta,shape,scale);

    % Resize output
    R_Gamma = reshape(R_Gamma,d1,d2);
    CI.low  = reshape(CI.low,d1,d2);
    CI.upp  = reshape(CI.upp,d1,d2);

end