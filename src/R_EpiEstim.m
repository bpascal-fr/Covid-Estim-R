% Estimation of the daily reproduction number R(t) of an epidemic under a
% Poisson epidemiological model
%
%         Z(t) | Z(1), ..., Z(t-1) = Poisson(R(t) * Zphi(t))
%
% with Zphi(t) = sum_s phi(s) * Z(t-s) where phi is the serial interval
% function.
%
% Maximum A Posteriori with a Gamma prior on R(t) and assuming that R(t) is
% constant across the past tau days:
%    
%         R_EpiEstim(t) = ( a + <Z>_tau(t)) / (1/b + <Zphi>_tau(t))
%
% from
% - Cori, A., Ferguson, N. M., Fraser, C., & Cauchemez, S. (2013).
% A new framework and software to estimate time-varying reproduction
% numbers during epidemics. American journal of epidemiology, 178(9),
% 1505-1512.



function R_Gamma = R_EpiEstim(Z,Zphi,tau)

    % Resize input
    [d1,d2]     = size(Z);
    if min(d1,d2) == 1
        Z       = reshape(Z,1,max(d1,d2));
        Zphi    = reshape(Zphi,1,max(d1,d2));
    end


    % default parameters of the prior Gamma(a,b)
    a       = 1; 
    b       = 5; 

    % define the window extracting past tau days
    win     = zeros(1,2*tau-1); win(end-tau+1:end) = 1;

    % compute the cumulative <Z>_tau and <Zphi>_tau
    sZ      = convn(Z,win,'same');
    sZphi   = convn(Zphi,win,'same');

    % explicit expression of the Maximum A Posteriori
    R_Gamma = (a+sZ)./(1/b + sZphi);

    % Resize output
    R_Gamma = reshape(R_Gamma,d1,d2);

end