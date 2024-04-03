% Proximity operator of the Kullback-Leibler divergence between
% Z and R * Zphi at R defined as
%
%      argmin_Q 1/2 * ||R - Q||^2 + gamma * DKL(Z | Q Zphi)
%
% and has closed form expression
%
%      prox_gamma*DKL (R) = (R - gamma.*Zphi + sqrt((R-gamma.*Zphi).^2 + 4*gamma.*Z))/2.
%
% Implementation N. Pustelnik, B. Pascal, C.-G. Lucas and P. Abry
% April 2020
%
% Updated and augmented by P. Abry and B. Pascal
% March 2024

function prox = prox_DKLw(R,Z,Zphi,gamma)

    % Inputs:  - R: reproduction number stored in a C x T matrix, C number of territories, T number of days
    %          - Z: new infection counts stored in a C x T matrix, C number of territories, T number of days
    %          - Zphi: infectiousness stored in a C x T matrix, C number of territories, T number of days
    %          - gamma: parameter of the proximity operator (descent step)  (gamma > 0)
    %
    % Outputs: - prox: proximal operator of parameter gamma of the Kullback-Leibler divergence between Z and R * Zphi

    prox  = (R - gamma.*Zphi + sqrt(abs(R-gamma.*Zphi).^2 + 4*gamma.*Z))/2;
    
    % Forces output to zero if Z(t) = Zphi(t) = 0
    index = find((Zphi==0)&(Z==0)) ;

    if ~isempty(index) 

        prox(index) = zeros(size(index));

    end

end



