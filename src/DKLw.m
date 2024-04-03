% Kullback-Leibler divergence between Z and R * Zphi
%
%      DKL(Z | R Zphi) = sum dkl(Z(t) | R(t) Zphi(t))
%
% with
%   dkl(z | p) = z log(z/p) + p -z  if z>0 & p>0,
%   dkl(z | p) = p                  if z = 0 & p >=0
%   dkl(z | p) = Inf                otherwise.
%
% Implementation N. Pustelnik, CNRS, ENS Lyon
% April 2020
%
% Updated and augmented by P. Abry and B. Pascal
% March 2024


function DKL = DKLw(R,Z,Zphi)

    % Inputs:  - R: reproduction number stored in a C x T matrix, C number of territories, T number of days
    %          - Z: new infection counts stored in a C x T matrix, C number of territories, T number of days
    %          - Zphi: infectiousness stored in a C x T matrix, C number of territories, T number of days
    %
    % Outputs: - DKL: Kullback-Leibler divergence between Z and R * Zphi

    % Find indices for which dkl is finite
    j = find(Z>0  & Zphi.*R>0 );
    k = find(Z==0 & R>=0);

    % Compute the sum of individual dkl 
    if length(j) + length(k) < length(R)

        DKL  = Inf;

    else

        DKLj      = Zphi(j) .* R(j) - Z(j) + Z(j).*log(Z(j)./(Zphi(j).*R(j)));
        DKLk      = Zphi(k) .* R(k);

        DKL  = sum(DKLj(:))+sum(DKLk(:));

    end

end