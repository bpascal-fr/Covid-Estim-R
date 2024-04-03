% Kullback-Leibler divergence between Z and R * Zphi + O
%
%      DKL(Z | R Zphi + O) = sum dkl(Z(t) | R(t) Zphi(t) + O(t))
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


function DKL = DKLw_outlier(X,Z,Zphi)

    % Inputs:  - X: reproduction number and corrective term stored in a C x 2T matrix, C number of territories, T number of days
    %               R stored in the first T columns
    %               O stored in the last T columns
    %          - Z: new infection counts stored in a C x T matrix, C number of territories, T number of days
    %          - Zphi: infectiousness stored in a C x T matrix, C number of territories, T number of days
    %
    % Outputs: - DKL: Kullback-Leibler divergence between Z and R * Zphi + O

    R = X(:,1:size(Zphi,2));
    O = X(:,size(Zphi,2)+1:end);
    P = Zphi .* R + O;

    j = find(Z>0  & P>0 );
    k = find(Z==0 & P>=0);

    if length(j) + length(k) < numel(P)

        DKL  = Inf;

    else

        DKLj      = P(j) - Z(j) + Z(j).*log(Z(j)./P(j));
        DKLk      = P(k);

        DKL  = sum(DKLj(:))+sum(DKLk(:));

    end


end
