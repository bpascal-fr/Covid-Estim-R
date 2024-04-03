% Proximity operator of the weighted squared l2 norm, defined as
%
%      argmin_Q 1/2 * ||R - Q||^2 + gamma * ||Z - Q * Zphi||^2.
%
% Implementation N. Pustelnik, B. Pascal, C.-G. Lucas and P. Abry
% April 2020
%
% Updated and augmented by P. Abry and B. Pascal
% March 2024

function  prox    = prox_L2w(R,Z,Zphi,gamma)

    % Inputs:  - R: reproduction number stored in a C x T matrix, C number of territories, T number of days
    %          - Z: new infection counts stored in a C x T matrix, C number of territories, T number of days
    %          - Zphi: infectiousness stored in a C x T matrix, C number of territories, T number of days
    %          - gamma: parameter of the proximity operator (descent step)  (gamma > 0)
    %
    % Outputs: - prox: proximal operator of parameter gamma of the squared l2 norm between Z and R * Zphi

    prox = (R + gamma.*Zphi.*Z)./(1 + gamma.*Zphi.^2);

end