% Proximity operator of the weighted squared l2 norm with erroneous count
% correction, defined as
%
%      argmin_(Q,U) 1/2 * ||R - Q||^2 + 1/2 * ||R - Q||^2 + gamma * ||Z - Q * Zphi - U||^2.
%
% % using the closed form expression of the proximity operator of the
% weighted squared l2 norm
% and the fact that (R,O) -> (Zphi * R + O) satisfies a frame property.
%
% Implementation N. Pustelnik, B. Pascal, C.-G. Lucas and P. Abry
% April 2020
%
% Updated and augmented by P. Abry and B. Pascal
% March 2024


function prox = prox_L2w_outlier(X,Z,Zphi,tau)

    % Inputs:  - X: reproduction number and corrective term stored in a C x 2T matrix, C number of territories, T number of days
    %               R stored in the first T columns
    %               O stored in the last T columns
    %          - Z: new infection counts stored in a C x T matrix, C number of territories, T number of days
    %          - Zphi: infectiousness stored in a C x T matrix, C number of territories, T number of days
    %          - tau: parameter of the proximity operator (descent step)  (tau > 0)
    %
    % Outputs: - prox: proximity operator of parameter tau of the squared l2 norm between Z and R * Zphi + O


    R = X(:,1:size(Zphi,2));
    O = X(:,size(Zphi,2)+1:end);
    P = Zphi'.*R' + O';

    % leverage the frame property to express the proximity operator
    prox_L2 = prox_L2w(P,Z',1,tau*(Zphi.^2+1)');
    prox = X - [Zphi.*(P'-prox_L2')./(Zphi.^2+1),...
        (P'-prox_L2')./(Zphi.^2+1)];


    % Forces output to zero if Z(t) = Zphi(t) = 0
    proxR                      = prox(:,1:size(Zphi,2));
    proxO                      = prox(:,size(Zphi,2)+1:end);
    proxR(Z == 0 & Zphi == 0)  = 0;
    proxO(Z == 0 & Zphi == 0)  = 0;
    prox(:,1:size(Zphi,2))     = proxR;
    prox(:,size(Zphi,2)+1:end) = proxO;

end



