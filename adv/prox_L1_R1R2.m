function prox = prox_L1_R1R2(LR,gamma,R1,R2)

    % Inputs:  - LR: partial Laplacian of R excluding the two first components of R at day 1 and day 2 stored as a C x (T-2) matrix
    %          - gamma: parameter of the proximity operator (descent step)  (gamma > 0)
    %          - R1: first component of R for each territory, stored as a C x 1 columns vector
    %          - R2: second component of R for each territory, stored as a C x 1 column vector
    %
    % Outputs: - prox: of the l1 norm with offset on the two first components
    % 
    %      argmin_Q ||Q - LR||^2 + gamma * |LR(1) - R2/2 + R1/4| + gamma * |LR(2) + R2/4| + gamma * ||LR_3:T||_1
    %
    % Implementation N. PUSTELNIK, CNRS, ENS Lyon
    % June 2019
    %
    % Updated to exclude the two first components
    % B. Pascal, CNRS, LS2N
    % April 2024

    LR(:,1)   = LR(:,1) - R2/2 + R1/4;
    LR(:,2)   = LR(:,2) + R2/4;

    prox      = max(abs(LR)-gamma,0).*sign(LR);

    prox(:,1) = prox(:,1) + R2/2 - R1/4;
    prox(:,2) = prox(2) - R2/4;

end
