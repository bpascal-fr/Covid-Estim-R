% Estimation of the daily reproduction number R(t) of an epidemic from new
% infection counts time series Z(t).
%
% Extended version keeping fixed the two first values R(1) and R(2).
%
% from
% - Cori, A., Ferguson, N. M., Fraser, C., & Cauchemez, S. (2013).
% A new framework and software to estimate time-varying reproduction
% numbers during epidemics. American journal of epidemiology, 178(9),
% 1505-1512.
%
% - Abry, P., Pustelnik, N., Roux, S., Jensen, P., Flandrin, P.,
% Gribonval, R., Lucas, C.-G., Guichard, É., Borgnat, P.,
% & Garnier, N. (2020). Spatial and temporal regularization to estimate
% COVID-19 reproduction number R(t): Promoting piecewise smoothness via
% convex optimization. PlosOne, 15(8), e0237901
%
% - Abry, P., Chevallier, J., Fort, G., & Pascal, B. (2023, December). 
% Pandemic intensity estimation from Stochastic Approximation-based algorithms. 
% In 2023 IEEE 9th International Workshop on Computational Advances in 
% Multi-Sensor Adaptive Processing (CAMSAP) (pp. 356-360). IEEE.


function [R,obj,incr,op] = R_Univariate_R1R2(Z,Zphi,lambda,opts)


    % Minimization of the Poisson penalized log-likelood
    %
    %   DKL(Z | R Zphi) + lambda * || (D2 R)_3:T ||_1 + lambda * | (D2 R)_1 - R2/2 + R1/4 | + lambda * | (D2 R)_2 + R2/4 |
    %
    % where DKL stands for the Kullback-Leibler divergence, D2 is the discrete
    % Laplacian operator, || . ||_1 the ell_1-norm defined as the sum of
    % absolute values, I is the indicative function of nonnegative real vectors
    % and lambda_T > 0 is a regularization parameter.
    %
    % The data fidelity term accounts for the epidemiological model proposed by
    % Cori et al., while the regularization term enforces a smooth,
    % piecewise-linear behavior for R(t).
    %
    %
    %
    % Inputs:  - Z: new infection counts
    %          - Zphi: global infectiousness defined as a weighted sum of past counts
    %          - lambda_T: regularization parameter for temporal regularity
    %          - opts: structure containing the properties of the
    %          regularizing functional and the parameters of the
    %          minimization algorithm
    %            - dataterm: 'DKL' (by default)  or 'L2'
    %            - R1: first component of R for each territory, stored as a C x 1 columns vector (default Z(1)/Zphi(1))
    %            - R2: second component of R for each territory, stored as a C x 1 column vector (default Z(2)/Zphi(2))
    %            - Ri: initialization of R (constant equal to 1 by default) of size C x T-2
    %            - iter: maximal number of iterations (1e6 by default)
    %            - incr: 'R' for increments on iterates, 'obj' increments on objective function
    %            - prec: tolerance for the stopping criterion (1e-7 by default)
    %            - stop: 'LimSup' (by default) smoothed increments over win past iterates, or 'Primal' pointwise increments
    %            - win: length of smoothing window (500 by default)
    %
    %
    % Outputs: - R: estimated regularized reproduction number excluding the fixed two first values
    %          - obj: values of the objective function w.r.t iterations
    %          - incr: normalized (smoothed) increments w.r.t iterations
    %          - op: linear direct and adjoint operators involved in the regularization term

    %% NORMALIZE INFECTION COUNTS AND INFECTIOUSNESS

    scale       = std(Z,[],2);   % scale of infection counts
    Z           = Z./scale;
    Zphi        = Zphi./scale;
    
    %% RESIZE INPUT 

    [d1,d2]     = size(Z);

    if min(d1,d2) == 1
        
        Z       = reshape(Z,1,max(d1,d2));
        Zphi    = reshape(Zphi,1,max(d1,d2));

    end

    %% DEFAULTS OPTIONS

    if nargin == 3
        opts     = struct;
    end

    % Regularizing functional
    if ~isfield(opts,'dataterm'),       opts.dataterm       = 'DKL'; end
    if ~isfield(opts,'regularization'), opts.regularization = 'L1'; end % no other option implemented for fixed two first components
    if ~isfield(opts,'prior'),          opts.prior          = 'laplacian'; end % no other option implemented for fixed two first components

    % Default first two values
    if ~isfield(opts,'R1'), opts.R1      = Z(:,1)./Zphi(:,1) ; end
    if ~isfield(opts,'R2'), opts.R2      = Z(:,2)./Zphi(:,2) ; end

    % Ensure no division by zero
    opts.R1(isnan(opts.R1))              = 0;
    opts.R2(isnan(opts.R2))              = 0;

    % Minimization algorithm
    if ~isfield(opts,'Ri'),    opts.Ri   = ones(size(Z(:,3:end))) ; end
    if ~isfield(opts,'iter'),  opts.iter = 1e6 ; end
    if ~isfield(opts,'incr'),  opts.incr = 'R' ; end
    if ~isfield(opts,'prec'),  opts.prec = 1e-7 ; end
    if ~isfield(opts,'stop'),  opts.stop ='LimSup'; end
    if ~isfield(opts,'win'),   opts.win = 500; end

    % Name of the estimator for displaying waiting bar
    opts.flag = 'Univariate, fixed R1, R2 (U-12)';


    %% OBJECTIVE FUNCTION, PROXIMITY OPERATORS AND LINEAR OPERATORS

    % number of days
    T                      = size(Zphi,2);
    Teff                   = T - 2;

    % Data fidelity term
    if strcmp(opts.dataterm,'DKL')

        opts.mu            = 0;
        objective.fidelity = @(y,Zeff) DKLw(y,Zeff,Zphi(:,3:T));
        prox.fidelity      = @(y,Zeff,tau) prox_DKLw(y,Zeff,Zphi(:,3:T),tau);

    elseif strcmp(opts.dataterm,'L2')

        opts.mu = 0;
        objective.fidelity = @(y,Zeff) 0.5*sum(sum((Zeff - Zphi(:,3:T).*y).^2));
        prox.fidelity      = @(y,Zeff,tau) prox_L2w(y,Zeff,Zphi(:,3:T),tau);
        
    end


    % Regularization term
    prox.regularization      = @(y,tau) prox_L1_R1R2(y,tau,lambda*opts.R1,lambda*opts.R2);
    objective.regularization = @(y,tau) tau*sum(abs(y(:,1) - lambda*opts.R2/2 + lambda*opts.R1/4)) + tau*sum(abs(y(:,2) + lambda*opts.R2/4)) + tau*sum(sum(abs(y(:,3:T))));


    % Linear operators
    filter_def   = opts.prior;
    computation  = 'direct';
    param.lambda = lambda;
    param.type   = '1D';
    param.op     = opts.prior;

    op.direct    = @(x) opL_R1R2(x, filter_def, computation, param);
    op.adjoint   = @(x) opLadj_R1R2(x, filter_def, computation, param);
    opts.normL   = lambda^2;
    

    %% RUN THE ALGORITHM AND PREPARE OUTPUTS

    % store the initialization of the primal-dual algorithm in correct form
    opts.xi           = opts.Ri;

    % Minimization of the functional with Chambolle-Pock algorithm
    Zeff              = Z(:,3:T);
    [R,obj,incr]      = PD_ChambollePock_Covid(Zeff, objective, op, prox, opts);

    % Linear operator involved in the regularization
    param.lambda      = 1;
    op.direct         = @(x)opL(x, filter_def, computation, param);
    op.adjoint        = @(x)opLadj(x, filter_def, computation, param);

    % Handle trivial estimates
    for c = 1:size(Z,1)
        if sum(isnan(Z(c,:))) == size(Z,2)
            R(c,:)    = 0;
        end
    end

    % Resize the output to fit input size
    if d2 == 1,     R = reshape(R,d1-2,d2); end

end