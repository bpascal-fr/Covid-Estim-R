% Multivariate simulatenous estimation of the daily reproduction number R(c,t) 
% of an epidemic and of corrective terms modeling misreported counts from new 
% infection counts time series Z(c,t) in C connected territories including
% a spatial graph-inspired regularization.
%
% from
% - Cori, A., Ferguson, N. M., Fraser, C., & Cauchemez, S. (2013).
% A new framework and software to estimate time-varying reproduction
% numbers during epidemics. American journal of epidemiology, 178(9),
% 1505-1512.
%
% - Pascal, B., Abry, P., Pustelnik, N., Roux, S., Gribonval, R., 
% & Flandrin, P. (2022). Nonsmooth convex optimization to estimate the 
% Covid-19 reproduction number space-time evolution with robustness against 
% low quality data. IEEE Transactions on Signal Processing, 70, 2859–2868.


function [R,O,obj,incr,op] = R_Multivariate_Correct(Z,Zphi,lambda_T,G,lambda_S,lambda_O,opts)


    % Minimization of the Poisson penalized log-likelood
    %
    %   DKL(Z | R Zphi + O) + lambda_T * || D2 R ||_1 + lambda_S * || G R ||_1 + lambda_O * || O ||_1 + I( R >= 0 )
    %
    % where DKL stands for the Kullback-Leibler divergence, D2 is the discrete
    % Laplacian operator, G is the discrete gradient operator on the graph 
    % of connected territories, || . ||_1 the ell_1-norm defined as the sum of
    % absolute values, I is the indicative function of nonnegative real vectors
    % and lambda_T, lambda_S, lambda_O > 0 are regularization parameters.
    %
    % The data fidelity term accounts for the epidemiological model proposed by
    % Cori et al., while the regularization term enforces a smooth,
    % piecewise linear temporal behavior and piecewise constant spatial 
    % behavior of the estimated R and sparsity of O.
    %
    %
    %
    % Inputs:  - Z: new infection counts
    %          - Zphi: global infectiousness defined as a weighted sum of past counts
    %          - lambda_T: regularization parameter for temporal regularity
    %          - G: matrix of the discrete gradient on the graph of connected territories
    %          - lambda_S: regularization parameter for spatial piecewise constancy of estimated reproduction numbers
    %          - lambda_O: regularization parameter for sparsity of the corrective term
    %          - opts: structure containing the properties of the
    %          regularizing functional and the parameters of the
    %          minimization algorithm
    %            - dataterm: 'DKL' (by default)  or 'L2'
    %            - regularization: 'L1' (by default) or 'L12'
    %            - prior: 'laplacian' (by default) or 'gradient'
    %            - spatial: 'L1' (by default) or 'L2' on G R
    %            - Ri: initialization of R (constant equal to 1 by default)
    %            - Oi: initialization of O (constant equal to 0 by default)
    %            - iter: maximal number of iterations (1e6 by default)
    %            - incr: 'R' for increments on iterates, 'obj' increments on objective function
    %            - prec: tolerance for the stopping criterion (1e-7 by default)
    %            - stop: 'LimSup' (by default) smoothed increments over win past iterates, or 'Primal' pointwise increments
    %            - win: length of smoothing window (500 by default)
    %            - flag: if 'none' not progression bar (optional)
    %
    %
    % Outputs: - R: estimated regularized reproduction number
    %          - O: estimated correction of misreported counts
    %          - obj: values of the objective function w.r.t iterations
    %          - incr: normalized (smoothed) increments w.r.t iterations
    %          - op: linear direct and adjoint operators invovled in the regularization term

    %% DEFAULTS OPTIONS

    if nargin == 6
        opts     = struct;
    end

    % Regularizing functional
    if ~isfield(opts,'dataterm'),       opts.dataterm       = 'DKL'; end
    if ~isfield(opts,'regularization'), opts.regularization = 'L1'; end
    if ~isfield(opts,'prior'),          opts.prior          = 'laplacian'; end
    if ~isfield(opts,'spatial'),        opts.spatial        = 'L1'; end

    % Minimization algorithm
    if ~isfield(opts,'Ri'),    opts.Ri   = ones(size(Z)) ; end
    if ~isfield(opts,'Oi'),    opts.Oi   = zeros(size(Z)) ; end
    if ~isfield(opts,'iter'),  opts.iter = 1e6 ; end
    if ~isfield(opts,'incr'),  opts.incr = 'R' ; end
    if ~isfield(opts,'prec'),  opts.prec = 1e-7 ; end
    if ~isfield(opts,'stop'),  opts.stop ='LimSup'; end
    if ~isfield(opts,'win'),   opts.win = 500; end

    % Name of the estimator for displaying waiting bar
    if isfield(opts.flag)
        if ~strcmp(opts.flag,'none')
            opts.flag = 'Multivariate corrected (M-C)';
        else
            opts = rmfield(opts,'flag');
        end
    end
    
    %% NORMALIZE INFECTION COUNTS AND INFECTIOUSNESS

    scale       = std(Z,[],2);   % scale of infection counts
    Z           = Z./scale;
    Zphi        = Zphi./scale;

    %% RESIZE INPUT 

    [d1,d2]     = size(Z);

    if min(d1,d2) == 1
        
        Z       = reshape(Z,1,max(d1,d2));
        Zphi    = reshape(Zphi,1,max(d1,d2));
        opts.Ri = reshape(opts.Ri,1,max(d1,d2));
        opts.Oi = reshape(opts.Oi,1,max(d1,d2));

    end

    %% OBJECTIVE FUNCTION, PROXIMITY OPERATORS AND LINEAR OPERATORS

    % number of days
    [C,T]                  = size(Zphi);

    % Data fidelity term
    if strcmp(opts.dataterm,'DKL')

        opts.mu            = 0;
        objective.fidelity = @(y,Z) DKLw_outlier(y,Z,Zphi);
        prox.fidelity      = @(y,Z,tau) prox_DKLw_outlier(y,Z,Zphi,tau);

    elseif strcmp(opts.dataterm,'L2')

        opts.mu = 0;
        objective.fidelity = @(y,Z) 0.5*sum(sum((Z - Zphi.*y(:,1:T) - y(:,T+1:end) ).^2));
        prox.fidelity      = @(y,Z,tau) prox_L2w_outlier(y,Z,Zphi,tau);
        
    end


    % Regularization term
    if strcmp(opts.regularization,'L1')

        if strcmp(opts.spatial,'L1')
            prox.regularization      = @(y,tau) [prox_L1(y(1:2*C,:),tau); max(y(2*C+1:3*C,:),0); prox_L1(y(3*C+1:end,:),tau)];
            objective.regularization = @(y,tau) tau*sum(sum(abs([y(1:2*C,:);y(3*C+1:end,:)])));
        elseif strcmp(opts.spatial,'L2')
            prox.regularization      = @(y,tau) [prox_L1(y(1:C,:),tau); prox_L2(y(C+1:2*C,:),tau); max(y(2*C+1:3*C,:),0); prox_L1(y(3*C+1:end,:),tau)];
            objective.regularization = @(y,tau) tau*sum(sum(abs([y(1:C,:);y(3*C+1:end,:)]))) + tau*sum(sum(y(C+1:2*C,:).^2));
        end

    elseif strcmp(opts.regularization,'L12')

        if strcmp(opts.spatial,'L1')
            prox.regularization      = @(y,tau) [prox_L12(y(1:C,:),tau); prox_L1(y(C+1:2*C,:),tau); max(y(2*C+1:3*C,:),0); prox_L1(y(3*C+1:end,:),tau)];
            objective.regularization = @(y,tau) tau*sum(sqrt(sum(y(1:C,:).^2,1))) + tau*sum(sum(abs([y(C+1:2*C,:);y(3*C+1:end,:)])));
        elseif strcmp(opts.spatial,'L2')
            prox.regularization      = @(y,tau) [prox_L12(y(1:C,:),tau); prox_L2(y(C+1:2*C,:),tau); max(y(2*C+1:3*C,:),0); prox_L1(y(3*C+1:end,:),tau)];
            objective.regularization = @(y,tau) tau*sum(sqrt(sum(y(1:C,:).^2,1))) + tau*sum(sum(y(C+1:2*C,:).^2)) + tau*sum(sum(abs(y(3*C+1:end,:))));
        end

    end


    % Linear operators
    filter_def   = opts.prior;
    computation  = 'direct';
    param.lambda = lambda_T;
    param.type   = '1D';
    param.op     = opts.prior;

    op.direct    = @(x) [opL(x(:,1:T), filter_def, computation, param); lambda_O * x(:,T+1:end); x(:,1:T); lambda_S * G * x(:,1:T)];
    op.adjoint   = @(x) [opLadj(x(1:C,:), filter_def, computation, param) + x(2*C+1:3*C,:) + lambda_S * G' * x(3*C+1:end,:), lambda_O * x(C+1:2*C,:)];
    opts.normL   = max(lambda_T^2 + 1 + lambda_S^2 * norm(G)^2, lambda_O^2);
    

    %% RUN THE ALGORITHM AND PREPARE OUTPUTS

    % store the initialization of the primal-dual algorithm in correct form
    opts.xi           = [opts.Ri, opts.Oi];

    % Minimization of the functional with Chambolle-Pock algorithm
    [x,obj,incr]      = PD_ChambollePock_Covid(Z, objective, op, prox, opts);

    % Linear operator involved in the regularization
    param.lambda      = 1;
    op.direct         = @(x)opL(x, filter_def, computation, param);
    op.adjoint        = @(x)opLadj(x, filter_def, computation, param);

    % Handle trivial estimates
    for c = 1:size(Z,1)
        if sum(isnan(Z(c,:))) == size(Z,2)
            x(c,:)    = 0;
        end
    end


    % Resize the output to fit input size
    R                 = reshape(x(:,1:T),d1,d2); 
    O                 = reshape(x(:,T+1:end),d1,d2).*scale; 

end