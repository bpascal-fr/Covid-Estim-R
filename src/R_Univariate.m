% Estimation of the daily reproduction number R(t) of an epidemic from new
% infection counts time series Z(t).
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


function [R,obj,incr,op] = R_Univariate(Z,Zphi,lambda,opts)


    % Minimization of the Poisson penalized log-likelood
    %
    %       DKL(Z | R Zphi) + lambda * || D2 R ||_1
    %
    % where DKL stands for the Kullback-Leibler divergence, D2 is the discrete
    % Laplacian operator, || . ||_1 the ell_1-norm defined as the sum of
    % absolute values and lambda > 0 is a regularization parameter.
    %
    % The data fidelity term accounts for the epidemiological model proposed by
    % Cori et al., while the regularization term enforces a smooth,
    % piecewise-linear behavior for R(t).
    %
    %
    %
    % Inputs:  - Z: new infection counts
    %          - Zphi: global infectiousness defined as a weighted sum of past counts
    %          - lambda: regularization parameter
    %          - opts: structure containing the properties of the
    %          regularizing functional and the parameters of the
    %          minimization algorithm
    %            - dataterm: 'DKL' (by default)  or 'L2'
    %            - regularization: 'L1' (by default) or 'L12'
    %            - prior: 'laplacian' (by default) or 'gradient'
    %            - Ri: initialization of R (constant equal to 1 by default)
    %            - strat: 'usual' (by default) or 'accelerated' (only for L2)
    %            - iter: maximal number of iterations (1e6 by default)
    %            - incr: 'R' for increments on iterates, 'obj' increments on objective function
    %            - prec: tolerance for the stopping criterion (1e-7 by default)
    %            - stop: 'LimSup' (by default) smoothed increments over win past iterates, or 'Primal' pointwise increments
    %            - win: length of smoothing window (500 by default)
    %
    %
    % Outputs: - R: estimated regularized reproduction number
    %          - obj: values of the objective function w.r.t iterations
    %          - incr: normalized (smoothed) increments w.r.t iterations
    %          - op: linear direct and adjoint operators involved in the regularization term
    
    %% DEFAULTS OPTIONS

    if nargin == 3
        opts     = struct;
    end

    % Regularizing functional
    if ~isfield(opts,'dataterm'),       opts.dataterm       = 'DKL'; end
    if ~isfield(opts,'regularization'), opts.regularization = 'L1'; end
    if ~isfield(opts,'prior'),          opts.prior          = 'laplacian'; end

    % Minimization algorithm
    if ~isfield(opts,'Ri'),    opts.Ri   = ones(size(Z)) ; end
    if ~isfield(opts,'strat'), opts.strat = 'usual' ; end
    if ~isfield(opts,'iter'),  opts.iter = 1e6 ; end
    if ~isfield(opts,'incr'),  opts.incr = 'R' ; end
    if ~isfield(opts,'prec'),  opts.prec = 1e-7 ; end
    if ~isfield(opts,'stop'),  opts.stop ='LimSup'; end
    if ~isfield(opts,'win'),   opts.win = 500; end

    % Name of the estimator for displaying waiting bar
    opts.flag = 'Univariate (U)';
    
    %% RESIZE INPUT 

    [d1,d2]     = size(Z);

    if min(d1,d2) == 1
        
        Z       = reshape(Z,1,max(d1,d2));
        Zphi    = reshape(Zphi,1,max(d1,d2));
        opts.Ri = reshape(opts.Ri,1,max(d1,d2));

    end

    %% OBJECTIVE FUNCTION, PROXIMITY OPERATORS AND LINEAR OPERATORS

    % Data fidelity term
    if strcmp(opts.dataterm,'DKL')

        opts.mu            = 0;
        objective.fidelity = @(y,Z) DKLw(y,Z,Zphi);
        prox.fidelity      = @(y,Z,tau) prox_DKLw(y,Z,Zphi,tau);

    elseif strcmp(opts.dataterm,'L2')

        if strcmp(opts.strat,'accelerated')
            opts.mu = min(Zphi.^2);
        elseif strcmp(opts.strat,'usual')
            opts.mu = 0;
        end

        objective.fidelity = @(y,Z) 0.5*sum((Z(:) - Zphi(:).*y(:)).^2);
        prox.fidelity      = @(y,Z,tau) prox_L2w(y,Z,Zphi,tau);

    end


    % Regularization term
    if strcmp(opts.regularization,'L1')

        prox.regularization      = @(y,tau) prox_L1(y,tau);
        objective.regularization = @(y,tau) tau*sum(abs(y(:)));

    elseif strcmp(opts.regularization,'L12')

        prox.regularization      = @(y,tau) prox_L12(y,tau);
        objective.regularization = @(y,tau) tau*sum(sqrt(sum(y.^2,1)));

    end


    % Linear operators
    filter_def   = opts.prior;
    computation  = 'direct';
    param.lambda = lambda;
    param.type   = '1D';
    param.op     = opts.prior;

    op.direct    = @(x)opL(x, filter_def, computation, param);
    op.adjoint   = @(x)opLadj(x, filter_def, computation, param);
    opts.normL   = lambda^2;
    

    %% RUN THE ALGORITHM AND PREPARE OUTPUTS

    % store the initialization of the primal-dual algorithm in correct form
    opts.xi             = opts.Ri;

    % Minimization of the functional with Chambolle-Pock algorithm
    [R,obj,incr]      = PD_ChambollePock_Covid(Z, objective, op, prox, opts);

    % Linear operator involved in the regularization
    param.lambda       = 1;
    op.direct          = @(x)opL(x, filter_def, computation, param);
    op.adjoint         = @(x)opLadj(x, filter_def, computation, param);

    % Resize the output to fit input size
    R                 = reshape(R,d1,d2);

end