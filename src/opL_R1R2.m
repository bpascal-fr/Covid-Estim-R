function xt = opL_R1R2(x,filter_def, computation,param)
    % define linear operator associated with the filter in the prior
    %
    % Implementation N. PUSTELNIK, CNRS, ENS Lyon
    % June 2019
    %
    % Updated to exclude the two first components
    % B. Pascal, CNRS, LS2N
    % April 2024

    dim = size(x);
    if strcmp(param.type,'1D')

        if ~exist('computation') %#ok<EXIST>
            computation = 'direct';
        end

        if ~exist('filter_def') %#ok<EXIST>
            filter_def = 'laplacian';
        end


        if strcmp(filter_def,'gradient')
            error('Gradient not available.')
        elseif strcmp(filter_def,'laplacian')
            if numel(param.lambda) == 1
                if strcmp(computation,'fourier')
                    error('Fourier not available use direct computation.')
                else
                    xt = param.lambda*[x(:,3:dim(2))/4 - x(:,2:dim(2)-1)/2 + x(:,1:dim(2)-2)/4,zeros(dim(1),2)];
                end
            else
                error('Multivariate not available.')
            end
        end
        
        xt = [x(:,1)/4, x(:,2)/4 - x(:,1)/2, xt];
    else
        error('Not relevant for 2D data.')
    end
