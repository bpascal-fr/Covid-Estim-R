function x = opLadj_R1R2(y,filter_def, computation,param)
    % Adjoint of the operator used in the penalization
    %
    % Implementation N. PUSTELNIK, CNRS, ENS Lyon
    % June 2019
    %
    % Updated to fix the two first components
    % B. Pascal, CNRS, LS2N
    % April 2024

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
                    y1 = y(:,1:2);
                    y  = y(:,3:end);
                    dim = size(y);
                    x  = param.lambda*([y(:,1:1)/4,-y(:,1)/2 + y(:,2)/4, y(:,3:dim(2)-2)/4 - y(:,2:dim(2)-3)/2 + y(:,1:dim(2)-4)/4, y(:,dim(2)-3)/4 - y(:,dim(2)-2)/2, y(:,dim(2)-2)/4]);
                    x(1)  = x(1) + param.lambda*(y1(:,1)/4 - y1(:,2)/2);
                    x(2)  = x(2) + param.lambda*(y1(:,2)/4);
                end
            else
                error('Multivariate not available.')
            end
        end


    else
        error('Not relevant for 2D data.')
    end