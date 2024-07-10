% Compare the new infection counts and estimated reproduction number for
% all French departments selected and during the analyzed time period.
%
% Implementation B. Pascal,
% March, 2024


function compare_Estim_France(Departments,Estimates,results)


    % Inputs:  - Departments: list of French departments referred to by their INSEE number that are to be compared
    %          - Estimates: estimated reproduction numbers for monitored French departments and the different implemented estimators
    %          - results: input data
    %                     Departments: list of the D French departments Names and associated Numbers and Indices
    %                     Dates: the T dates corresponding to the considered time period
    %                     Z: infection count times series stored as a matrix of size D x T
    %                     Zphi: global infectiousness times series stored as a matrix of size D x T
    %                     FontSize: desired FontSize for the plots
    %                     Estimates: list of estimators computed
    %                     MLE: Maximum likelihood estimator stored as a matrix of size D x T (if computed)
    %                     Gamma: EpiEstim Mean a posteriori estimator stored as a matrix of size D x T (if computed)
    %                     U: Univariate piecewise linear in time estimator stored as a matrix of size D x T (if computed)
    %                     U_C: Univariate piecewise linear in time estimator with sparse corrective terms stored as a matrix of size D x T (if computed)
    %                     M: Multivariate piecewise linear in time and piecewise constant in space estimator stored as a matrix of size D x T (if computed)
    %                     M_C: Multivariate piecewise linear in time and piecewise constant in space  estimator with sparse corrective terms stored as a matrix of size D x T (if computed)

    if ~isfield(results,'FontSize'), results.FontSize = 22.5; end
    if isempty(Departments),         Departments = results.Departments.Numbers; end

    AllEstimates = ["MLE","Gamma","U","U-C","M","M-C"];
    AllKeys      = ["MLE","Gamma","U","U_C","M","M_C"];
    AllTitles    = ["Maximum Likelihood Estimator (MLE)",
        "Bayesian estimator with Gamma prior (Gamma)",
        "Unviariate estimator (U)",
        "Univariate estimator with correction (U-C)",
        "Multiariate estimator (M)",
        "Multivariate estimator with correction (M-C)"
        ];

    Colors       = [[0, 0, 1]; % bleu
        [0, 0.5, 0]; % green
        [1, 0.84, 0]; % yellow
        [1, 0.65, 0]; % orange
        [1, 0, 0] % red
        ];

    % Check that all Departments are existing and included in the analysis
    DisplayDepartments = [];
    for d = 1:length(Departments)
        
        if ~sum(strcmp(results.Departments.Numbers,Departments(d)))
            warning(strcat(Departments(d)," will be ignored in the plots, either because it is not a valid French department INSEE number or because it was not included in the analysis."))
        else
            DisplayDepartments = [DisplayDepartments, Departments(d)];
        end
    end

    % By default plot all available estimates and discard invalid estimates
    valid_est    = 0;
    if isempty(Estimates)
        Estimates = results.Estimates;
    else
        for est = Estimates
            if isempty(find(strcmp(AllEstimates,est),1))
                warning(strcat("Estimator ",est," will be ignored because it is not a valid estimator name."))
            else
                if isempty(find(strcmp(results.Estimates,est),1))
                    warning(strcat("Estimator ",est," will be ignored because it has not been computed."))
                else
                    valid_est = valid_est + 1;
                end
            end
        end
    end

    if ~valid_est

        warning('No valid estimator name in the list. R estimates not displayed. Valid estimators names are: MLE, Gamma, U, U-C, M and M-C.')

    else

        % Manage the case when Z is univariate and stored in a columns vector of size T x 1 instead of 1 x T
        [d1,d2] = size(results.Z);
        if min(d1,d2) == 1
            results.Z    = reshape(results.Z,1,max(d1,d2));
            results.Zphi = reshape(results.Zphi,1,max(d1,d2));
            if isfield(results,'MLE'),   results.MLE = reshape(results.MLE,1,max(d1,d2));     end
            if isfield(results,'Gamma'), results.Gamma = reshape(results.Gamma,1,max(d1,d2)); end
            if isfield(results,'U'),     results.U = reshape(results.U,1,max(d1,d2));         end
            if isfield(results,'U_C'),   results.U_C = reshape(results.U_C,1,max(d1,d2));     end
            if isfield(results,'M'),     results.M = reshape(results.M,1,max(d1,d2));         end
            if isfield(results,'M_C'),   results.M_C = reshape(results.M_C,1,max(d1,d2));     end
        end

        iEst      = 0;
        for Estim = AllEstimates

            iEst = iEst + 1;
            ine  = find(strcmp(Estim,Estimates),1);

            if ~isempty(ine) && ~isempty(find(strcmp(results.Estimates,Estim),1))

                Q = [];
                L = cell(1,min(5,length(DisplayDepartments)));
                f                = figure(20+iEst);
                plot(results.Dates,ones(size(results.Z(1,:))),'k-','LineWidth',1);
                grid on; hold on
                for d = 1:min(5,length(DisplayDepartments))
                    ind          = find(strcmp(results.Departments.Numbers,DisplayDepartments(d)));
                    q            = plot(results.Dates, results.(AllKeys(iEst))(ind,:),'-','linewidth',2,'color',Colors(d,:)) ;
                    Q            = [Q,q];
                    L{d}         = results.Departments.Names(ind);
                end
                leg              = legend(Q,L,'location','best');
                leg.Interpreter  = 'Latex'; leg.Color ='none'; leg.FontSize = results.FontSize;
                set(gca,'ticklabelinterpreter','Latex','fontsize',results.FontSize,'color','None')
                VX               = [results.Dates(1) results.Dates(end)] ; xlim(VX)
                VY               = [0, 3]; ylim(VY)
                title(AllTitles(iEst),'Interpreter','Latex')
                ylabel('reproduction number','Interpreter','Latex')
                f.Position       = [211 287 943 314];

            end

        end

        % Do not display more than five departments simultaneously
        if length(DisplayDepartments) > 5
            ind          = find(strcmp(results.Departments.Numbers,DisplayDepartments(6)));
            ignored      = results.Departments.Names_Table(ind);
            for d        = 7:length(DisplayDepartments)
                ind      = find(strcmp(results.Departments.Numbers,DisplayDepartments(d)));
                ignored  = strcat(ignored,", ",results.Departments.Names_Table(ind));
            end

            if length(ignored) == 1
                warning(strcat("For sake of readability only five departments are compared. ",ignored," has been ignored in the plots."))
            else
                warning(strcat("For sake of readability only five departments are compared. ",ignored," have been ignored in the plots."))
            end
        end


    end

