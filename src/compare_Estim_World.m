% Compare the new infection counts and estimated reproduction number for
% all countries selected and during the analyzed time period.
%
% Implementation B. Pascal,
% March, 2024


function compare_Estim_World(Countries,Estimates,results)


        % Inputs:  - Countries: list of countries that are to be compared
    %          - Estimates: estimated reproduction numbers for all monitored countries and the different implemented estimators
    %          - results: input data and output estimates for the different estimators computed
    %                     Countries: list of the C countries monitored
    %                     Dates: the T dates corresponding to the considered time period
    %                     Z: infection count times series stored as a matrix of size C x T
    %                     Zphi: global infectiousness times series stored as a matrix of size C x T
    %                     FontSize: desired FontSize for the plots
    %                     Estimates: list of estimators computed
    %                     MLE: Maximum likelihood estimator stored as a matrix of size C x T (if computed)
    %                     Gamma: EpiEstim Mean a posteriori estimator stored as a matrix of size C x T (if computed)
    %                     U: Univariate piecewise linear estimator stored as a matrix of size C x T (if computed)
    %                     U_C: Univariate piecewise linear estimator with sparse corrective terms stored as a matrix of size C x T (if computed)

    if ~isfield(results,'FontSize'), results.FontSize = 22.5; end
    if isempty(Countries),           Countries = results.Countries; end
    
    AllEstimates = ["MLE","Gamma","U","U-C"];
    AllKeys      = ["MLE","Gamma","U","U_C"];
    AllTitles    = ["Maximum Likelihood Estimator (MLE)",
        "Bayesian estimator with Gamma prior (Gamma)",
        "Unviariate estimator (U)",
        "Univariate estimator with correction (U-C)",
        ];

    Colors       = [[0, 0, 1]; % bleu
        [0, 0.5, 0]; % green
        [1, 0.84, 0]; % yellow
        [1, 0.65, 0]; % orange
        [1, 0, 0] % red
        ];

    % Check that all Countries are existing and included in the analysis
    DisplayCountries = [];
    for n = 1:length(Countries)
        
        if ~sum(strcmp(results.Countries,Countries(n)))
            warning(strcat(Countries(n)," will be ignored in the plots, either because it was not found in JHU repository or because it has not been indicated by the user in the preamble."))
        else
            DisplayCountries = [DisplayCountries, Countries(n)];
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
                L = cell(1,min(5,length(DisplayCountries)));
                f                = figure(20+iEst);
                plot(results.Dates,ones(size(results.Z(1,:))),'k-','LineWidth',1);
                grid on; hold on
                for n = 1:min(5,length(DisplayCountries))
                    ind          = find(strcmp(results.Countries,DisplayCountries(n)));
                    q            = plot(results.Dates, results.(AllKeys(iEst))(ind,:),'-','linewidth',2,'color',Colors(n,:)) ;
                    Q            = [Q,q];
                    L{n}         = results.Countries(ind);
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
        if length(DisplayCountries) > 5
            ind          = find(strcmp(results.Countries,DisplayCountries(6)));
            ignored      = results.Countries(ind);
            for n        = 7:length(DisplayCountries)
                ind      = find(strcmp(results.Countries,DisplayCountries(n)));
                ignored  = strcat(ignored,", ",results.Countries(ind));
            end

            if length(ignored) == 1
                warning(strcat("For sake of readability only five countries are compared. ",ignored," has been ignored in the plots."))
            else
                warning(strcat("For sake of readability only five countries are compared. ",ignored," have been ignored in the plots."))
            end
        end


    end

