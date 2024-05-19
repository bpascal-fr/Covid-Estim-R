% Display the new infection counts and estimated reproduction number for
% all French departments selected and during the analyzed time period.
%
% Implementation B. Pascal,
% March, 2024


function display_Estim_France(Departments,Estimates,results)


    % Inputs:  - Departments: list of French departments referred to by their INSEE number that are to be displayed
    %          - Estimates: estimated reproduction numbers for monitored French departments and the different implemented estimators
    %          - results: input data 
    %                     Departments: list of the D French departments Names and associated Numbers and Indices
    %                     Dates: the T dates corresponding to the considered time period
    %                     Z: infection count times series stored as a matrix of size D x T
    %                     Zphi: global infectiousness times series stored as a matrix of size D x T
    %                     FontSize: desired FontSize for the plots
    %                     Estimates: list of estimators computed
    %                     MLE: Maximum likelihood estimator stored as a matrix of size D x T (if computed)
    %                     Gamma: Maximum a posteriori estimator stored as a matrix of size D x T (if computed)
    %                     U: Univariate piecewise linear in time estimator stored as a matrix of size D x T (if computed)
    %                     U_C: Univariate piecewise linear in time estimator with sparse corrective terms stored as a matrix of size D x T (if computed)
    %                     M: Multivariate piecewise linear in time and piecewise constant in space estimator stored as a matrix of size D x T (if computed)
    %                     M_C: Multivariate piecewise linear in time and piecewise constant in space  estimator with sparse corrective terms stored as a matrix of size D x T (if computed)

    if isempty(Departments),         Departments      = results.Departments.Numbers; end
    if ~isfield(results,'FontSize'), results.FontSize = 22.5; end
    if ~isfield(results,'Dates'),    results.Dates    = 1:size(Z,2); end
    

    AllEstimates = ["MLE","Gamma","U","U-C","M","M-C"];

    % By default plot all available estimates and discard invalid estimates
    if isempty(Estimates)
        Estimates = results.Estimates;
    else
        for est = Estimates
            if isempty(find(strcmp(AllEstimates,est),1)) 
                warning(strcat("Estimator ",est," will be ignored because it is not a valid estimator name."))
            else
                if isempty(find(strcmp(results.Estimates,est),1))
                    warning(strcat("Estimator ",est," will be ignored because it has not been computed."))
                end
            end
        end
    end


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



    for d = 1:length(Departments)

        if ~sum(strcmp(results.Departments.Numbers,Departments(d)))
            warning(strcat(Departments(d)," will be ignored in the plots, either because it is not a valid French department INSEE number or because it was not included in the analysis."))
        else
            ind             = find(strcmp(results.Departments.Numbers,Departments(d)));
            f2              = figure(2000 + d); clf
            plot(results.Dates,ones(size(results.Z(ind,:))),'k-','LineWidth',1);
            grid on; hold on
            Q               = [];
            iEst            = 1;
            if ~isempty(find(strcmp(Estimates,'MLE'),1)) && ~isempty(find(strcmp(results.Estimates,'MLE'),1))
                q1          = plot(results.Dates, results.MLE(ind,:),'-','linewidth',2,'color',[0.7,0.7,0.7]) ;
                Q           = [Q, q1];
                L{iEst}     = '$\mathrm{R}^{\mathrm{MLE}}_t$';
                iEst        = iEst + 1;
            end
            if ~isempty(find(strcmp(Estimates,'Gamma'),1)) && ~isempty(find(strcmp(results.Estimates,'Gamma'),1))
                q2          = plot(results.Dates, results.Gamma(ind,:),'-','linewidth',2,'color',[0,0.5,0]) ;
                Q           = [Q, q2];
                L{iEst}     = '$\mathrm{R}_{\tau,t}^\Gamma$';
                iEst        = iEst + 1;
            end
            if ~isempty(find(strcmp(Estimates,'U'),1)) && ~isempty(find(strcmp(results.Estimates,'U'),1))
                q3          = plot(results.Dates, results.U(ind,:),'--','linewidth',2,'color','blue') ;
                Q           = [Q, q3];
                L{iEst}     = '$\mathrm{R}^{\mathrm{U}}_t$';
                iEst        = iEst + 1;
            end
            if ~isempty(find(strcmp(Estimates,'U-C'),1)) && ~isempty(find(strcmp(results.Estimates,'U-C'),1))
                q4          = plot(results.Dates, results.U_C(ind,:),'--','linewidth',2,'color','red') ;
                Q           = [Q, q4];
                L{iEst}     = '$\mathrm{R}^{\mathrm{U-C}}_t$';
                iEst        = iEst + 1;
            end
            if ~isempty(find(strcmp(Estimates,'M'),1)) && ~isempty(find(strcmp(results.Estimates,'M'),1))
                q5          = plot(results.Dates, results.M(ind,:),'-','linewidth',2,'color','blue') ;
                Q           = [Q, q5];
                L{iEst}     = '$\mathrm{R}^{\mathrm{M}}_t$';
                iEst        = iEst + 1;
            end
            if ~isempty(find(strcmp(Estimates,'M-C'),1)) && ~isempty(find(strcmp(results.Estimates,'M-C'),1))
                q6          = plot(results.Dates, results.M_C(ind,:),'-','linewidth',2,'color','red') ;
                Q           = [Q, q6];
                L{iEst}     = '$\mathrm{R}^{\mathrm{M-C}}_t$';
                iEst        = iEst + 1;
            end
            if isempty(Q)

                close(2000 + d)
                warning('No valid estimator name in the list. R estimates not displayed. Valid estimators names are: MLE, Gamma, U, U-C, M and M-C.')

            else

                leg             = legend(Q,L,'location','best');
                leg.Interpreter = 'Latex'; leg.Color ='none'; leg.FontSize = results.FontSize;
                set(gca,'ticklabelinterpreter','Latex','fontsize',results.FontSize,'color','None')
                VX              = [results.Dates(1) results.Dates(end)] ; xlim(VX)
                VY              = [0, 3]; ylim(VY)
                title(results.Departments.Names(ind),'Interpreter','Latex')
                ylabel('reproduction number','Interpreter','Latex')
                f2.Position     = [211 287 943 314];



                f3 = figure(3000 + d); clf
                subplot(211)
                p               = plot(results.Dates, results.Z(ind,:),'-','linewidth',2,'color','black') ;
                grid on ; hold on
                q               = plot(results.Dates, results.Zphi(ind,:),'-.','linewidth',2,'color','black') ;
                leg             = legend([p,q],'$\mathrm{Z}_t$','$\Phi_t^{\mathrm{Z}}$','location','best');
                leg.Interpreter = 'Latex';
                leg.Color       = 'none';
                leg.FontSize    = results.FontSize;
                VX              = [results.Dates(1) results.Dates(end)] ;
                xlim(VX)
                xticklabels({})
                title(results.Departments.Names(ind),'Interpreter','Latex')
                ylabel('infection counts','Interpreter','Latex')
                set(gca,'ticklabelinterpreter','Latex','fontsize',results.FontSize,'color','None')
                subplot(212)
                plot(results.Dates,ones(size(results.Z(ind,:))),'k-','LineWidth',1);
                grid on; hold on
                Q               = [];
                iEst            = 1;
                if ~isempty(find(strcmp(Estimates,'MLE'),1))  && ~isempty(find(strcmp(results.Estimates,'MLE'),1))
                    q1          = plot(results.Dates, results.MLE(ind,:),'-','linewidth',2,'color',[0.7,0.7,0.7]) ;
                    Q           = [Q, q1];
                    L{iEst}     = '$\mathrm{R}^{\mathrm{MLE}}_t$';
                    iEst        = iEst + 1;
                end
                if ~isempty(find(strcmp(Estimates,'Gamma'),1)) && ~isempty(find(strcmp(results.Estimates,'Gamma'),1))
                    q2          = plot(results.Dates, results.Gamma(ind,:),'-','linewidth',2,'color',[0,0.5,0]) ;
                    Q           = [Q, q2];
                    L{iEst}     = '$\mathrm{R}_{\tau,t}^\Gamma$';
                    iEst        = iEst + 1;
                end
                if ~isempty(find(strcmp(Estimates,'U'),1)) && ~isempty(find(strcmp(results.Estimates,'U'),1))
                    q3          = plot(results.Dates, results.U(ind,:),'--','linewidth',2,'color','blue') ;
                    Q           = [Q, q3];
                    L{iEst}     = '$\mathrm{R}^{\mathrm{U}}_t$';
                    iEst        = iEst + 1;
                end
                if ~isempty(find(strcmp(Estimates,'U-C'),1)) && ~isempty(find(strcmp(results.Estimates,'U-C'),1))
                    q4          = plot(results.Dates, results.U_C(ind,:),'--','linewidth',2,'color','red') ;
                    Q           = [Q, q4];
                    L{iEst}     = '$\mathrm{R}^{\mathrm{U-C}}_t$';
                    iEst        = iEst + 1;
                end
                if ~isempty(find(strcmp(Estimates,'M'),1)) && ~isempty(find(strcmp(results.Estimates,'M'),1))
                    q5          = plot(results.Dates, results.M(ind,:),'-','linewidth',2,'color','blue') ;
                    Q           = [Q, q5];
                    L{iEst}     = '$\mathrm{R}^{\mathrm{M}}_t$';
                    iEst        = iEst + 1;
                end
                if ~isempty(find(strcmp(Estimates,'M-C'),1)) && ~isempty(find(strcmp(results.Estimates,'M-C'),1))
                    q6          = plot(results.Dates, results.M_C(ind,:),'-','linewidth',2,'color','red') ;
                    Q           = [Q, q6];
                    L{iEst}     = '$\mathrm{R}^{\mathrm{M-C}}_t$';
                    iEst        = iEst + 1;
                end
                leg             = legend(Q,L,'location','best');
                leg.Interpreter = 'Latex'; leg.Color ='none'; leg.FontSize = results.FontSize;
                set(gca,'ticklabelinterpreter','Latex','fontsize',results.FontSize,'color','None')
                VX              = [results.Dates(1) results.Dates(end)] ; xlim(VX)
                VY              = [0, 3]; ylim(VY)
                ylabel('reproduction number','Interpreter','Latex')
                f3.Position     = [211 87 943 710];
            end
        end

    end

end
