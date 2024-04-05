% Display the new infection counts and estimated reproduction number for
% all countries selected and during the analyzed time period.
%
% Implementation B. Pascal,
% March, 2024


function display_Estim_World_R1R2(Countries,Estimates,results)


    % Inputs:  - Countries: list of countries that are to be displayed
    %          - Estimates: estimated reproduction numbers for all monitored countries and the different implemented estimators
    %          - results: input data and output estimates for the different estimators computed
    %                     Countries: list of the C countries monitored
    %                     Dates: the T dates corresponding to the considered time period
    %                     Z: infection count times series stored as a matrix of size C x T
    %                     Zphi: global infectiousness times series stored as a matrix of size C x T
    %                     FontSize: desired FontSize for the plots
    %                     Estimates: list of estimators computed
    %                     MLE: Maximum likelihood estimator stored as a matrix of size C x T (if computed)
    %                     Gamma: Maximum a posteriori estimator stored as a matrix of size C x T (if computed)
    %                     U: Univariate piecewise linear estimator stored as a matrix of size C x T (if computed)
    %                     U_C: Univariate piecewise linear estimator with sparse corrective terms stored as a matrix of size C x T (if computed)

    if ~isfield(results,'FontSize'), results.FontSize = 22.5; end
    if isempty(Countries),           Countries = results.Countries; end
    
    AllEstimates = ["MLE","Gamma","U-12","U-C-12"];

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
        if isfield(results,'MLE'),      results.MLE = reshape(results.MLE,1,max(d1,d2));     end
        if isfield(results,'Gamma'),    results.Gamma = reshape(results.Gamma,1,max(d1,d2)); end
        if isfield(results,'U_12'),     results.U = reshape(results.U_12,1,max(d1,d2));      end
        if isfield(results,'U_C_12'),   results.U_C = reshape(results.U_C_12,1,max(d1,d2));  end
    end

    % Number of days
    T     = size(Z,2);
    Teff  = T-2;

    for n = 1:length(Countries)

        if ~sum(contains(results.Countries,Countries(n)))
            warning(strcat(Countries(n)," was not analysed, either because it was not found in JHU repository or because it has not been indicated by the user in the preamble. Hence, ",Countries(n)," will not appear in the plots."))

        else

            f2              = figure(2000 + n); clf
            plot(results.Dates,ones(size(results.Z(n,:))),'k-','LineWidth',1);
            grid on; hold on
            Q               = [];
            iEst            = 1;
            if ~isempty(find(strcmp(Estimates,'MLE'),1))  && ~isempty(find(strcmp(results.Estimates,'MLE'),1))
                q1          = plot(results.Dates, results.MLE(n,:),'-','linewidth',2,'color',[0.7,0.7,0.7]) ;
                Q           = [Q, q1];
                L{iEst}     = '$\mathrm{R}^{\mathrm{MLE}}_t$';
                iEst        = iEst + 1;
            end
            if ~isempty(find(strcmp(Estimates,'Gamma'),1)) && ~isempty(find(strcmp(results.Estimates,'Gamma'),1))
                q2          = plot(results.Dates, results.Gamma(n,:),'-','linewidth',2,'color',[0,0.5,0]) ;
                Q           = [Q, q2];
                L{iEst}     = '$\mathrm{R}_{\tau,t}^\Gamma$';
                iEst        = iEst + 1;
            end
            if ~isempty(find(strcmp(Estimates,'U-12'),1)) && ~isempty(find(strcmp(results.Estimates,'U-12'),1))
                q3          = plot(results.Dates(3:Teff), results.U_12(n,:),'-','linewidth',2,'color','blue') ;
                Q           = [Q, q3];
                L{iEst}     = '$\mathrm{R}^{\mathrm{U-12}}_t$';
                iEst        = iEst + 1;
            end
            if ~isempty(find(strcmp(Estimates,'U-C-12'),1)) && ~isempty(find(strcmp(results.Estimates,'U-C-12'),1))
                q4          = plot(results.Dates(3:Teff), results.U_C(n,:),'-','linewidth',2,'color','red') ;
                Q           = [Q, q4];
                L{iEst}     = '$\mathrm{R}^{\mathrm{U-C-12}}_t$';
                iEst        = iEst + 1;
            end
            if isempty(Q)

                close(2000 + n)
                warning('No valid estimator name in the list. R estimates not displayed. Valid estimators names are: MLE, Gamma, U and U-C.')

            else

                leg             = legend(Q,L,'location','best');
                leg.Interpreter = 'Latex'; leg.Color ='none'; leg.FontSize = results.FontSize;
                set(gca,'ticklabelinterpreter','Latex','fontsize',results.FontSize,'color','None')
                VX              = [results.Dates(1) results.Dates(end)] ; xlim(VX)
                VY              = [0, 3]; ylim(VY)
                title(Countries(n),'Interpreter','Latex')
                ylabel('reproduction number','Interpreter','Latex')
                f2.Position     = [211 287 943 314];



                f3 = figure(3000 + n); clf
                subplot(211)
                p               = plot(results.Dates, results.Z(n,:),'-','linewidth',2,'color','black') ;
                grid on ; hold on
                q               = plot(results.Dates, results.Zphi(n,:),'-.','linewidth',2,'color','black') ;
                leg             = legend([p,q],'$\mathrm{Z}_t$','$\Phi_t^{\mathrm{Z}}$','location','best');
                leg.Interpreter = 'Latex';
                leg.Color       = 'none';
                leg.FontSize    = results.FontSize;
                VX              = [results.Dates(1) results.Dates(end)] ;
                xlim(VX)
                xticklabels({})
                title(Countries(n),'Interpreter','Latex')
                ylabel('infection counts','Interpreter','Latex')
                set(gca,'ticklabelinterpreter','Latex','fontsize',results.FontSize,'color','None')
                subplot(212)
                plot(results.Dates,ones(size(results.Z(n,:))),'k-','LineWidth',1);
                grid on; hold on
                Q               = [];
                iEst            = 1;
                if ~isempty(find(strcmp(Estimates,'MLE'),1))  && ~isempty(find(strcmp(results.Estimates,'MLE'),1))
                    q1          = plot(results.Dates, results.MLE(n,:),'-','linewidth',2,'color',[0.7,0.7,0.7]) ;
                    Q           = [Q, q1];
                    L{iEst}     = '$\mathrm{R}^{\mathrm{MLE}}_t$';
                    iEst        = iEst + 1;
                end
                if ~isempty(find(strcmp(Estimates,'Gamma'),1)) && ~isempty(find(strcmp(results.Estimates,'Gamma'),1))
                    q2          = plot(results.Dates, results.Gamma(n,:),'-','linewidth',2,'color',[0,0.5,0]) ;
                    Q           = [Q, q2];
                    L{iEst}     = '$\mathrm{R}_{\tau,t}^\Gamma$';
                    iEst        = iEst + 1;
                end
                if ~isempty(find(strcmp(Estimates,'U-12'),1)) && ~isempty(find(strcmp(results.Estimates,'U-12'),1))
                    q3          = plot(results.Dates(3:Teff), results.U_12(n,:),'-','linewidth',2,'color','blue') ;
                    Q           = [Q, q3];
                    L{iEst}     = '$\mathrm{R}^{\mathrm{U-12}}_t$';
                    iEst        = iEst + 1;
                end
                if ~isempty(find(strcmp(Estimates,'U-C-12'),1)) && ~isempty(find(strcmp(results.Estimates,'U-C-12'),1))
                    q4          = plot(results.Dates(3:Teff), results.U_C_12(n,:),'-','linewidth',2,'color','red') ;
                    Q           = [Q, q4];
                    L{iEst}     = '$\mathrm{R}^{\mathrm{U-C-12}}_t$';
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
