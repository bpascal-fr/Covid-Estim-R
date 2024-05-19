% Display the new infection counts and estimated reproduction number for
% all French departments selected, possibly including overseas territories, 
% during the analyzed time period.
%
% Implementation B. Pascal,
% March, 2024


function display_Counts_France(Departments,results)


    % Inputs:  - Departments: list of French departments, referred to by their INSEE numbers, that are to be displayed 
    %                         (default: if Departments = [] all monitored departments)
    %          - results: input data 
    %                     Departments: list of the D French departments Names and associated Numbers and Indices
    %                     Dates: the T dates corresponding to the considered time period
    %                     Z: infection count times series stored as a matrix of size D x T
    %                     Zphi: global infectiousness times series stored as a matrix of size D x T
    %                     FontSize: desired FontSize for the plots
    
    if isempty(Departments),         Departments      = results.Departments.Numbers; end
    if ~isfield(results,'Dates'),    results.Dates    = 1:size(Z,2); end
    if ~isfield(results,'FontSize'), results.FontSize = 22.5; end

    for d = 1:length(Departments)
        
        if ~sum(strcmp(results.Departments.Numbers,Departments(d)))
            warning(strcat(Departments(d)," will be ignored in the plots, either because it is not a valid French department INSEE number or because it was not included in the analysis."))
        else
            ind             = find(strcmp(results.Departments.Numbers,Departments(d)));
            f1              = figure(1000 + d); clf
            p               = plot(results.Dates, results.Z(ind,:),'-','linewidth',2,'color','black') ;
            grid on ; hold on
            q               = plot(results.Dates, results.Zphi(ind,:),'-.','linewidth',2,'color','black') ;
            leg             = legend([p,q],'$\mathrm{Z}_t$','$\Phi_t^{\mathrm{Z}}$','location','best');
            leg.Interpreter = 'Latex';
            leg.Color       = 'none';
            leg.FontSize    = results.FontSize;
            VX              = [results.Dates(1) results.Dates(end)] ;
            xlim(VX)
            title(results.Departments.Names(ind),'Interpreter','Latex')
            ylabel('infection counts','Interpreter','Latex')
            set(gca,'ticklabelinterpreter','Latex','fontsize',results.FontSize,'color','None')
            f1.Position     = [211 287 943 314];
        end

    end

end
