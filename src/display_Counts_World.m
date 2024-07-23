% Display the new infection counts for all countries selected and during 
% the analyzed time period.
%
% Implementation B. Pascal,
% March, 2024


function display_Counts_World(Countries,results)


    % Inputs:  - Countries: list of countries that are to be displayed 
    %                       (default: if Countries = [] all monitored countries)
    %          - results: input data 
    %                     Countries: list of the C countries monitored
    %                     Dates: the T dates corresponding to the considered time period
    %                     Z: infection count times series stored as a matrix of size C x T
    %                     Zphi: global infectiousness times series stored as a matrix of size C x T
    %                     FontSize: desired FontSize for the plots
    
    if isempty(Countries),           Countries        = results.Countries; end
    if ~isfield(results,'Dates'),    results.Dates    = 1:size(Z,2); end
    if ~isfield(results,'FontSize'), results.FontSize = 22.5; end
    

    for n = 1:length(Countries)

        if ~sum(contains(results.Countries,Countries(n)))
            warning(strcat(Countries(n)," was not analysed, either because it was not found in JHU repository or because it has not been indicated by the user in the preamble. Hence, ",Countries(n)," will not appear in the plots."))
        else
            ine             = find(strcmp(results.Countries,Countries(n)),1);
            f1              = figure(1000 + n); clf; 
            p               = plot(results.Dates, results.Z(ine,:),'-','linewidth',2,'color','black') ;
            grid on ; hold on
            q               = plot(results.Dates, results.Zphi(ine,:),'-.','linewidth',2,'color','black') ;
            leg             = legend([p,q],'$\mathrm{Z}_t$','$\Phi_t^{\mathrm{Z}}$','location','best');
            leg.Interpreter = 'Latex';
            leg.Color       = 'none';
            leg.FontSize    = results.FontSize;
            VX              = [results.Dates(1) results.Dates(end)] ;
            xlim(VX)
            title(Countries(n),'Interpreter','Latex')
            ylabel('infection counts','Interpreter','Latex')
            set(gca,'ticklabelinterpreter','Latex','fontsize',results.FontSize,'color','None')
            f1.Position     = [211 287 943 314];
        end

    end

end
