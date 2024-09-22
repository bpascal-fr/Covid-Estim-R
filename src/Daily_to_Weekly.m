% Turn multivariate daily infection counts into weekly infection 
% counts by aggregating counts over seven days in each territory separately
% and compute the associated  weekly global infectiousness.
% 
% The daily discretized serial interval function is also weekly discretized
% using the left-rectangle integration method with constant step size of 
% one day length.
%
% References
%
% - Nash, R. K., Bhatt, S., Cori, A., & Nouvellet, P. (2023). Estimating the 
% epidemic reproduction number from temporally aggregated incidence data: 
% A statistical modelling approach and software tool. 
% PLOS Computational Biology, 19(8), e1011439.
%
% - Pascal, B., Vaiter, S. (2024, September). Risk Estimate under a
% Nonstationary Autoregressive Model for Data-Driven Reproduction Number
% Estimation. Preprint. arXiv:.
%
% B. Pascal and P. Abry, May 2024.


function [Z_Week, Zphi_Week, M_Week] = Daily_to_Weekly(Z_Day, opts)

    % Inputs: - Z_Day: multivariate daily infection count time series 
    %         - opts: parameters of the model, structure containing
    %                   opts.FontSize: font size in the plots (default FontSize = 22.5)
    %                   opts.Dates: abstract dates in datetime format for display
    %                   opts.Phi: daily discretized serial interval function to be used (default: daily discretized Covid19 serial interval function)
    %                   opts.R: ground truth daily reproduction number
    %                   opts.Countries: list of the C countries monitored
    %
    %
    % Outputs: - Z_Week: weekly aggregated infection counts
    %          - Zphi_Week: associated global infectiousness
    %          - M_Week: weekly model parameters, structure containing
    %                   M_Week.Phi: weekly discretized serial interval function
    %                   M_Week.Phi_Day: daily discretized serial interval function
    %                   M_Week.R: weekly reproduction number (If daily reproduction number is provided.)
    %                   M_Week.Dates: dates corresponding to the last days of the weeks over which infection counts are aggreagted.  (If dates are provided provided.)

    %% RESIZE INPUT 

    [d1,d2]     = size(Z_Day);

    if min(d1,d2) == 1
        
        Z_Day   = reshape(Z_Day,1,max(d1,d2));

        if isfield(opts,'R')
            [e1,e2]     = size(opts.R);
            if (e1 == d1)&(e2 == d2)
                opts.R  = reshape(opts.R,1,max(d1,d2));
            else
                warning('The reproduction number should have the same has Z. Since it is not the case it will be ignored.')
            end
        end

    end

    if nargin < 2

        % plot font size
        FontSize = 22.5;

        % temporal axis
        Dates     = 1:size(Z_Day,2);

    else

        if ~isfield(opts,'FontSize'); opts.FontSize = 22.5;        end
        if ~isfield(opts,'Dates');    opts.Dates = 1:length(GT.R); end

        % plot font size
        FontSize = opts.FontSize;

        % temporal axis
        Dates    = opts.Dates;

    end

    % Aggregation of daily infection counts over seven days to get weekly counts
    W            = floor(length(Z_Day)/7);
    Z_Week       = zeros(size(Z_Day,1),W);
    for w = 1:W
        tmp_Z         = Z_Day(:,7*(w-1)+1:7*w);
        Z_Week(:,w)   = sum(tmp_Z,2);
    end

    % Corresponding dates
    Dates_Day    = Dates(1:7*W); 
    Dates_Week   = Dates_Day(7:7:end);

    % Construct the discretized serial interval function
    if ~isfield(opts,'Phi') % default: Covid19 serial interval function
        tau_day  = 25;                            % memory horizon
        shape    = 1/0.28;                        % shape parameter
        scale    = 1.87;                          % scale parameter
        Phi      = gampdf(0:tau_day,shape,scale); % discretized Gamma probability density function
    else
        % the serial interval function should start with a zero and should be a row vector
        if ~(opts.Phi(1) == 0)
            tau_day = length(opts.Phi);
            Phi     = reshape(opts.Phi,1,tau_day);
            Phi     = [0, Phi];
        else
            tau_day = length(opts.Phi) - 1;
            Phi     = reshape(opts.Phi,1,tau_day+1);
        end
    end
    % tau - 1 should be a multiple of seven
    if mod(tau_day,7)
        add_day = 7 - mod(tau_day,7);
        tau_day = tau_day + add_day;
        Phi     = [Phi, zeros(1,add_day)];
    end
    Phi     = Phi/sum(Phi);             % normalize the serial interval function (optional)

    % Aggregation over seven days to get a weekly discretized interval function
    Phi_Week = zeros(1,tau_day/7);
    for w = 1:tau_day/7
        tmp_Phi       = Phi(7*(w-1)+2:7*w+1);
        Phi_Week(w+1) = sum(tmp_Phi);
    end

    % Weekly global infectiousness
    [Zphi_Week,Z_Week] = Phi_normal(Z_Week,Phi_Week);

    % One week cropped because infectiousness not defined at week 1 due to border effect induced by absence of any past count
    Z_Day              = Z_Day(:,8:end);
    Dates              = Dates(8:end);
    Dates_Week         = Dates_Week(2:end);
   

    %% DISPLAY DAILY VS. WEEKLY INFECTION COUNTS

    for n = 1:size(Z_Week,1)
        fn                = figure(100 + n); clf
        q                 = plot(Dates_Week,Z_Week(n,:),'-^','linewidth',1.5,'markersize',10,'color','black');
        grid on ; hold on
        ylabel('weekly counts','Interpreter','Latex')
        yyaxis right
        p                 = plot(Dates,Z_Day(n,:),'linewidth',1.5,'color',[0.5,0.5,0.5]);
        ax                = gca;
        ax.YAxis(2).Color = [0.5,0.5,0.5];
        ax.YAxis(1).Color = 'black';
        xlim([Dates(1) Dates(end)])
        leg             = legend([p,q],'$\mathrm{Z}_t^{\mathrm{Day}}$','$\mathrm{Z}_w^{\mathrm{Week}}$','location','best');
        leg.Interpreter = 'Latex';
        leg.Color       = 'none';
        leg.FontSize    = FontSize;
        ylabel('daily counts','Interpreter','Latex')
        if isfield(opts,'Countries')
            title(strcat("Fom daily to weekly counts in ",opts.Countries(n)),'Interpreter','Latex')
        else
            title(strcat("Fom daily to weekly counts in territory ",num2str(n)),'Interpreter','Latex')
        end
        set(gca,'ticklabelinterpreter','Latex','fontsize',FontSize,'color','None')
        fn.Position     = [211 287 943 314];
    end

     %% STORE RESULTS

    if d2 == 1
        Z_Week     = reshape(Z_Week,W,d2);
        Zphi_Week  = reshape(Zphi_Week,W,d2);
    end

    M_Week.Phi_Day = Phi(1:end-add_day);
    M_Week.Phi     = Phi_Week;
    M_Week.Dates   = Dates_Week;

end