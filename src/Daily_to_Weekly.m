% Turn daily infection counts Z(t) into weekly infection counts by
% aggregating Z(t) over seven days and compute the associated weekly global
% infectiousness.
% 
% The daily discretized serial interval function is also weekly discretized
% using the left-rectangle integration method with constant step size of 
% one day length.
%
% - XXX
%
% B. Pascal and P. Abry, May 2024.


function [Z_Week, Zphi_Week, M_Week] = Daily_to_Weekly(Z_Day, opts)

    % Inputs: - Z_Day: daily infection count time series
    %         - opts: parameters of the model, structure containing
    %                   opts.FontSize: font size in the plots (default FontSize = 22.5)
    %                   opts.Dates: abstract dates in datetime format for display
    %                   opts.Phi: daily discretized serial interval function to be used (default: weekly discretized Covid19 serial interval function)
    %                   opts.R: ground truth daily reproduction number
    %
    %
    % Outputs: - Z_Week: weekly aggregated infection counts
    %          - Zphi_Week: associated global infectiousness
    %          - M_Week: weekly model parameters, structure containing
    %                   M_Week.Phi: weekly discretized serial interval function
    %                   M_Week.R: ground truth weekly reproduction number (If ground truth daily reproduction number is provided.)
    %                   M_Week.Dates: dates corresponding to the last day of the week over which infection counts are aggreagted.  (If dates are provided provided.)

    if nargin < 2

        % plot font size
        FontSize = 22.5;

        % temporal axis
        Dates     = 1:length(Z_Day);

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
    Z_Week       = zeros(1,W);
    for w = 1:W
        tmp_Z         = Z_Day(7*(w-1)+1:7*w);
        Z_Week(w)     = sum(tmp_Z);
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
    for w = 1:4
        tmp_Phi       = Phi(7*(w-1)+2:7*w+1);
        Phi_Week(w+1) = sum(tmp_Phi);
    end

    % Weekly global infectiousness
    Zphi_Week  = conv(Z_Week,Phi_Week);
    Zphi_Week  = Zphi_Week(1:length(Z_Week));

    %% STORE RESULTS

    M_Week.Phi   = Phi_Week;
    M_Week.Dates = Dates_Week;

    %% DISPLAY DAILY VS. WEEKLY INFECTION COUNTS
    
    f4                = figure(4); clf
    q                 = plot(Dates_Week,Z_Week,'-^','linewidth',1.5,'markersize',10,'color','black');
    grid on ; hold on
    % ylim([0 1.1*max(Z_Day)])
    ylabel('weekly counts','Interpreter','Latex')
    yyaxis right
    p                 = plot(Dates,Z_Day,'linewidth',1.5,'color',[0.5,0.5,0.5]);
    ax                = gca;
    ax.YAxis(2).Color = [0.5,0.5,0.5];
    ax.YAxis(1).Color = 'black';
    xlim([Dates(1) Dates(end)])
    leg             = legend([p,q],'$\mathrm{Z}_t^{\mathrm{Day}}$','$\mathrm{Z}_w^{\mathrm{Week}}$','location','best');
    leg.Interpreter = 'Latex';
    leg.Color       = 'none';
    leg.FontSize    = FontSize;
    ylabel('daily counts','Interpreter','Latex')
    title('Fom daily to weekly infection counts','Interpreter','Latex')
    set(gca,'ticklabelinterpreter','Latex','fontsize',FontSize,'color','None')
    f4.Position     = [211 287 943 314];

end