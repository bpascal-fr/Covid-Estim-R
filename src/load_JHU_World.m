% Load Covid19 new infections counts in 200+ countries worldwide during
% ~ 3 years of Covid19 pandemic made available by the Johns Hopkins
% University on the repository: https://coronavirus.jhu.edu/ from confirmed
% cases collected by National Health Agencies.
%
% Last release of JHU repository on October 3, 2023.
%
% On March 31, 2024, data are still available at https://coronavirus.jhu.edu/
% the present code proposes to either:
% - download them from the Internet,
% - load a .csv file stored in the folder data.
%
% See data/README.md for more detail.
%
% B. Pascal and P. Abry, March 2024.



function [Z_User, Zphi_User, Dates_User, Countries] = load_JHU_World(User_Countries,opts)

    % Inputs:  - User_Countries: vector of strings containing the names of the countries which are to be monitored, 
    %                 e.g., User_Countries = ["France","United Kingdom","Spain"].
    %            If User_Countries = [] all 201 available countries are selected.
    %          - opts: structure indicating the time period selected by the user containing (optional)
    %                   opts.LastDay: last day of the time period in format 'YYYY-MM-DD' (default '2023-03-09', latest day possible)
    %                   opts.T: number of day of the time period (default: total number of days available in JHU repository).
    %                   opts.Download: 0 for loading the .csv in folder data, 1 for downloading data from https://coronavirus.jhu.edu/
    %                   opts.Outlier: if true apply a sliding median smoothing to discard outlier values (default: true)
    %                   opts.win: window size for the sliding median preprocessing (default value: 7)
    %
    %
    % Outputs: - Z_User: new infection counts in the C countries in Countries and during T days, stored in a C x T matrix.
    %          - Zphi_User: global infectiousness, defined as a weighted sum of past counts, in the C countries in Countries and during T days, stored in a C x T matrix.
    %          - Dates_User: T dates of the time period in datetime format
    %          - Countries: list of monitored countries; if some countries required by user are not found in JHU repository they are ignored.


    if nargin < 2

        % Last day of the time period
        ChosenTime = [2023,03,09];

        % Length of the studied period
        TUser      = -1;

        % Load data from local folder
        Download   = 0;

        % sliding median smoothing
        Outlier    = false;

        % Window size for the sliding median
        win        = 7;

    else
        
        if ~isfield(opts,'LastDay'),  opts.LastDay  = '2023-03-09'; end
        if ~isfield(opts,'T'),        opts.T        = -1; end
        if ~isfield(opts,'Download'), opts.Download = 0; end
        if ~isfield(opts,'Outlier');  opts.Outlier  = true; end
        if ~isfield(opts,'win'),      opts.win    = 7; end
        
        % Last day of the period
        LastDay    = opts.LastDay;
        ChosenTime = [str2double(LastDay(1:4)),str2double(LastDay(6:7)),str2double(LastDay(9:10))];

        % Length of the studied period
        TUser      = opts.T;

        % Load or download data
        Download   = opts.Download;

        % sliding median smoothing
        Outlier  = opts.Outlier;

        % Window size for the sliding median
        win        = opts.win;

    end
   

    %% LOAD DATA MADE AVAILABLE BY JOHNS HOPKINS UNIVERSITY

    % Download or load from stored data
    if Download
        url         = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv';
        filename    ='COVID-19-JHU_downloaded.csv';
        outfilename = websave(filename,url);
    else
        outfilename = 'data/COVID-19-JHU_stored.csv';
    end


    % Handle different encoding of the dates depending on Matlab version
    [~,d_MAT]      = version;      % Matlab release date
    y_MAT          = d_MAT(end-3:end); % year of release


    % Load Covid19 confirmed cases in all countries
    if strcmp(y_MAT, 2019)
        T       = readtable(outfilename);
        Dates   = datenum(table2array(T(1,5:end)));
        T       = readtable(outfilename,'ReadVariableNames',true);
    else
        opts                    = detectImportOptions(outfilename);
        opts.VariableNamesLine  = 1;
        opts.VariableNamingRule = 'preserve';
        T                       = readtable(outfilename,opts);
        datesAscii              = T.Properties.VariableNames ;
        Dates                   = datenum(datesAscii(5:end));
    end


    % Store data in correct format for processing
    CountriesTMP   = table2array(T(:,2)) ;
    Provinces      = table2array(T(:,1)) ;
    TConfirmedTMP  = table2array(T(:,5:end));
    Compteur       = 1 ;
    k              = 2 ;

    while k <= length(CountriesTMP)
        All_Countries{Compteur}            = CountriesTMP{k} ;
        if strcmp(Provinces{k},'')
            TConfirmed(Compteur,:)     = TConfirmedTMP(k,:) ;
            k                          = k+1 ;
        else
            TConfirmed(Compteur,:)     = TConfirmedTMP(k,:) ;
            k2                         = 1 ;
            while strcmp(CountriesTMP{k+k2},CountriesTMP{k})
                TConfirmed(Compteur,:) = TConfirmed(Compteur,:) + TConfirmedTMP(k+k2,:) ;
                k2                     = k2+1 ;
            end
            k                          = k+k2 ;
        end
        Compteur                       = Compteur + 1 ;
    end

    %% STORE INDICES CORRESPONDING TO COUNTRIES REQUIRED BY THE USER
    % If the country is not found in the list or misspelled it is ignored and
    % not included in the list of analyzed countries

    N_Countries = [];
    Countries   = [];

    for n = 1:length(User_Countries)

        N_n            = find(strcmp(All_Countries,User_Countries(n)));

        if isempty(N_n)

            warning(strcat(User_Countries(n),' not found, hence removed from Countries and ignored in the extraction.'))

        else

            N_Countries = [N_Countries, N_n];
            Countries   = [Countries, User_Countries(n)];

        end


    end


    %% EXTRACT NEW INFECTION COUNT MULTIVARIATE TIME SERIES

    % Daily cumulated number of cases
    Confirmed = TConfirmed(N_Countries,:);

    % Daily new infection counts
    Zdata     = diff(Confirmed,1,2);
    Dates     = Dates(2:end);

    % Store dates
    aa        = datevec(Dates(end)) ; aa = aa(1:3) ;
    bb        = datevec(Dates(1))   ; bb = bb(1:3) ;
    tdate     = datetime(bb(1),bb(2),bb(3)):datetime(aa(1),aa(2),aa(3)) ;

    % Start at the beginning of the epidemic for the studied countries
    begin_countries = zeros(size(Countries));

    for n = 1:length(Countries)

        ind_begin          = find(Zdata(n,:) > 0);
        begin_countries(n) = ind_begin(1);

    end

    begin_pandemic  = min(begin_countries);

    Zdata           = Zdata(:,begin_pandemic:end);
    tdate           = tdate(begin_pandemic:end);
    Dates           = Dates(begin_pandemic:end);

    %% PREPROCESSING

    % Remove negative counts
    Zdata(Zdata < 0) = 0 ;

    % Apply a sliding median to discard outlier values
    alpha            = 7 ;
    if Outlier
        Zalpha       = sliding_median(Zdata,alpha,win);
    else
        Zalpha       = Zdata;
    end

    %% COMPUTE GLOBAL INFECTIOUSNESS

    % Serial interval function modeled by a Gamma distribution with
    % - mean: 6.6 days
    % - standard deviation 3.5 days
    % truncated at 25 days and normalized.
    shape    = 1/0.28;
    scale    = 1.87;
    tau_phi  = 25;
    Phi      = gampdf(0:tau_phi,shape,scale);
    Phi      = Phi/sum(Phi); % normalize the weights applied to past tau_phi infection counts

    % Infectiousness: weighted sum of past infection counts
    [Zphi,Zalpha] = Phi_normal(Zalpha,Phi);
    tdate         = tdate(2:end);
    Dates         = Dates(2:end);
    % one day cropped because infectiousness not defined at day 1 due to border effect induced by absence of any past count
    


    %% FOCUS ON THE SPECIFIED TIME PERIOD

    % Find the TUser dates
    tChosenTime = datenum(ChosenTime);
    index       = find(Dates == tChosenTime);

    % Define last day
    if isempty(index)
        warning('Required last day not covered by JHU repository, replaced by the last availble day: March 9, 2023.')       
        index = size(Zalpha,2);
    end

    % Define the length of the time period
    if TUser == -1
        date = 1:index;
    else
        if index-TUser+1 <= 0
            warning(strcat('Time period not fully included in JHU repository hence adjusted to T = '," ",num2str(index),' days.'))
            date        = 1:index;
        else
            date        = index-TUser+1:index;
        end
    end

    %% INFECTION COUNTS, INFECTIOUSNESS AND DATES

    % Extract infection counts, infectiousness and dates for the specified time period
    Z_User          = Zalpha(:,date);
    Zphi_User       = Zphi(:,date);
    Dates_User      = tdate(date);


    % Messages to the user
    monitored = User_Countries(1); 
    for n     = 2:length(User_Countries)
        monitored = strcat(monitored,", ",User_Countries(n));
    end

    disp('---------------------------------------------------------')
    disp('COUNTRIES AND TIME PERIOD CHOSEN FOR ANALYSIS:')
    disp(strcat("Monitored countries: ",monitored,"."))
    disp(strcat("Time period: from ",string(Dates_User(1))," to ",string(Dates_User(end))))
    disp('---------------------------------------------------------')

end