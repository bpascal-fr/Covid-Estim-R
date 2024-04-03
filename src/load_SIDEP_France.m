% Load Covid19 new infections counts in the French departments, including 
% overseas territories, during ~ 3 years of Covid19 pandemic made available
% by the French government on the repository: 
% https://www.data.gouv.fr/fr/datasets/donnees-de-laboratoires-pour-le-depistage-a-compter-du-18-05-2022-si-dep/.
% Cases are collected by Regional Health Agencies.
%
% Last release of the SIDEP repository on June 30, 2023.
%
% On March 31, 2024, data are still available at https://coronavirus.jhu.edu/
% the present code proposes to either:
% - download them from the Internet,
% - load a .csv file stored in the folder data.
%
% See data/README.md for more detail.
%
% B. Pascal and P. Abry, March 2024.



function [Z_User, Zphi_User, Dates_User, Departments_User] = load_SIDEP_France(Departments,opts)

    % Inputs:  - Departments: list of French departments referred to by their INSEE number that are to be displayed 
    %                         (default: if Departments = [] all monitored departments)
    %          - opts: structure indicated the time period selected by the user containing (optional)
    %                   opts.LastDay: last day of the time period in format 'YYYY-MM-DD' (default '2023-06-27', latest day possible)
    %                   opts.T: number of day of the time period (default: total number of days available in SIDEP repository).
    %                   opts.Download: 0 for loading the .csv in folder data, 1 for downloading data from https://coronavirus.jhu.edu/
    %
    %
    % Outputs: - Z_User: new infection counts in the D French departments and during T days, stored in a D x T matrix.
    %          - Zphi_User: global infectiousness, defined as a weighted sum of past counts, in the D French departments and during T days, stored in a D x T matrix.
    %          - Dates_User: T dates of the time period in datetime format
    %          - Departments_User: structure containing information about the selected departments
    %                    Names: list of French departments names
    %                    Number: official INSEE numerotation of French departments
    %                    Indices: index of each department in the storage from 1 to D


    if nargin < 1

        % Last day of the time period
        ChosenTime = [2023,06,27];

        % Length of the studied period
        TUser      = -1;

        % Load data from local folder
        Download   = 0;

    else
        
        if ~isfield(opts,'LastDay'),  opts.LastDay  = '2023-06-27'; end
        if ~isfield(opts,'T'),        opts.T        = -1; end
        if ~isfield(opts,'Download'), opts.Download = 0; end
        
        % Last day of the period
        LastDay    = opts.LastDay;
        ChosenTime = [str2double(LastDay(1:4)),str2double(LastDay(6:7)),str2double(LastDay(9:10))];

        % Length of the studied period
        TUser      = opts.T;

        % Load or download data
        Download   = opts.Download;

    end
   
    % List of all French departments, including overseas territories
    AllDepartments      = get_FrenchDepartments(0); 

    % By default select all departments
    if isempty(Departments)
        
        Departments = AllDepartments.Numbers;

    end

    %% LOAD DATA MADE AVAILABLE BY THE FRENCH GOVERNMENT

    % Download or load from stored data
    if Download
        url         = 'https://www.data.gouv.fr/fr/datasets/r/426bab53-e3f5-4c6a-9d54-dba4442b3dbc';
        filename    = 'data/Covid_French_Department_downloaded.csv';
        outfilename = websave(filename,url);
    else
        outfilename = 'data/Covid_French_Department_stored.csv';
    end
    
    % Set up the Import Options 
    import                  = delimitedTextImportOptions("NumVariables", 9);
    import.DataLines        = [2, Inf]; % Specify range 
    import.Delimiter        = ";"; % Specify delimiter
    import.VariableNames    = ["dep", "jour", "pop", "P", "T", "Ti", "Tp", "Td", "cl_age90"]; % Specify column names 
    import.VariableTypes    = ["string", "datetime", "double", "double", "double", "double", "double", "double", "double"]; % Specify column types
    import.ExtraColumnsRule = "ignore"; % Specify file level properties
    import.EmptyLineRule    = "read";
    import                  = setvaropts(import, "jour", "InputFormat", "yyyy-MM-dd"); % Specify variable properties
    import                  = setvaropts(import, ["pop", "P", "T", "Ti", "Tp", "Td", "cl_age90"], "DecimalSeparator", ",");
    import                  = setvaropts(import, ["pop", "P", "T", "Ti", "Tp", "Td", "cl_age90"], "ThousandsSeparator", ".");

    % Load Covid19 confirmed cases in all departments
    SIDEP                   = readtable(outfilename, import);

    % Store data in correct format for processing
    Zdata_repeat            = SIDEP.P;
    tdate_repeat            = SIDEP.jour;
    tdate                   = unique(tdate_repeat); 
    T_All                   = length(tdate); % Total number of days 

    Departments_repeat      = SIDEP.dep;
 
    % Select departments chosen for analysis
    D_User                   = length(Departments); % number of selected departments
    N_User                   = zeros(D_User,1); % indices of the selected departments
    Zdata                    = zeros(D_User,T_All) ; 
    for d = 1:D_User
        index                = find(strcmp(Departments_repeat,Departments(d)));
        N_User(d)            = find(strcmp(AllDepartments.Numbers,Departments(d)));
        if length(index)     == length(tdate)
            Zdata(d,:)       = Zdata_repeat(index)' ;
        else
            datestmp         = tdate_repeat(index) ;
            [~,~,index]      = intersect(datestmp,tdate) ;
            Zdata(d,index)   = Zdata_repeat(index)' ;
        end
    end
    Departments_User.Names_Table = AllDepartments.Names_Table(N_User);
    Departments_User.Names       = AllDepartments.Names(N_User);
    Departments_User.Numbers     = AllDepartments.Numbers(N_User);
    Departments_User.Indices     = AllDepartments.Indices(N_User);


    %% PREPROCESSING

    % Remove negative counts
    Zdata(Zdata < 0) = 0 ;

    % Apply a sliding median to discard outlier values
    alpha            = Inf;
    Zalpha           = sliding_median(Zdata,alpha,7);

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
    % one day cropped because infectiousness not defined at day 1 due to border effect induced by absence of any past count
    


    %% FOCUS ON THE SPECIFIED TIME PERIOD

    % Find the TUser dates
    tChosenTime = datenum(ChosenTime);
    index       = find(datenum(tdate) == tChosenTime);

    % Define last day
    if isempty(index)
        warning('Required last day not covered by SIDEP repository, replaced by the last availble day: June 27, 2023.')       
        index = size(Zalpha,2);
    end

    % Define the length of the time period
    if TUser == -1
        date = 1:index;
    else
        if index-TUser+1 <= 0
            warning(strcat('Time period not fully included in SIDEP repository hence adjusted to T = '," ",num2str(index),' days.'))
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
    monitored = Departments_User.Names_Table(1); 
    for d     = 2:length(Departments_User.Names_Table)
        monitored = strcat(monitored,", ",Departments_User.Names_Table(d));
    end

    disp('---------------------------------------------------------')
    disp('DEPARTMENTS AND TIME PERIOD CHOSEN FOR ANALYSIS:')
    disp(strcat("Monitored French departments: ",monitored,"."))
    disp(strcat("From ",string(Dates_User(1))," to ",string(Dates_User(end))))
    disp('---------------------------------------------------------')

end