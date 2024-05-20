clear all
close all
clc

addpath(genpath(pwd))

% All estimators process infection count times series of size C x T
% - C: number of territories monitored
% - T: number of days in the time period.
%
% Univariate estimators are computed in parallel for each territory.


%% LOAD AND DISPLAY NEW INFECTION COUNTS

% List all countries referenced in JHU repository
list             = 1;
get_AllCountries(list)
% list = 0: do not display countries
% list = 1: display all countries names

% List of countries to monitor
User_Countries             = ["France","United Kingdom","Spain"]; % If [] all 201 countries referenced in the JHU repository are selected.
% Countries are stored in the order indicated in User_Countries or by alphabetic order if User_Countries = [].

% Time period
opts_load.LastDay           = '2023-02-05'; % Last day of the time period in format 'YYYY-MM-DD'. If not provided set to March 9, 2023.
opts_load.T                 = 490;           % Length of the time period. If -1 or not provided entire available time period.

% Load data
[Z, Zphi, Dates, Countries] = load_JHU_World(User_Countries,opts_load);

% Store countries, dates, infection counts and infectiousness
results.Countries = Countries;
results.Dates     = Dates;
results.Z         = Z;
results.Zphi      = Zphi;

% Countries to be displayed
DisplayCountries  = ["France","United Kingdom"];

% Display infection counts time series for monitored countries
results.FontSize  = 22.5; % Adjust font size in plots
display_Counts_World(DisplayCountries,results) % If DisplayCountries = [] all monitored countries are displayed.

%% AGGREGATE DAILY INFECTION COUNTS TO GET WEEKLY INFECTION COUNTS

% Options for daily to weekly infection count aggregation
opts_week.Dates     = Dates;
opts_week.Countries = Countries;

% From daily to weekly infection counts
[Z_Week, Zphi_Week, M_Week] = Daily_to_Weekly(Z, opts_week);