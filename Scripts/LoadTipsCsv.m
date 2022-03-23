function tipsData = LoadTipsCsv(file)%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: C:\home\UMICH\ascites\data\post_tips_us_14d_evb.csv
%
% Auto-generated by MATLAB on 21-Mar-2022 21:35:13

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 15);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "StudyID", "pre_tte_date", "tips_date", "Stenttypeused", "StentInternalDiameter", "GraftLinedLength", "GraftUnlinedLength", "PRA_base", "PHV_base", "PPV_base", "PHV_base", "PRA_TIPS", "PHV_TIPS", "PPV_TIPS"];
opts.VariableTypes = ["double", "double", "datetime", "datetime", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Stenttypeused", "EmptyFieldRule", "auto");
opts = setvaropts(opts, "pre_tte_date", "InputFormat", "yyyy-MM-dd");
opts = setvaropts(opts, "tips_date", "InputFormat", "yyyy-MM-dd");
opts = setvaropts(opts, ["PHV_base", "PHV_TIPS"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["PHV_base", "PHV_TIPS"], "ThousandsSeparator", ",");

% Import the data
tipsData = readtable(file, opts);


%% Clear temporary variables
clear opts