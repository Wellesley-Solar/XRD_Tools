% This function will open data exported from ocean-optics spectrometers and
% save it as a stacked csv in the original folder with sample names as
% headers

function data = oceanoptics
[filename, pathname]=uigetfile('*.txt', 'Open any PL file in folder you want to extract from');
Tfiles=dir(strcat(pathname, '*.txt'));

% list of all files in folder which are not directories
Tfiles=Tfiles([Tfiles.isdir]==0);

for j=1:length(Tfiles)
    
    CurrentFile=convertCharsToStrings(strcat(pathname,Tfiles(j).name));
    
    % Specify range and delimiter
    opts = delimitedTextImportOptions("NumVariables", 2);
    opts.DataLines = [15, Inf];
    opts.Delimiter = "\t";

    % Specify column names and types
    samplename = convertCharsToStrings(Tfiles(j).name);
    opts.VariableNames = ["Wavelength", samplename];
    opts.VariableTypes = ["double", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Import the data
    temp = readtable(CurrentFile, opts);
    data(:,1) = temp(:,1);
    data(:,j+1) = temp(:,2);
    %To help me interpret when I've turned the light on (a.k.a. what is my
    %first frame)
    intensity(j) = sum(table2array(temp(:,2)));
    frame(j) = j+1;

%% Clear temporary variables
clear opts
end 

%cut dead frames at beginning and end
plot(frame,intensity);
slope = gradient(intensity(:)) ./ gradient(frame(:));
[a,b] = max(slope);
data_sub = [data(:,1) data(:,b:end-5)];

%export the data to the desired folder
pathname = uigetdir;
SavedName =convertCharsToStrings(strcat(pathname, Tfiles(j).name));
writetable(data_sub, SavedName);
end


% work on using outputted data in new function