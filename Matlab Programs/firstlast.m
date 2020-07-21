% To be run after oceanoptics.m This function opens an exported .txt plots
% the first and last files within that table

%function [i1, iend] = firstlast
[filename, pathname]=uigetfile('*.txt', 'Open any PL file in folder you want to extract from');
Tfiles=dir(strcat(pathname, '*.txt'));

% list of all files in folder which are not directories
Tfiles=Tfiles([Tfiles.isdir]==0);
figure
hold on

for j=1:length(Tfiles)
    
    CurrentFile=convertCharsToStrings(strcat(pathname,Tfiles(j).name));
    % Specify column names and types
    samplename = convertCharsToStrings(Tfiles(j).name);
    % Import the data
    temp = table2array(readtable(CurrentFile));
    first(:,1) = temp(:,1);
    first(:,j+1) = temp(:,2);
    last(:,1) = temp(:,1);
    last(:,j+1) = temp(:,end);
    plot(first(:,1),first(:,j+1),last(:,1),last(:,j+1))
end 
%end