%Jill
% XRD analysis file
%get files
% Prompts user for a file, then it looks up all files in that folder
[filename, pathname]=uigetfile('*.txt', 'Open any XRD file in folder you want to extract from');
Tfiles=dir(strcat(pathname, '*.txt'));

% list of all files in folder which are not directories
Tfiles=Tfiles([Tfiles.isdir]==0);
%graph start
hold on 

for count=1:length(files)
    CurrentFile=strcat(pathname,Tfiles(count).name);
    data=dlmread(CurrentFile);
    angle = data(:,1);
    intensity = data(:,2);
    t = t + deltat
    plot(angle,intensity, 'DisplayName', txt, 'color', [(1-count*.01),0,count*.01]) 
end

save('Tdata.mat')


