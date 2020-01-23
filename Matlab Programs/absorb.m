function absorb(deltat)
%xvals=xvals'; 

% Prompts user for a file, then it looks up all files in that folder
[filename, pathname]=uigetfile('*.txt', 'Open any PL file in folder you want to extract from');
Tfiles=dir(strcat(pathname, '*.txt'));

% list of all files in folder which are not directories
Tfiles=Tfiles([Tfiles.isdir]==0);

%initialize data holders
%headernames='Wavelength (nm),';
t = 0 %initialize time

% Get data
%subplot(1,2,1)
hold on
%legend('-DynamicLegend')

for j=1:length(Tfiles)
    CurrentFile=strcat(pathname,Tfiles(j).name);
    NewData=dlmread(CurrentFile,'\t',14,0);
    wave = NewData(:,1);
    tdata = -log10(NewData(:,2)./100); %change %transmission to absorbance in OD
    Adata(:,1) = wave;
    Adata(:,j+1) = tdata; 
    t = t + deltat;
    txt = ['t=',num2str(t)];
    plot(wave,Adata(:,j+1), 'DisplayName', txt, 'color', [(1-j*.01),0,j*.01]) 
end

save('Tdata.mat', 'Adata')
