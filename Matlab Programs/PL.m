function PL(skip)
%xvals=xvals';
%deltat is time interval between scans
%wait is the time inbetween intervals I want (these need to divide to get
%an integer)
%lambda is a value for reference wavelength you want to hold intensity
%constant between samples (program will normalize to)


% Prompts user for a file, then it looks up all files in that folder
[filename, pathname]=uigetfile('*.txt', 'Open any file in folder you want to extract from');
Tfiles=dir(strcat(pathname, '*.txt'));

% list of all files in folder which are not directories
Tfiles=Tfiles([Tfiles.isdir]==0);

%initialize data holders
%headernames='Wavelength (nm),';
t = 0; %initialize time

% Get data
%subplot(1,2,1)
figure(1)
hold on
%legend('-DynamicLegend')
n = 1;
for j=1:floor(length(Tfiles)/skip)
    %Get Data
    CurrentFile=strcat(pathname,Tfiles(1+(j-1)*skip).name);
    NewData=dlmread(CurrentFile,'\t',14,0);
    wave = NewData(:,1);
    data = NewData(:,2); 
    PLdata(:,1) = wave;
    PLdata(:,j+1) = data; %save raw PL data from spectrometer in huge array
    
    %Determine Integration time
    %fid = fopen(CurrentFile); %opens up file as text
    %temp = textscan(fid,'%s'); %scans text file
    %integ = str2num(temp{1}{22}); %pulls integration time in seconds from file
    %fclose(fid)
    integ = .4;
    %Determine PL center of mass
    com(j) = sum(PLdata(:,1).*PLdata(:,j+1)/sum(PLdata(:,j+1)));
    time(j) = t + j*integ;
    n= n+1;
end

%remove data before I turn the light on
M = find(com == min(com(9:end)));
time2 = time - time(M);
plot(time2,com)
figure(1)
xlabel('Time [s]')
ylabel('PL Center of Mass [nm]')
save('PLdata.mat', 'PLdata')

%make an additional plot of PL spectra at start and end of light exposure
figure(2)
hold on
plot(PLdata(:,1),PLdata(:,M)*10,PLdata(:,1),PLdata(:,end))