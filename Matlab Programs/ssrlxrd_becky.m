%This program takes an imported xrd data set (stacked csv) and attempts to quassian fit up to three peaks.
%assumes the first column of a csv is diffraction angle in 2_theta

%used constants
ang = 0.982381; %wavelength of x-ray in angstroms

%open type of file in folder - will open all files of that type
fileNames = {};
pathNames = {};
data = {};
[addfileNames, addpathNames] = uigetfile('*.csv*','MultiSelect', 'on');
fileNames = [fileNames; addfileNames];
pathNames = [pathNames; addpathNames];

%convert file and path to full file name and import data into numeric array

for i = 1:length(fileNames);
    combinedName{i} = strcat(pathNames{i},fileNames{i});
    input = readtable(combinedName{i});
    headers = input.Properties.VariableNames;
    data{i} = table2array(readtable(combinedName{i}));   
end

%remove background and plot

%change from 2_that to q
data{1}(:,1) = (4*pi*sin(pi/180*data{1}(:,1)/2)/ang); 

%chose domain for plotting
lowerq = 0.5;
upperq = 4;

[d, lim1] = min(abs(data{1}(:,1)-lowerq));
[d, lim2] = min(abs(data{1}(:,1)-upperq));

data_sub = data{1}(lim1:lim2,:);

%plot file for subtracting background and manually chose points for
%background line
plot(data_sub(:,1),data_sub(:,2))
xlabel('Q [�^-^1]');
ylabel('Intensity [a.u.]');
[x,y] = getpts;
back = interp1(x,y,data_sub(:,1),'spline');

%remove background
it = size(data_sub);

for i = 1:(it(2)-1)
    data_sub(:,i+1) = data_sub(:,i+1)-back;
end

%save edited data
name = strcat(pathNames{1},'Edit',fileNames{1});
final = array2table(data_sub);
final.Properties.VariableNames = headers;
writetable(final,name)

%below is code to fit a narrow range of data to a guassian (for determining
%(100) peak positions
% 
% dt = 20 ;
% t(1) = 0;
% 
% for i = 1:30;
%     %i is length of xrd data set
%     f = fit(q(lim1:lim2),xrd(lim1:lim2,i+1)-back(lim1:lim2),'gauss2');
%     
%     %pull parameters for first peak
%     int1(i) = f.a1;
%     a1(i) = 2*pi./f.b1;
%     strain1(i) = f.c1;
%     
%     %pull parameters for second peak
%     int2(i) = f.a2;
%     a2(i) = 2*pi./f.b2;
%     strain2(i) = f.c2;
%    
%     %pull parameters for third peak
%     %int3(i) = f.a3;
%     %a3(i) = 2*pi./f.b3;
%     %strain3(i) = f.c3;
%     
%     t(i) = 0+dt*(i-1);
% end
% 
% figure
% yyaxis left
% plot(t,int1)
% yyaxis right
% plot (t,a1)
% 
% figure
% yyaxis left
% plot(t,int2)
% yyaxis right
% plot (t,a3)
%     
