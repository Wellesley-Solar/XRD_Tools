
%Made by Becky to analyze UV-Vis spectra taken on Chemistry Integrating
%sphere

fileNames = {};
pathNames = {};
data = {};
[addfileNames, addpathNames] = uigetfile('*.csv*','MultiSelect', 'on');
fileNames = [fileNames; addfileNames];
pathNames = [pathNames; addpathNames];
 
figure 
hold on
xlabel('Wavelength [nm]')
ylabel('Absorbance [OD]')

for i = 1:length(fileNames)
    combinedName = strcat(pathNames{1},fileNames{i});
    fid = fopen(combinedName,'r')
    C = textscan(fid, '%f %f', 'Delimiter', ',', 'headerLines', 2)
    data{i}.N = fileNames{i}
    data{i}.W = C{1};
    data{i}.A = C{2};
    plot(data{i}.W,data{i}.A,'DisplayName',data{i}.N);
end

figure
hold on 

%remove background
for i = 1:(length(fileNames)-1)
    newdata{i}.A = data{i}.A-data{2}.A;
    plot(data{i}.W,newdata{i}.A,'DisplayName',data{i}.N);
end


% %Trying to fituvvis data from Cary
% %need to integrate into absorption_chemistry code
% 
% figure
% hold on
% 
% lowernm = 250;
% uppernm = 550;
% lowerback = 600;
% upperback = 800;
% 
% uv_switch = 349;
% 
% [d, lim1] = min(abs(ABSData(:,1)-lowernm));
% [d, lim2] = min(abs(ABSData(:,1)-uppernm));
% [d, lim3] = min(abs(ABSData(:,1)-lowerback));
% [d, lim4] = min(abs(ABSData(:,1)-upperback));
% 
% wave = ABSData(lim2:lim1,1);
% 
% %correct for bump where the uv-light turns on
% [d, lim5] = min(abs(ABSData(:,1)-uv_switch));
% 
% for i = 1:8
%     shift(i) = ABSData(lim5-1,(i+2))-ABSData(lim5,i+2);
%     data_shift(:,i) = [ABSData(1:lim5-1,i+2); (ABSData(lim5:end,i+2)+shift(i))]; 
%     back(i) = mean(data_shift(lim4:lim3,i));
%     data_sub(:,i) = data_shift(lim2:lim1,i)-back(i);
%     plot(wave,data_sub(:,i))
% end
