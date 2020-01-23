
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



