%Jill Mankoff 1/22/2020
%File for reading and plotting temperatures from txt files from SSRL
%use directory with all txt files to get a plot
myDir = uigetdir;
files= dir(fullfile(myDir,'*.txt'));
temps = []
for k=1:length(files)
    fname = files(k).name
    fullfname = fullfile(myDir, fname);
    txt = importdata(fullfname);
    %find temp
    tempc = txt{4,1}
    tempstr = tempc(1:end-1)
    temp = str2double(tempstr)
    temps = [temps temp]
end

figure
plot(temps)
