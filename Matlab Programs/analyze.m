
dat = importdata('D1_MAPbIBr2_Xraydeg.csv');    %Imports csv
new_dat=sub_bg(dat);    %Process to remove background
plot_ontop(new_dat, 0, 5);    %Plot

