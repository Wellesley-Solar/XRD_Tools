function [] = plot_ontop(dat, offset, spacing)
%   plot_ontop plots data from different runs on top of each other
%   No output parameter
%   dat: data array with 1st column being Q and the rest being
%   intensity data
%   offset: vertical spacing between each job on the plot
%   spacing: the spacing between jobs

q=dat(:,1);

[~,c]=size(dat);

figure
hold on
for i=1:spacing:c-1
    plot(q, dat(:,i+1)+offset*i, 'LineWidth',1,'DisplayName',['Job ' num2str(i)]);
end
hold off
legend
xlabel('Q (A^{-1})')
ylabel('Intensity')
end
