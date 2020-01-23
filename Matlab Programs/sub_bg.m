function [n_dat] = sub_bg(dat,lambda)
%   sub_bg process data struct, convert to Q, and subtract background
%   Output: 
%       n_dat: data array with 1st column being Q and the rest being
%   intensity data
%   Input:
%       dat: data struct with data imported from stack csv
%       lambda: OPTIONAL, wavelength of the beam

if ~exist('lambda','var')
 % third parameter does not exist, so default it to something
  lambda = .9744;
end

twotheta=dat.data(:,1)/180*pi;

q = 4*pi*sin(twotheta/2)./lambda;

dat1=dat.data(:,2);


bg=medfilt1(dat1, 100);

[r,c]=size(dat.data);
n_dat=zeros(r,c);
n_dat(:,1)=q;

for i=1:c-1

    temp_dat = dat.data(:,i+1)-bg;
    
    n_dat(:,i+1) = temp_dat;
end

end

