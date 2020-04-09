%This program takes an imported xrd data set (stacked csv) and attempts to quassian fit up to three peaks. Only valide for Matlab r2016b
%This program needs to be improved for curve fitting
dt = 20 ;
t(1) = 0;
lowerq = 0.95;
upperq = 1.1;
[d, lim1] = min(abs(q-lowerq));
[d, lim2] = min(abs(q-upperq));

for i = 1:30;
    %i is length of xrd data set
    f = fit(q(lim1:lim2),xrd(lim1:lim2,i+1)-back(lim1:lim2),'gauss2');
    
    %pull parameters for first peak
    int1(i) = f.a1;
    a1(i) = 2*pi./f.b1;
    strain1(i) = f.c1;
    
    %pull parameters for second peak
    int2(i) = f.a2;
    a2(i) = 2*pi./f.b2;
    strain2(i) = f.c2;
   
    %pull parameters for third peak
    %int3(i) = f.a3;
    %a3(i) = 2*pi./f.b3;
    %strain3(i) = f.c3;
    
    t(i) = 0+dt*(i-1);
end

figure
yyaxis left
plot(t,int1)
yyaxis right
plot (t,a1)

figure
yyaxis left
plot(t,int2)
yyaxis right
plot (t,a3)
    