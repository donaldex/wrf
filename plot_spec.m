%%script 2

function x=plot_spec(u)

%%%%%energyspec_wrf.m is used for getting energy spectrum k=1:N/2
%%%this script is for plotting....k has to *N/distance


count=1;
a=energyspec_wrf(u,370,4,21);
b=energyspec_wrf(u,421,5,29);
c=energyspec_wrf(u,377,4,23);
d=energyspec_wrf(u,421,5,25);
e=energyspec_wrf(u,409,4,27);
for l=1:5
    for m=1:5
a=a+energyspec_wrf(u,370+l,4+m,21);
b=b+energyspec_wrf(u,421+l,5+m,29);
c=c+energyspec_wrf(u,377+l,4+m,23);
d=d+energyspec_wrf(u,395+l,4+m,25);
e=e+energyspec_wrf(u,409+l,4+m,27);
count=count+1;
    end
end
x=(a+b+c+d+e)/count/5;

end