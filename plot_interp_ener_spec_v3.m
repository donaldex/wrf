function x=plot_interp_ener_spec_v3


a=interp_energy_spec_v3(18);
b=interp_energy_spec_v3(19);
c=interp_energy_spec_v3(20);
d=interp_energy_spec_v3(21);
e=interp_energy_spec_v3(22);
f=interp_energy_spec_v3(23);

y=(a+b+c+d+e+f)/6;

N=size(y);
N=N(2);

k=[0:N/2];
loglog(k,y)

title('energy spectrum');
xlabel('k')
ylabel('Ek_w')


end