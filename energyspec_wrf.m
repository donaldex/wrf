%%%script 1

function a=energyspec_wrf(u,x1,z1,t)

uy1=u(x1,:,z1,t);
N=size(uy1);
N=N(2);

uhat1 = fft(uy1)/N;
x=sum(uy1.^2)/N-sum(abs(uhat1).^2)
%x=sum(uy1.^2)-sum(abs(uhat1).^2)/N;
Ek1 = 0.5*abs(uhat1(1:N/2+1)).^2;
Ek1(2:N/2) = Ek1(2:N/2) + 0.5*abs(uhat1(N:-1:N/2+2)).^2;
a=Ek1;
%%plot(abs(uhat));
%%%%plot(logk(1:N/2+1),log(Ek1));
end