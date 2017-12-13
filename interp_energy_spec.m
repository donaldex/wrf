function x=interp_energy_spec
%%indepedent program

w_in=ncread('wrfout_d01_0001-01-01_00:00:00','W');
ph=ncread('wrfout_d01_0001-01-01_00:00:00','PH');
phb=ncread('wrfout_d01_0001-01-01_00:00:00','PHB');
z=(ph+phb)/9.8;
time=21;
w=w_in(:,:,:,time);
z=z(:,:,:,time);
wq=interp3(w);
zq=interp3(z);
[nx,ny,nz,nt]=size(wq);
avgy_wq=zeros(nx,nz);
flu_wq=zeros(nx,ny,nz);
sumy_wq=zeros(nx,nz);
%%%%%%%%%%%%%%points at 5km
pt = find(zq>5000 & zq<5025);
k_index=floor(pt/nx/ny)+1;
j_index=floor((pt-(k_index-1)*nx*ny)/ny)+1;
i_index=pt-(k_index-1)*nx*ny-(j_index-1)*ny;



for i=1:nx
for k=1:nz
    avgy_wq(i,k)=sum(wq(i,:,k))/ny;
end 
end

for j=1:ny
    for i=1:nx
        for k=1:nz
   %flu_w(:,j,:)=w(:,j,:)-avgy_w
   flu_wq(i,j,k)=wq(i,j,k)-avgy_wq(i,k);
        end
    end
end

%%%%%%%%%%%%%finding the max flu point.
flu_w_line=flu_wq.^2;

for i=1:nx
    for k=1:nz
sumy_wq(i,k)=sum(flu_w_line(i,:,k));
    end
end

maxi=0;
maxk=0;
max_sumy_wq=0;

for i=1:size(i_index)
   for k=1:size(k_index)       
      if sumy_wq(i_index(i),k_index(k))>max_sumy_wq
          maxi=i_index(i);
          maxk=k_index(k);
          max_sumy_wq=sumy_wq(i_index(i),k_index(k));
      end       
   end    
end



%%%%%%%%%%%%%%%%%%%%%%%FFT
uy1=wq(maxi,:,maxk);
N=size(uy1);
N=N(2);

uhat1 = fft(uy1)/N;
x=sum(uy1.^2)/N-sum(abs(uhat1).^2);
%x=sum(uy1.^2)-sum(abs(uhat1).^2)/N;
Ek1 = 0.5*abs(uhat1(1:N/2+1)).^2;
Ek1(2:N/2) = Ek1(2:N/2) + 0.5*abs(uhat1(N:-1:N/2+2)).^2;
a=Ek1;

k=[0:N/2];
loglog(k,a)
%logkhalf,b,'b',logkhalf,c,'g');
title('energy spectrum');
xlabel('k')
ylabel('Ek_w')




end