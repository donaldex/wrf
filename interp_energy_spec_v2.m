function x=interp_energy_spec_v2
%%indepedent program

w_in=ncread('wrfout_d01_0001-01-01_00:00:00','W',[1,1,1,21],[576,144,101,1]);
ph=ncread('wrfout_d01_0001-01-01_00:00:00','PH',[1,1,1,21],[576,144,101,1]);
phb=ncread('wrfout_d01_0001-01-01_00:00:00','PHB',[1,1,1,21],[576,144,101,1]);
z=(ph+phb)/9.81;
%time=21;
deltax=1000;
[nx,ny,nz,nt]=size(w_in);
%N=[1:nx]
%M=[1:ny]
%P=[1:nz]
w=w_in(:,:,:);
z=z(:,:,:);
z=double(z);
w=double(w);

x1d=[0:576-1]'*1000;
y1d=[0:144-1]'*1000;
z1d=[0:101-1]'*100;
[X1,Y1,Z1]=meshgrid(x1d,y1d,z1d);


XQ=X1(:,:,1);
YQ=Y1(:,:,1);
ZQ=X1(:,:,1);
ZQ(:,:)=5000;

wq=interp3(X1,Y1,z,w,XQ,YQ,ZQ);
%[nx,ny,nz]=size(wq);
avgy_wq=zeros(nx);
flu_wq=zeros(nx,ny);
sumy_wq=zeros(nx);

%%%%%%%%%%%%%%points at 5km
%pt = find(zq>5000 & zq<5025);
%k_index=floor(pt/nx/ny)+1;
%j_index=floor((pt-(k_index-1)*nx*ny)/ny)+1;
%i_index=pt-(k_index-1)*nx*ny-(j_index-1)*ny;



for i=1:nx
    avgy_wq(i)=sum(wq(i,:))/ny;
end 

for j=1:ny
    for i=1:nx
   %flu_w(:,j,:)=w(:,j,:)-avgy_w
   flu_wq(i,j)=wq(i,j)-avgy_wq(i);
    end
end

%%%%%%%%%%%%%finding the max flu point.
flu_w_line=flu_wq.^2;

for i=1:nx
sumy_wq(i)=sum(flu_w_line(i,:));
end

maxi=0;
max_sumy_wq=0;

for i=1:nx    
      if sumy_wq(i)>max_sumy_wq
          maxi=i;
          max_sumy_wq=sumy_wq(i);
      end
end



%%%%%%%%%%%%%%%%%%%%%%%FFT
uy1=wq(maxi,:);
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