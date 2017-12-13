%%script 2.2

function x=plot_spec2(w)

%%%%%energyspec_wrf.m is used for getting energy spectrum k=1:N/2
%%%this script is for plotting....k has to *N/distance
%%%time?
time=21;
ph=ncread('wrfout_d01_0001-01-01_00:00:00','PH');
phb=ncread('wrfout_d01_0001-01-01_00:00:00','PHB');
z=(ph+phb)/9.8;
z=z(:,:,:,time);  %%time 21
n=size(z);
nx=n(1);
nz=n(3);
ny=n(2);
%[row,col,dep] = find(z>4999 & z<5001);
pt = find(z>4000 & z<5000)

avgy_w=zeros(nx,nz);
flu_w=zeros(nx,ny,nz);
for i=1:nx
for k=1:nz
    avgy_w(i,k)=sum(w(i,:,k))/ny;
end 
end

for j=1:ny
    for i=1:nx
        for k=1:nz
   %flu_w(:,j,:)=w(:,j,:)-avgy_w
   flu_w(i,j,k)=w(i,j,k)-avgy_w(i,k);
        end
    end
end
%maxind=[0,0,0];
%ncol=size(col);
%ncol=ncol(2);
%nrow=size(row);
%nrow=nrow(2);
%ndep=size(dep);
%ndep=ndep(2);
%maxflu=0;
%for r=1:nrow
%    for c=1:ncol
%        for d=1:ndep
%        if (flu_w(row(r),col(c),dep(d)))^2>maxflu
%        maxind=[r,c,d];
%        maxflu=(flu_w(row(r),col(c),dep(d)))^2;
%        end       
%        end
%    end
%end

%[r,c,d]=find(flu_w(row,col,dep)==max(max(max(flu_w(row,col,dep)))));


flu_w_line=flu_w.^2;
flu_w_line=flu_w_line(:);
pt2=find(flu_w_line(pt)==max(flu_w_line(pt)));
%%max is flu_w_line(pt(pt2))
%ind=find(flu_w.^2==flu_w_line(pt(pt2(1))));
ind=pt(pt2(1))

dep=floor(ind/nx/ny)+1
col=floor((ind-(dep-1)*nx*ny)/ny)+1
row=ind-(dep-1)*nx*ny-(col-1)*ny
x=energyspec_wrf(w,row,dep,time);


end