clear all; close all; clc;
n1=174;
n2=500;
o1=0.0;
o2=0.0;
d1=20.0;
d2=20.0;
npml=10;
unit=ones(1,npml);
n1pml=n1+npml;
n2pml=n2;
%filein=input("Input file: ","%s");
filein='marmousi_II_marine.rho';
%fileout=input("Outnput file: ","%s");
fileout='marmousi_II_marine_exttop.rho';

modin=zeros(n1,n2);
fid=fopen(filein,"rb");
modin=fread(fid,size(modin),'float');
figure(1)
vmax1=max(max(modin))
vmin1=min(min(modin))
imagesc(modin);
colorbar;
colormap("jet");
caxis([vmin1 vmax1]);

modout=zeros(n1pml,n2pml);
modout(npml+1:n1pml,1:n2pml)=modin;
modout(1:npml,:)=modin(1,1);

figure(1)
vmax2=max(max(modout))
vmin2=min(min(modout))
imagesc(modout);
colorbar;
colormap("jet");
caxis([vmin2 vmax2]);
axis equal;

fid=fopen(fileout,'wb');
fwrite(fid,modout,'float');
fclose(fid);
