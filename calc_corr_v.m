
close all
clear

Nx=2048;
Ny=512;
Nz=1536;
jcond=105
jct=Ny-jcond+1

nf=1;
phivv=          single(zeros(Nz,Nx,Ny/2));
phivu=          single(zeros(Nz,Nx,Ny/2));
phivw=          single(zeros(Nz,Nx,Ny/2));

phivdudx=       single(zeros(Nz,Nx,Ny/2));
phivdvdx=       single(zeros(Nz,Nx,Ny/2));
phivdwdx=       single(zeros(Nz,Nx,Ny/2));

phivdudy=       single(zeros(Nz,Nx,Ny/2));
phivdvdy=       single(zeros(Nz,Nx,Ny/2));
phivdwdy=       single(zeros(Nz,Nx,Ny/2));

phivdudz=       single(zeros(Nz,Nx,Ny/2));
phivdvdz=       single(zeros(Nz,Nx,Ny/2));
phivdwdz=       single(zeros(Nz,Nx,Ny/2));

%phivfx=        single( zeros(Nz,Nx,Ny/2));
phivvoz=        single(zeros(Nz,Nx,Ny/2));
phivwoy=        single(zeros(Nz,Nx,Ny/2));

tstart=1;
tend=10;
tstep=1;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')

j=1:Ny/2;
jcr=j+Ny/2;
for time=tstart:tstep:tend
%for time=1:1
        time
        fvel=sprintf("../data/velfieldpar_%02d.mat",time);
        m=matfile(fvel)
        vfj=single(m.vF(:,:,jcond));
        vfj(1,1)=0;
        vfjt=single(m.vF(:,:,jct));
        vfjt(1,1)=0;
        tic
%       phivv(:,:,j)=phivv(:,:,j)+conj(vfj).*single(m.vF(1:Nz,1:Nx,j));
%        phivu(:,:,j)=phivu(:,:,j)+conj(vfj).*single(m.uF(1:Nz,1:Nx,j));
%        phivw(:,:,j)=phivw(:,:,j)+conj(vfj).*single(m.wF(1:Nz,1:Nx,j));
        phivv(:,:,j)=phivv(:,:,j)+conj(vfj).*single(m.vF(1:Nz,1:Nx,j))+flip( conj(vfjt).*single(m.vF(1:Nz,1:Nx,jcr)),3);
        phivu(:,:,j)=phivu(:,:,j)+conj(vfj).*single(m.uF(1:Nz,1:Nx,j))-flip( conj(vfjt).*single(m.uF(1:Nz,1:Nx,jcr)),3);
        phivw(:,:,j)=phivw(:,:,j)+conj(vfj).*single(m.wF(1:Nz,1:Nx,j))-flip( conj(vfjt).*single(m.wF(1:Nz,1:Nx,jcr)),3);
        toc
        clear m

        fvelgx=sprintf("../data/velgradx_%03d.mat",time);
        mgx=matfile(fvelgx)
%       phivdudx(:,:,j)=phivdudx(:,:,j)+conj(vfj).*single(mgx.dudxF(1:Nz,1:Nx,j));
%       phivdvdx(:,:,j)=phivdvdx(:,:,j)+conj(vfj).*single(mgx.dvdxF(1:Nz,1:Nx,j));
%       phivdwdx(:,:,j)=phivdwdx(:,:,j)+conj(vfj).*single(mgx.dwdxF(1:Nz,1:Nx,j));
        phivdudx(:,:,j)=phivdudx(:,:,j)+conj(vfj).*single(mgx.dudxF(1:Nz,1:Nx,j))-flip(conj(vfjt).*single(mgx.dudxF(1:Nz,1:Nx,jcr)),3);
        phivdvdx(:,:,j)=phivdvdx(:,:,j)+conj(vfj).*single(mgx.dvdxF(1:Nz,1:Nx,j))+flip(conj(vfjt).*single(mgx.dvdxF(1:Nz,1:Nx,jcr)),3);
        phivdwdx(:,:,j)=phivdwdx(:,:,j)+conj(vfj).*single(mgx.dwdxF(1:Nz,1:Nx,j))-flip(conj(vfjt).*single(mgx.dwdxF(1:Nz,1:Nx,jcr)),3);
        toc
        clear mgx

        fvelgy=sprintf("../data/velgrady_%03d.mat",time);
        mgy=matfile(fvelgy)
%       phivdudy(:,:,j)=phivdudy(:,:,j)+conj(vfj).*single(mgy.dudyF(1:Nz,1:Nx,j));
%       phivdvdy(:,:,j)=phivdvdy(:,:,j)+conj(vfj).*single(mgy.dvdyF(1:Nz,1:Nx,j));
%       phivdwdy(:,:,j)=phivdwdy(:,:,j)+conj(vfj).*single(mgy.dwdyF(1:Nz,1:Nx,j));
        phivdudy(:,:,j)=phivdudy(:,:,j)+conj(vfj).*single(mgy.dudyF(1:Nz,1:Nx,j))+flip(conj(vfjt).*single(mgy.dudyF(1:Nz,1:Nx,jcr)),3);
        phivdvdy(:,:,j)=phivdvdy(:,:,j)+conj(vfj).*single(mgy.dvdyF(1:Nz,1:Nx,j))-flip(conj(vfjt).*single(mgy.dvdyF(1:Nz,1:Nx,jcr)),3);
        phivdwdy(:,:,j)=phivdwdy(:,:,j)+conj(vfj).*single(mgy.dwdyF(1:Nz,1:Nx,j))+flip(conj(vfjt).*single(mgy.dwdyF(1:Nz,1:Nx,jcr)),3);
        toc
        clear mgy

        fvelgz=sprintf("../data/velgradz_%03d.mat",time);
        mgz=matfile(fvelgz)

%       phivdudz(:,:,j)=phivdudz(:,:,j)+conj(vfj).*single(mgz.dudzF(1:Nz,1:Nx,j));
%       phivdvdz(:,:,j)=phivdvdz(:,:,j)+conj(vfj).*single(mgz.dvdzF(1:Nz,1:Nx,j));
%       phivdwdz(:,:,j)=phivdwdz(:,:,j)+conj(vfj).*single(mgz.dwdzF(1:Nz,1:Nx,j));

        phivdudz(:,:,j)=phivdudz(:,:,j)+conj(vfj).*single(mgz.dudzF(1:Nz,1:Nx,j))-flip(conj(vfjt).*single(mgz.dudzF(1:Nz,1:Nx,jcr)),3);
        phivdvdz(:,:,j)=phivdvdz(:,:,j)+conj(vfj).*single(mgz.dvdzF(1:Nz,1:Nx,j))+flip(conj(vfjt).*single(mgz.dvdzF(1:Nz,1:Nx,jcr)),3);
        phivdwdz(:,:,j)=phivdwdz(:,:,j)+conj(vfj).*single(mgz.dwdzF(1:Nz,1:Nx,j))-flip(conj(vfjt).*single(mgz.dwdzF(1:Nz,1:Nx,jcr)),3);
        toc
        clear mgz

        ft=sprintf("../data/Transfer_%03d.mat",time);
        mt=matfile(ft)

%       phivvoz(:,:,j)=phivvoz(:,:,j)+conj(vfj).*single(mt.vozF(1:Nz,1:Nx,j));
%       phivwoy(:,:,j)=phivwoy(:,:,j)+conj(vfj).*single(mt.woyF(1:Nz,1:Nx,j));

        phivvoz(:,:,j)=phivvoz(:,:,j)+conj(vfj).*single(mt.vozF(1:Nz,1:Nx,j))-flip(conj(vfjt).*single(mt.vozF(1:Nz,1:Nx,jcr)),3);
        phivwoy(:,:,j)=phivwoy(:,:,j)+conj(vfj).*single(mt.woyF(1:Nz,1:Nx,j))-flip(conj(vfjt).*single(mt.woyF(1:Nz,1:Nx,jcr)),3);

        clear mt vfj vfjt
        toc
end


N=single((Nx*Nz)/(2*nf));

fn=sprintf('../data/corr_v_reflect_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true)
%mf=matfile('testmem.mat','Writable',true)
%for
jw =1:Ny/2;
%jw
        mf.Rvv(1:Nz,1:Nx,jw)=ifft2(single((phivv(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivv;
        mf.Rvu(1:Nz,1:Nx,jw)=ifft2(single((phivu(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivu;
        mf.Rvw(1:Nz,1:Nx,jw)=ifft2(single((phivw(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivw;

        mf.Rvdudx(1:Nz,1:Nx,jw)=ifft2(single((phivdudx(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivdudx;
        mf.Rvdvdx(1:Nz,1:Nx,jw)=ifft2(single((phivdvdx(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivdvdx;
        mf.Rvdwdx(1:Nz,1:Nx,jw)=ifft2(single((phivdwdx(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivdwdx;

        mf.Rvdudy(1:Nz,1:Nx,jw)=ifft2(single((phivdudy(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivdudy;
        mf.Rvdvdy(1:Nz,1:Nx,jw)=ifft2(single((phivdvdy(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivdvdy;
        mf.Rvdwdy(1:Nz,1:Nx,jw)=ifft2(single((phivdwdy(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivdwdy;

        mf.Rvdudz(1:Nz,1:Nx,jw)=ifft2(single((phivdudz(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivdudz;
        mf.Rvdvdz(1:Nz,1:Nx,jw)=ifft2(single((phivdvdz(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivdvdz;
        mf.Rvdwdz(1:Nz,1:Nx,jw)=ifft2(single((phivdwdz(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivdwdz;

        mf.Rvvoz(1:Nz,1:Nx,jw)= ifft2(single((phivvoz(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivvoz;
        mf.Rvwoy(1:Nz,1:Nx,jw)= ifft2(single((phivwoy(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivwoy;
%end
mf.j=jcond;
toc
~
~
~
~
~
~
~
~
~
