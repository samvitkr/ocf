kcond=64;
f=sprintf('../data/jointpdf_uw_k_%03d.mat',kcond);
load(f);
fl=sprintf('../data/lse_coeff_k_%03d.mat',kcond);
mf=matfile(fl);

fn=sprintf('../data/fieldsq2_lse_k_%03d.mat',kcond);
mn=matfile(fn,"Writable",true);
mn.lse =uq2.*mf.L1 + wq2.*mf.L2;

fn4=sprintf('../data/fieldsq4_lse_k_%03d.mat',kcond);
mn4=matfile(fn4,"Writable",true);
mn4.lse=uq4.*mf.L1 + wq4.*mf.L2;


