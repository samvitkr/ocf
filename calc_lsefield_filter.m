kcond=64;
load('../data/siren_filters.mat');
fn=sprintf('../data/fieldsq2_lse_k_%03d.mat',kcond);
fnp=sprintf('../data/fieldsq2pos_lse_k_%03d.mat',kcond);
fnn=sprintf('../data/fieldsq2neg_lse_k_%03d.mat',kcond);

m = matfile(fn);

mp = matfile(fnp,'Writable',true);
mp.lse= ifft2(filpos.*fft2(m.lse) ,'symmetric');
clear mp

mn = matfile(fnn,'Writable',true);
mn.lse= ifft2(filneg.*fft2(m.lse) ,'symmetric');
clear mn
%fn4=sprintf('../data/fieldsq4_lse_k_%03d.mat',kcond);
