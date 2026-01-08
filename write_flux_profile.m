
load('ocf_test.mat')
vozm =squeeze(mean(voz,[ 1 2]));
woym =squeeze(mean(woy,[ 1 2]));
viscm=squeeze(mean(visc,[ 1 2]));
fid = fopen('out.txt','w');
fprintf(fid, '%.5f %.5f %.5f %.5f\n', [voz, woy, visc, woy-voz+visc].');
fclose(fid);

