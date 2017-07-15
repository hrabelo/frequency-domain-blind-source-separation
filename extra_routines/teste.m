% Para testes aleatórios

SIR = cell(7,18);
SDR = cell(7,18);
SAR = cell(7,18);
ITER = cell(7,18);

rtimes = [0 0.05 0.1 0.15 0.2 0.25 0.5];

frr = 1;
K = 4096;
windlen = 2048;
smoothFlag = 0; 
dist_scrmic = 1;                        
ang_src(1) = 50;      
ang_src(2) = 95;      
ang_src(3) = 135;    
N = 3;
M = 3;
winda = wind{2};
J = Jc{2};
mics_struct = 'line2d';
reverb_time = rtimes(frr);
perm_meth = 'supervised';
corr_env = 'PowValue2';
selected_combs = [8];
mixing_mode = 'ismbuild';
tmpi = 1;
bss_test;

perm_meth = 'tdoa';
tmpi = 2;
bss_test;

perm_meth = 'doa_conjcorr';
tmpi = 3;
bss_test;

perm_meth = 'doa_harmcorr';
tmpi = 4;
bss_test;

perm_meth = 'doa_allcorr';
tmpi = 5;
bss_test;

perm_meth = 'conjcorr';
tmpi = 6;
bss_test;

perm_meth = 'allcorr';
tmpi = 7;
bss_test;

for frr = 2:7
    
    reverb_time = rtimes(frr);
    mixing_mode = 'ismbuild';
    perm_meth = 'supervised';
    tmpi = 1;
    bss_test;

    perm_meth = 'tdoa';
    tmpi = 2;
    bss_test;

    perm_meth = 'doa_conjcorr';
    tmpi = 3;
    bss_test;

    perm_meth = 'doa_harmcorr';
    tmpi = 4;
    bss_test;

    perm_meth = 'doa_allcorr';
    tmpi = 5;
    bss_test;

    perm_meth = 'conjcorr';
    tmpi = 6;
    bss_test;

    perm_meth = 'allcorr';
    tmpi = 7;
    bss_test;
    
end
   
%{
perm_meth = 'conjcorr';
corr_env = 'PowValue2';
tmpi = 3;
bss_test;

perm_meth = 'harmcorr';
corr_env = 'AbsValue';
tmpi = 4;
bss_test;

perm_meth = 'harmcorr';
corr_env = 'PowValue2';
tmpi = 5;
bss_test;

perm_meth = 'globalcorr';
corr_env = 'AbsValue';
tmpi = 6;
bss_test;

perm_meth = 'globalcorr';
corr_env = 'PowValue2';
tmpi = 7;
bss_test;

perm_meth = 'localcorr';
corr_env = 'AbsValue';
tmpi = 8;
bss_test;

perm_meth = 'localcorr';
corr_env = 'PowValue2';
tmpi = 9;
bss_test;

perm_meth = 'allcorr';
corr_env = 'AbsValue';
tmpi = 10;
bss_test;

perm_meth = 'allcorr';
corr_env = 'PowValue';
tmpi = 11;
bss_test;

perm_meth = 'allcorr';
corr_env = 'SPDValue';
tmpi = 12;
bss_test;

perm_meth = 'allcorr';
corr_env = 'PowValue2';
tmpi = 13;
bss_test;
%}


%{
winda = wind{8};
J = Jc{8};
tmpi = 2;
bss_test;

winda = wind{2};
J = Jc{2};
mics_struct = 'line2d';
reverb_time = 0.15;
perm_meth = 'tdoa';
selected_combs = [13 22 45];
mixing_mode = 'ismbuild';
tmpi = 1;
bss_test;
%}

%{
reverb_time = 0.1;
mixing_mode = 'ismbuild';
tmpi = 2;
bss_test;

reverb_time = 0.15;
mixing_mode = 'ismbuild';
tmpi = 3;
bss_test;

reverb_time = 0.2;
mixing_mode = 'ismbuild';
tmpi = 4;
bss_test;

reverb_time = 0.25;
mixing_mode = 'ismbuild';
tmpi = 5;
bss_test;

reverb_time = 0.3;
mixing_mode = 'ismbuild';
tmpi = 6;
bss_test;

reverb_time = 0.35;
mixing_mode = 'ismbuild';
tmpi = 7;
bss_test;

reverb_time = 0.4;
mixing_mode = 'ismbuild';
tmpi = 8;
bss_test;

reverb_time = 0.45;
mixing_mode = 'ismbuild';
tmpi = 9;
bss_test;

reverb_time = 0.5;
mixing_mode = 'ismbuild';
tmpi = 10;
bss_test;

reverb_time = 0.55;
mixing_mode = 'ismbuild';
tmpi = 11;
bss_test;

reverb_time = 0.6;
mixing_mode = 'ismbuild';
tmpi = 12;
bss_test;

reverb_time = 0.65;
mixing_mode = 'ismbuild';
tmpi = 13;
bss_test;

reverb_time = 0.7;
mixing_mode = 'ismbuild';
tmpi = 13;
bss_test;
%}

%{
K = 2048;
windlen = 2048;
smoothFlag = 0; 
dist_scrmic = 1;                        
winda = wind{2};
J = Jc{2};
ang_src(1) = 50;      
ang_src(2) = 95;      
ang_src(3) = 135;    
N = 3;
M = 3;
mics_struct = 'line2d';
reverb_time = 0.15;
perm_meth = 'doa_conjcorr';
corr_env = 'PowValue2';
selected_combs = [8 17 25];
mixing_mode = 'ismbuild';
tmpi = 10;
bss_test;

perm_meth = 'doa_globalcorr';
tmpi = 11;
bss_test;

mics_struct = 'cluster2d';
perm_meth = 'tdoa';
mixing_mode = 'ismbuild';
tmpi = 12;
bss_test;

perm_meth = 'supervised';
tmpi = 13;
bss_test;


smoothFlag = 1;
smooth_filter = [0.25 0.5 0.25];
tmpi = 14;
bss_test;

smooth_filter = [0.003 0.0602 0.2516 0.3902 0.2516 0.0602 0.003];
tmpi = 15;
bss_test;

smooth_filter = [0.01 0.0817 0.24 0.3363 0.24 0.0817 0.01];
tmpi = 16;
bss_test;

smooth_filter = [0.0092 0.0795 0.2407 0.3409 0.2407 0.0795 0.0092];
tmpi = 17;
bss_test;

smooth_filter = [0.0014 0.0032 0.0129 0.9787 0.0129 0.0032 0.0014];
tmpi = 18;
bss_test;
%}


%{
tmp = [y(1, 1:10000); y(1, 501:10500)];
mR = mean(tmp, 2);
vR = var(tmp, 0, 2);

RR = fast_corr2by2(tmp, tmp, mR, mR, vR, vR);
disp(RR)
%}

%{
tmp1 = squeeze(S(1:K/2+1, 1, :));
mR1 = mean(tmp1, 2);
vR1 = var(tmp1, 0, 2);

tmp2 = squeeze(S(1:K/2+1, 2, :));
mR2 = mean(tmp2, 2);
vR2 = var(tmp2, 0, 2);

RR = fast_corr2by2(tmp1, tmp2, mR1, mR2);

figure;
colormap('gray');
imagesc(RR);  set(gca, 'Ydir', 'normal');
%}

%{
TeTa = 0:pi/360:pi;
AnGlE = 0:0.5:180;
freq = [0 : fs/K : (fs/K)*(K/2-1), -fs/2 : fs/K : -fs/K];

testing_freqs = floor(rand(1,9)*256);
for cf = 1:9
    Wtmp = squeeze(W(testing_freqs(cf), :, :));
    U1 = abs(sum(diag(Wtmp(1,:)) * exp(i*2*pi*freq(testing_freqs(cf))*(1/343)*[0; 0.04]*cos(TeTa)), 1));
    [dummy ind] = min(U1);
    dirpat(1,1) = AnGlE(ind);
  
    U2 = abs(sum(diag(Wtmp(2,:)) * exp(i*2*pi*freq(testing_freqs(cf))*(1/343)*[0; 0.04]*cos(TeTa)), 1));
    [dummy ind] = min(U2);
    dirpat(2,1) = AnGlE(ind);
    
    Atmp = inv(Wtmp);
    mak = acos(angle(Atmp(1,:) ./ Atmp(2,:)) / (2*pi*freq(testing_freqs(cf))*(1/343)*0.04))*180/pi;
    mak2 = acos(angle(-Wtmp(:,1) ./ Wtmp(:,2)) / (2*pi*freq(testing_freqs(cf))*(1/343)*0.04)).'*180/pi;
    
    subplot(3, 3, cf);
    plot(AnGlE, U1);    xlabel(freq(testing_freqs(cf)));    title(mak2);
    ylabel(mak);                                            hold;       
    plot(AnGlE, U2, 'r');                                   hold;
end

clear Wtmp Atmp mak mak2 dirpat U1 U2 TeTa AnGlE freq testing_freqs

%}