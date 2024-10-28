clc
clear
%% Connectivity using ISPC & PLI
load EEGdata

% compute Laplacian
EEG.lap = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
% two channels for computing synchronization
chan1 = 'FCz';
chan2 = 'POz';
minfreq = 8;
maxfreq = 30;
numfrex = 50;

fwhm_range = [ 2/minfreq 4/maxfreq ];
frex  = linspace(minfreq,maxfreq,numfrex);
fwhms = linspace(fwhm_range(1),fwhm_range(2),numfrex);

wtime = -1:1/EEG.srate:1;
nWave = length(wtime);
nData = EEG.pnts*EEG.trials;
nConv = nWave+nData-1;
halfw = (length(wtime)-1)/2;


chan1_fft = fft( reshape(EEG.data(strcmpi(chan1,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);
chan2_fft = fft( reshape(EEG.data(strcmpi(chan2,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);
chan1Lap_fft = fft( reshape(EEG.lap (strcmpi(chan1,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);
chan2Lap_fft = fft( reshape(EEG.lap (strcmpi(chan2,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);

[ispc,pli] = deal( zeros(2,numfrex,EEG.pnts) );

for fi=1:numfrex
    cmw  = exp(2*1i*pi*frex(fi).*wtime) .* exp(-4*log(2)*wtime.^2./fwhms(fi)^2);
    cmwX = fft(cmw,nConv);

    % chan1
    sig1 = ifft(cmwX.*chan1_fft,nConv);
    sig1 = sig1(halfw+1:end-halfw);
    sig1 = reshape(sig1,EEG.pnts,EEG.trials);
    
    % chan2
    sig2 = ifft(cmwX.*chan2_fft,nConv);
    sig2 = sig2(halfw+1:end-halfw);
    sig2 = reshape(sig2,EEG.pnts,EEG.trials);
    
    phasediffVOLT = exp(1i*( angle(sig1)-angle(sig2) ));
    
    sig1 = ifft(cmwX.*chan1Lap_fft,nConv);
    sig1 = sig1(halfw+1:end-halfw);
    sig1 = reshape(sig1,EEG.pnts,EEG.trials);
 
    sig2 = ifft(cmwX.*chan2Lap_fft,nConv);
    sig2 = sig2(halfw+1:end-halfw);
    sig2 = reshape(sig2,EEG.pnts,EEG.trials);
    
    phasediff = exp(1i*( angle(sig1)-angle(sig2) ));
    
    % ISPC and PLI for voltage
    ispc(1,fi,:) = abs(mean(phasediffVOLT,2));
    pli(1,fi,:)  = abs(mean(sign(imag(phasediffVOLT)),2));
    
    % ISPC and PLI for laplacian
    ispc(2,fi,:) = abs(mean(phasediff,2));
    pli(2,fi,:)  = abs(mean(sign(imag(phasediff)),2));
    
end

%% plotting

figure(1), clf
colormap hot
clim = [.1 .4];
labels = {'Voltage';'Laplacian'};

for i=1:2
    
    % ISPC
    subplot(2,2,i)
    contourf(EEG.times,frex,squeeze(ispc(i,:,:)),40,'linecolor','none')
    set(gca,'xlim',[-300 1000],'clim',clim)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title([ 'ISPC: ' labels{i} ])
    
    
    % PLI
    subplot(2,2,i+2)
    contourf(EEG.times,frex,squeeze(pli(i,:,:)),40,'linecolor','none')
    set(gca,'xlim',[-300 1000],'clim',clim)
    caxis(clim)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title([ 'PLI: ' labels{i} ])
end

