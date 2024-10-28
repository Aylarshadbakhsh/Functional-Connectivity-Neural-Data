clc
clear
close all
%% Connectivity analyses
load v1_laminar.mat
minfreq = 10;
maxfreq = 80; 
numfrex = 100;
frex = linspace(minfreq,maxfreq,numfrex);
widt = linspace(5,10,numfrex);
% time window
tidx=dsearchn(timevec',[0,1]');
[corrmat1,corrmat2] = deal(zeros(numfrex,size(csd,1),size(csd,1)));
for fi=1:numfrex
    supdat = reshape(csd,size(csd,1),[]);
    hdat = hilbert( filterFGx( supdat,srate,frex(fi),widt(fi) )' ).';
    hdat = reshape(hdat,size(csd));
    hdat = hdat(:,tidx(1):tidx(2),:);
    powdat = abs(hdat);
    phsdat = angle(hdat);
    for triali=1:size(powdat,3)
        
        %Power correlation
        corrmat1(fi,:,:) = squeeze(corrmat1(fi,:,:)) + corr(squeeze(powdat(:,:,triali))','type','s');
        
        %phase synchronization
        for chani=1:size(csd,1)
            for chanj=1:size(csd,1)
            
            % Extract the angles (phase data) for each channel and trial
                 tmpAi = phsdat(chani,:,triali); % phase of channel i
                 tmpAj = phsdat(chanj,:,triali); % phase of channel j
            
            % Compute phase synchronization for each trial
                 trialsynch = abs(mean(exp(1i*(tmpAi - tmpAj)), 2)); 
            
            % Add synchronization value to the matrix
                  corrmat2(fi,chani,chanj) = squeeze(corrmat2(fi,chani,chanj)) + trialsynch;
            end
       end
    end
end 
corrmat1 = corrmat1 / triali;
corrmat2 = corrmat2 / triali;
% 3 specific frequencies to display
frex2plot = [ 12 41 56 ];


figure(1),
subplot(311), hold on
plot(frex,mean( reshape(corrmat1(:,:,:),numfrex,[]).^2 ,2),'ks-','linewid',2,'markerfacecolor','w')
plot(frex,mean( reshape(corrmat2(:,:,:),numfrex,[]).^2 ,2),'ro-','linewid',2,'markerfacecolor','w')
xlabel('Frequency (Hz)')
ylabel('Correlation or synchronization')
legend({'Power Correlation';'Phase Synchronization'})

for fi=1:3
    
    fidx = dsearchn(frex',frex2plot(fi));
    
    % power correlations
    subplot(3,3,3+fi)
    imagesc(squeeze(corrmat1(fidx,:,:)))
    xlabel('Channels'), ylabel('Channels')
    set(gca,'clim',[-1 1])
    axis square
    title([ 'Power corrs.: ' num2str(round(frex(fidx))) ' Hz' ])
    
    
    % phase synchronization
    subplot(3,3,6+fi)
    imagesc(squeeze(corrmat2(fidx,:,:)))
    xlabel('Channels'), ylabel('Channels')
    set(gca,'clim',[0 1])
    axis square
    title([ 'Phase synch.: ' num2str(round(frex(fidx))) ' Hz' ])
end

