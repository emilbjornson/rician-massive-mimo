%This Matlab script can be used to generate Figure 6 in the article:
%
%Ozgecan Ozdogan, Emil Bjornson, Erik G. Larsson, “Massive MIMO with
%Spatially Correlated Rician Fading Channels,” IEEE Transactions on
%Communications, To appear.
%
%Download article: https://arxiv.org/abs/1805.07972
%
%This is version 1.0 (Last edited: 2019-02-01)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


%Empty workspace and close figures
close all;
clear;
clc

%Number of BSs
L = 16;

%Number of UEs per BS
K = 10;

%Define the range of BS antennas
Mrange = 100;

%Extract maximum number of BS antennas
Mmax = max(Mrange);

%Define the range of pilot reuse factors
fRange = 1;

%Select the number of setups with random UE locations
nbrOfSetups = 10;

%Select the number of channel realizations per setup
nbrOfRealizations = 40;

%Communication bandwidth
B = 20e6;

%Total uplink transmit power per UE (W)
p = 0.1;

%Noise figure at the BS (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
sigma2=db2pow(noiseVariancedBm-30); %W

%Select length of coherence block
tau_c = 200;
%Angular standard deviation per path in the local scattering model (in degrees)
ASDdeg = 5;
%Power control parameter delta
deltadB=10;
%Pilot Reuse parameter
f=[1,2,4]; %Possible values are 1,2,4,16



%Prepare to save simulation results
sumSE_MonteCarlo_MMSE = zeros(length(Mrange),nbrOfSetups);
sumSE_MonteCarlo_EWMMSE = zeros(length(Mrange),nbrOfSetups);
sumSE_MonteCarlo_LS = zeros(length(Mrange),nbrOfSetups);


%Preallocate matrices for storing the covariance matrices and mean vectors
RNormalized=zeros(Mmax,Mmax,K,L,L);
HMeanNormalized=zeros(Mmax,K,L,L);

%% Go through all setups
for n = 1:nbrOfSetups
    
    %R and HMean is normalized i.e., norm(HMeanNormalized(:,k,l,j))^2=Mmax
    [RNormalized,HMeanNormalized,channelGaindB,ricianFactor,probLOS] = functionExampleSetup(L,K,Mmax,ASDdeg,0);
    
    %Controlled UL channel gain over noise
    channelGainUL= functionPowerControl( channelGaindB,noiseVariancedBm,deltadB,L);
    
    %If all UEs have same power
    % channelGainOverNoise=channelGaindB-noiseVariancedBm;
    
    for m=1:length(Mrange)
        
        
        %Generate the channel realizations and scale the normalized R and HMean
        [R_UL,HMean_UL,H_UL,H_UL_Rayleigh] = functionChannelGeneration( RNormalized(1:Mrange(m),1:Mrange(m),:,:,:)...
            ,HMeanNormalized(1:Mrange(m),:,:,:),channelGainUL,ricianFactor,probLOS,K,L,Mrange(m),nbrOfRealizations);
        
        for s=1:length(f)
            %MMSE Rician SE
            tau_p=f(s)*K;
            %Prelog factor assuming only UL transmission
            prelogFactor=(tau_c -tau_p)/tau_c;
            
            %The UL SE with MMSE for Rician fading
            %Generate the MMSE channel estimates
            Hhat_MMSE_UL = functionChannelEstimateMMSE(R_UL,HMean_UL,H_UL,nbrOfRealizations,Mrange(m),K,L,p,f(s),tau_p);
            %Compute UL SE using use and then forget bound in (19) with MMSE
            SE_MonteCarlo_MMSE = functionMonteCarloSE_UL(Hhat_MMSE_UL,H_UL,prelogFactor,nbrOfRealizations,Mrange(m),K,L,p);
            
            %The UL SE with EW-MMSE for Rician fading
            %Generate the EWMMSE channel estimates
            Hhat_EWMMSE = functionChannelEstimateEWMMSE(R_UL,HMean_UL,H_UL,nbrOfRealizations,Mrange(m),K,L,p,f(s),tau_p);
            %Compute UL SE using use and then forget bound in (19) with EWMMSE
            SE_MonteCarlo_EWMMSE= functionMonteCarloSE_UL(Hhat_EWMMSE,H_UL,prelogFactor,nbrOfRealizations,Mrange(m),K,L,p);
            
            %The UL SE with LS for Rician fading
            %Generate the LS channel estimates
            Hhat_LS = functionChannelEstimateLS(H_UL,nbrOfRealizations,Mrange(m),K,L,p,f(s),tau_p);
            %Compute UL SE using use and then forget bound in (19) with LS
            SE_MonteCarlo_LS = functionMonteCarloSE_UL(Hhat_LS,H_UL,prelogFactor,nbrOfRealizations,Mrange(m),K,L,p);
            
            
            sumSE_MonteCarlo_MMSE(m,n,s) = mean(sum(SE_MonteCarlo_MMSE,1));
            sumSE_MonteCarlo_EWMMSE(m,n,s) = mean(sum(SE_MonteCarlo_EWMMSE,1));
            sumSE_MonteCarlo_LS(m,n,s) = mean(sum(SE_MonteCarlo_LS,1));
            
        end
        
        
        
        %Output simulation progress
        disp([num2str(Mrange(m)) ' antennas of ' num2str(Mmax)]);
        clear R_UL HMean_UL
    end
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    clear RNormalized HMeanNormalized
    
end

%% Plot the figure
MMSE=reshape(mean(sumSE_MonteCarlo_MMSE,2),1,length(f));
EW_MMSE=reshape(mean(sumSE_MonteCarlo_EWMMSE,2),1,length(f));
LS=reshape(mean(sumSE_MonteCarlo_LS,2),1,length(f));

SE=[MMSE(1:3)' EW_MMSE(1:3)' LS(1:3)'];
bar(SE')
ylabel('Average sum SE [bit/s/Hz/cell]');
set(gca,'xticklabel',['  MMSE '; 'EW-MMSE'; '   LS  '])
legend('f=1','f=2','f=4','Location','NorthEast');
colormap(hot);
