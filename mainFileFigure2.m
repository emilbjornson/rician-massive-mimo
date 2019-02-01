%This Matlab script can be used to generate Figure 2 in the article:
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
nbrOfSetups = 50;

%Select the number of channel realizations per setup
nbrOfRealizations = 100;

%Communication bandwidth
B = 20e6;

%Total uplink transmit power per UE (W)
p = 0.1;
%Noise figure at the BS (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
%sigma2=db2pow(noiseVariancedBm-30); %W

%Select length of coherence block
tau_c = 200;
%Angular standard deviation per path in the local scattering model (in degrees)
ASDdeg = 5;
%Power control parameter delta
deltadB=10;
%Pilot Reuse parameter
f=1;
%Pilot Length
tau_p=K;
%Prelog factor assuming only UL transmission
prelogFactor=(tau_c -tau_p)/tau_c;


%Prepare to save simulation results
userSE_MMSE= zeros(K*L,nbrOfSetups,1);
userSE_EWMMSE= zeros(K*L,nbrOfSetups,1);
userSE_LS= zeros(K*L,nbrOfSetups,1);

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
        
        %The UL SE with MMSE for Rician fading
        %Generate the MMSE channel estimates
        Hhat_MMSE_UL = functionChannelEstimateMMSE(R_UL,HMean_UL,H_UL,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute UL SE using use and then forget bound in (19) with MMSE
        SE_MonteCarlo_MMSE = functionMonteCarloSE_UL(Hhat_MMSE_UL,H_UL,prelogFactor,nbrOfRealizations,Mrange(m),K,L,p);
        
        
        
        %The UL SE with EW-MMSE for Rician fading
        %Generate the EWMMSE channel estimates
        Hhat_EWMMSE = functionChannelEstimateEWMMSE(R_UL,HMean_UL,H_UL,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute UL SE using use and then forget bound in (19) with EWMMSE
        SE_MonteCarlo_EWMMSE= functionMonteCarloSE_UL(Hhat_EWMMSE,H_UL,prelogFactor,nbrOfRealizations,Mrange(m),K,L,p);
        
        
        
        %The UL SE with LS for Rician fading
        %Generate the LS channel estimates
        Hhat_LS = functionChannelEstimateLS(H_UL,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute UL SE using use and then forget bound in (19) with LS
        SE_MonteCarlo_LS = functionMonteCarloSE_UL(Hhat_LS,H_UL,prelogFactor,nbrOfRealizations,Mrange(m),K,L,p);
        
        
        %Store SEs for all UEs
        userSE_MMSE(:,n) = SE_MonteCarlo_MMSE(:);
        userSE_EWMMSE(:,n) = SE_MonteCarlo_EWMMSE(:);
        userSE_LS(:,n) = SE_MonteCarlo_LS(:);
        
        
        %Output simulation progress
        disp([num2str(Mrange(m)) ' antennas of ' num2str(Mmax)]);
        clear R_UL HMean_UL
    end
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    clear RNormalized HMeanNormalized
    
end

%% Plot the figure
CDFnumbers = linspace(0,1,K*L*nbrOfSetups);
figure(1);
hold on; box on;
plot(sort(userSE_MMSE(:)),CDFnumbers,'r-');
plot(sort(userSE_EWMMSE(:)),CDFnumbers,'k--');
plot(sort(userSE_LS(:)),CDFnumbers,'b-.');
xlabel('SE per UE [bit/s/Hz]');
ylabel('Cumulative Distribution Function (CDF)');
legend('MMSE estimator','EW-MMSE estimator','LS estimator','Location','SouthEast');
