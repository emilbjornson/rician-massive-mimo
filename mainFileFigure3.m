%This Matlab script can be used to generate Figure 3 in the article:
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
Mrange = 10:10:100;

%Extract maximum number of BS antennas
Mmax = max(Mrange);

%Define the range of pilot reuse factors
fRange = 1;

%Select the number of setups with random UE locations
nbrOfSetups = 20;

%Select the number of channel realizations per setup
nbrOfRealizations = 50;

%Communication bandwidth
B = 20e6;

%Total uplink pilot power per UE (W)
p = 0.1;

%Downlink pilot power per UE (W)
rho=0.1;

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
tau_p=f*K;

%Compute the prelog factor assuming only downlink transmission
prelogFactor=(tau_c -tau_p)/tau_c;



%Prepare to save simulation results for Theoretical and Monte-Carlo
sumSE_MonteCarlo_MMSE = zeros(length(Mrange),nbrOfSetups);
sumSE_MonteCarlo_EWMMSE = zeros(length(Mrange),nbrOfSetups);
sumSE_MonteCarlo_LS = zeros(length(Mrange),nbrOfSetups);


sumSE_Theoretical_MMSE = zeros(length(Mrange),nbrOfSetups);
sumSE_Theoretical_EWMMSE = zeros(length(Mrange),nbrOfSetups);
sumSE_Theoretical_LS = zeros(length(Mrange),nbrOfSetups);


sumSE_MonteCarlo_MMSE_Rayleigh = zeros(length(Mrange),nbrOfSetups);
sumSE_MonteCarlo_EWMMSE_Rayleigh = zeros(length(Mrange),nbrOfSetups);
sumSE_MonteCarlo_LS_Rayleigh = zeros(length(Mrange),nbrOfSetups);

sumSE_Theoretical_MMSE_Rayleigh = zeros(length(Mrange),nbrOfSetups);
sumSE_Theoretical_EWMMSE_Rayleigh = zeros(length(Mrange),nbrOfSetups);
sumSE_Theoretical_LS_Rayleigh = zeros(length(Mrange),nbrOfSetups);


%Preallocate matrices for storing the covariance matrices and mean vectors
RNormalized=zeros(Mmax,Mmax,K,L,L);
HMeanNormalized=zeros(Mmax,K,L,L);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %R and HMean is normalized i.e., norm(HMeanNormalized(:,k,l,j))^2=Mmax
    [RNormalized,HMeanNormalized,channelGaindB,ricianFactor,probLOS] = functionExampleSetup(L,K,Mmax,ASDdeg,0);
    
    %Controlled DL channel gain over noise
    channelGainDL= functionPowerControl( channelGaindB,noiseVariancedBm,deltadB,L);
    
    %If all UEs have same power
    % channelGainOverNoise=channelGaindB-noiseVariancedBm;
    
    for m=1:length(Mrange)
        
        %Generate the channel realizations and scale the normalized R and HMean
        [R_DL,HMean_DL,H_DL,H_DL_Rayleigh] = functionChannelGeneration( RNormalized(1:Mrange(m),1:Mrange(m),:,:,:)...
            ,HMeanNormalized(1:Mrange(m),:,:,:),channelGainDL,ricianFactor,probLOS,K,L,Mrange(m),nbrOfRealizations);
        
        %Creating zero mean vectors for Rayleigh fading i.e., blocking all the LoS paths
        HMeanZero=zeros(Mrange(m),K,L,L);
        
        %The DL SE with MMSE for Rician fading
        %Generate the MMSE channel estimates
        Hhat_MMSE_DL = functionChannelEstimateMMSE(R_DL,HMean_DL,H_DL,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute DL SE using use and then forget bound in (40) with MMSE
        SE_MonteCarlo_MMSE  = functionMonteCarloSE_DL(Hhat_MMSE_DL,H_DL,prelogFactor,nbrOfRealizations,K,L,rho);
        %Compute DL SE using Theorem 4
        SE_Theoretical_MMSE  = functionTheoreticalSE_DL_MMSE( R_DL,HMean_DL,Mrange(m),K,L,p,tau_p,prelogFactor,rho);
        
        %Save average sum SE per cell
        sumSE_MonteCarlo_MMSE(m,n) = mean(sum(SE_MonteCarlo_MMSE,1));
        sumSE_Theoretical_MMSE(m,n) = mean(sum(SE_Theoretical_MMSE,1));
        
        
        %The DL SE with MMSE for Rayleigh fading
        %Generate the MMSE channel estimates
        Hhat_MMSE_DL_Rayleigh = functionChannelEstimateMMSE(R_DL,HMeanZero,H_DL_Rayleigh,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute DL SE using use and then forget bound in (40) with MMSE
        SE_MonteCarlo_MMSE_Rayleigh = functionMonteCarloSE_DL(Hhat_MMSE_DL_Rayleigh,H_DL_Rayleigh,prelogFactor,nbrOfRealizations,K,L,rho);
        %Compute DL SE using Theorem 4
        SE_Theoretical_MMSE_Rayleigh  = functionTheoreticalSE_DL_MMSE( R_DL,HMeanZero,Mrange(m),K,L,p,tau_p,prelogFactor,rho);
        
        %Save average sum SE per cell
        sumSE_MonteCarlo_MMSE_Rayleigh(m,n) = mean(sum(SE_MonteCarlo_MMSE_Rayleigh,1));
        sumSE_Theoretical_MMSE_Rayleigh(m,n) = mean(sum(SE_Theoretical_MMSE_Rayleigh,1));
        
        
        
        
        %The DL SE with EW-MMSE for Rician fading
        %Generate the EWMMSE channel estimates
        Hhat_EWMMSE = functionChannelEstimateEWMMSE(R_DL,HMean_DL,H_DL,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute DL SE using use and then forget bound in (40) with EW-MMSE
        SE_MonteCarlo_EWMMSE = functionMonteCarloSE_DL(Hhat_EWMMSE,H_DL,prelogFactor,nbrOfRealizations,K,L,rho);
        %Compute DL SE using Theorem 5
        SE_Theoretical_EWMMSE = functionTheoreticalSE_DL_EWMMSE( R_DL,HMean_DL,Mrange(m),K,L,p,tau_p,prelogFactor,rho);
        
        
        %Save average sum SE per cell
        sumSE_MonteCarlo_EWMMSE(m,n) = mean(sum(SE_MonteCarlo_EWMMSE,1));
        sumSE_Theoretical_EWMMSE(m,n) = mean(sum(SE_Theoretical_EWMMSE,1));
        
        %The DL SE with EW-MMSE for Rayleigh fading
        %Generate the EWMMSE channel estimates
        Hhat_EWMMSE_Rayleigh = functionChannelEstimateEWMMSE(R_DL,HMeanZero,H_DL_Rayleigh,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute DL SE using use and then forget bound in (40) with EW-MMSE
        SE_MonteCarlo_EWMMSE_Rayleigh = functionMonteCarloSE_DL(Hhat_EWMMSE_Rayleigh,H_DL_Rayleigh,prelogFactor,nbrOfRealizations,K,L,rho);
        %Compute DL SE using Theorem 5
        SE_Theoretical_EWMMSE_Rayleigh = functionTheoreticalSE_DL_EWMMSE( R_DL,HMeanZero,Mrange(m),K,L,p,tau_p,prelogFactor,rho);
        
        
        %Save average sum SE per cell
        sumSE_MonteCarlo_EWMMSE_Rayleigh(m,n) = mean(sum(SE_MonteCarlo_EWMMSE_Rayleigh,1));
        sumSE_Theoretical_EWMMSE_Rayleigh(m,n) = mean(sum(SE_Theoretical_EWMMSE_Rayleigh,1));
        
        
        %The DL SE with LS for Rician fading
        %Generate the LS channel estimates
        Hhat_LS = functionChannelEstimateLS(H_DL,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute DL SE using use and then forget bound in (40) with LS
        SE_MonteCarlo_LS = functionMonteCarloSE_DL(Hhat_LS,H_DL,prelogFactor,nbrOfRealizations,K,L,rho);
        %Compute DL SE using Theorem 6
        SE_Theoretical_LS  = functionTheoreticalSE_DL_LS( R_DL,HMean_DL,Mrange(m),K,L,p,tau_p,prelogFactor,rho);
        
        %Save average sum SE per cell
        sumSE_MonteCarlo_LS(m,n) = mean(sum(SE_MonteCarlo_LS,1));
        sumSE_Theoretical_LS(m,n) = mean(sum(SE_Theoretical_LS,1));
        
        
        %The DL SE with LS for Rayleigh fading
        %Generate the LS channel estimates
        Hhat_LS_Rayleigh = functionChannelEstimateLS(H_DL_Rayleigh,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute DL SE using use and then forget bound in (40) with LS
        SE_MonteCarlo_LS_Rayleigh = functionMonteCarloSE_DL(Hhat_LS_Rayleigh,H_DL_Rayleigh,prelogFactor,nbrOfRealizations,K,L,rho);
        %Compute DL SE using Theorem 6
        SE_Theoretical_LS_Rayleigh = functionTheoreticalSE_DL_LS( R_DL,HMeanZero,Mrange(m),K,L,p,tau_p,prelogFactor,rho);
        
        %Save average sum SE per cell
        sumSE_MonteCarlo_LS_Rayleigh(m,n) = mean(sum(SE_MonteCarlo_LS_Rayleigh,1));
        sumSE_Theoretical_LS_Rayleigh(m,n) = mean(sum(SE_Theoretical_LS_Rayleigh,1));
        
        
        %Output simulation progress
        disp([num2str(Mrange(m)) ' antennas of ' num2str(Mmax)]);
    end
    
    
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    clear RNormalized HMeanNormalized
    
end

%% Plot the figure
figure;
c1=plot(Mrange,mean(sumSE_MonteCarlo_MMSE,2),'ks','MarkerSize',7);
hold on
m1=plot(Mrange,mean(sumSE_Theoretical_MMSE,2),'r-+','LineWidth',1);
hold on
c2=plot(Mrange,mean(sumSE_MonteCarlo_EWMMSE,2),'ks','MarkerSize',7);
hold on
m2=plot(Mrange,mean(sumSE_Theoretical_EWMMSE,2),'b-x','LineWidth',1);
hold on
c3=plot(Mrange,mean(sumSE_MonteCarlo_LS,2),'ks','MarkerSize',7);
hold on
m3=plot(Mrange,mean(sumSE_Theoretical_LS,2),'k.-','LineWidth',1);
hold on
c1r=plot(Mrange,mean(sumSE_MonteCarlo_MMSE_Rayleigh,2),'ks','MarkerSize',7);
hold on
m1r=plot(Mrange,mean(sumSE_Theoretical_MMSE_Rayleigh,2),'r+:','LineWidth',1);
hold on
c2r=plot(Mrange,mean(sumSE_MonteCarlo_EWMMSE_Rayleigh,2),'ks','MarkerSize',7);
hold on
m2r=plot(Mrange,mean(sumSE_Theoretical_EWMMSE_Rayleigh,2),'bx:','LineWidth',1);
hold on
c3r=plot(Mrange,mean(sumSE_MonteCarlo_LS_Rayleigh,2),'ks','MarkerSize',7);
hold on
m3r=plot(Mrange,mean(sumSE_Theoretical_LS_Rayleigh,2),'k.:','LineWidth',1);



xlabel('Number of antennas (M)');
ylabel('Average sum SE [bit/s/Hz/cell]');
legend([m1 m2 m3 m1r m3r ],'MMSE with Rician fading',' EW-MMSE with Rician fading'...
    ,'LS with Rician fading','MMSE with Rayleigh fading', 'EW-MMSE/LS with Rayleigh','Location','SouthEast');




