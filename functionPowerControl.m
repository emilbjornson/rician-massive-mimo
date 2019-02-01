function [ channelGainOverNoise ] = functionPowerControl( channelGaindB,noiseVariancedBm,deltadB,L)
%Uplink power control part based on [1] Section 7, Eq. 7.11
%
%[1] Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017),
%"Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency",
%Foundations and Trends in Signal Processing: Vol. 11, No. 3-4,
%pp. 154-655. DOI: 10.1561/2000000093.
%
%This Matlab function was developed to generate simulation results to:
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


%Compute the normalized average channel gain, where the normalization
%is based on the noise power
channelGainOverNoiseOriginal = channelGaindB-noiseVariancedBm;
channelGainOverNoise = channelGainOverNoiseOriginal;

%Go through all cells
for j = 1:L
    %Scale the average channel gains by applying the power control
    %policy of (7.11). Note that we are including this power
    %control here so that we can then view it as if all UEs of
    %transmitting at maximum power
    
    betajMin = min(channelGainOverNoiseOriginal(:,j,j));
    differenceSNR = channelGainOverNoiseOriginal(:,j,j)-betajMin;
    backoff = differenceSNR-deltadB;
    backoff(backoff<0) = 0;
    channelGainOverNoise(:,j,:) = channelGainOverNoiseOriginal(:,j,:)-repmat(backoff,[1 1 L]);
    
end


end

