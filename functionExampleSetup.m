function [R,HMean,channelGaindB,ricianFactor,probLOS] = functionExampleSetup(L,K,M,ASDdeg,allLoS)
%This function generates the channel statistics between UEs at random
%locations.
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


%Model parameters
%Set the length in meters of the total square area
squareLength = 1000;

%Number of BSs per dimension
nbrBSsPerDim = sqrt(L);

%Standard deviation of shadow fading in dB
sigma_sf_NLOS=10; %for NLOS
sigma_sf_LOS=4;   %for LOS
%Minimum distance between BSs and UEs
minDistance = 35;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance
%Number of the paths
N=6;
%Distance between BSs in vertical/horizontal direction
interBSDistance = squareLength/nbrBSsPerDim;

%Deploy BSs on the grid
locationsGridHorizontal = repmat(interBSDistance/2:interBSDistance:squareLength-interBSDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);

%Compute all nine alternatives of the BS locations when using wrap around
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

%Prepare to put out UEs in the cells
UEpositions = zeros(K,L);
perBS = zeros(L,1);

%Prepare to store normalized spatial covariance matrices and mean vectors
R = zeros(M,M,K,L,L);
HMean=zeros(M,K,L,L);

%Prepare to store channel gain (in dB), Rician Factor \kappa and
%probabaility of LoS for each UE
channelGaindB=zeros(K,L,L);
probLOS=zeros(K,L,L);
ricianFactor=zeros(K,L,L);

%The maximum distance for LoS
maxdistLOS=300;

% Go through all the cells
for l = 1:L
    
    %Put out K UEs in the cell, uniformly at random. The procedure is
    %iterative since UEs that do not satisfy the minimum distance are
    %replaced with new UEs
    while perBS(l)<K
        
        %Put out new UEs
        UEremaining = K-perBS(l);
        posX = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
        posY = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
        posXY = posX + 1i*posY;
        
        %Keep those that satisfy the minimum distance
        posXY = posXY(abs(posXY)>=minDistance);
        
        %Store new UEs
        UEpositions(perBS(l)+1:perBS(l)+length(posXY),l) = posXY + BSpositions(l);
        perBS(l) = perBS(l)+length(posXY);
        
    end
    
    
    %Go through all BSs
    for j = 1:L
        
        %Compute the distance from the UEs in cell l to BS j with a wrap
        %around topology, where the shortest distance between a UE and the
        %nine different locations of a BS is considered
        
        [distancesBSj,~] = min(abs( repmat(UEpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);
        
        %Compute the probability of LoS and Rician Factor (in linear scale)
        probLOS(:,l,j)=(rand<((maxdistLOS-distancesBSj)./maxdistLOS));
        ricianFactor(:,l,j)=db2pow(13-0.03*distancesBSj);
        %If all UEs on the network have LoS paths
        if allLoS==1
            probLOS=ones(K,L,L);
        end
        
        %Compute average channel gain using the large-scale fading
        %model based on 3GPP model while neglecting the shadow fading for each UE
        for k = 1:K
            
            if probLOS(k,l,j)==1
                channelGaindB(k,l,j)= -30.18-26*log10(distancesBSj(k));
            else
                channelGaindB(k,l,j)= -34.53-38*log10(distancesBSj(k));
            end
        end
        
    end
    
    
    %Go through all UEs in cell l and generate shadow fading realizations
    for k = 1:K
        
        if probLOS(k,l,l)==1
            shadowing = sigma_sf_LOS*randn(1,1,L);
            channelGainShadowing = channelGaindB(k,l,:) + shadowing;
        else
            shadowing = sigma_sf_NLOS*randn(1,1,L);
            channelGainShadowing = channelGaindB(k,l,:) + shadowing;
        end
        
        %Check if another BS has a larger average channel gain to the UE
        %than BS l
        while channelGainShadowing(l) < max(channelGainShadowing)
            
            if probLOS(k,l,l)==1
                shadowing = sigma_sf_LOS*randn(1,1,L);
                channelGainShadowing = channelGaindB(k,l,:) + shadowing;
            else
                shadowing = sigma_sf_NLOS*randn(1,1,L);
                channelGainShadowing = channelGaindB(k,l,:) + shadowing;
            end
            
        end
        
        %Store average channel gains with shadowing fading
        channelGaindB(k,l,:) = channelGainShadowing;
    end
    
    
    for j=1:L
        
        [~,whichpos] = min(abs( repmat(UEpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);
        
        for k=1:K
            %Compute nominal angles between UE k in cell l and BS j, and
            %generate spatial correlation matrices for the channels using the
            %local scattering model
            
            angleBSj = angle(UEpositions(k,l)-BSpositionsWrapped(j,whichpos(k)));
            
            
            HMean(:,k,l,j)=(exp(1i*2*pi.*(0:(M-1))*sin(angleBSj)*antennaSpacing)); %Normalized Mean vector
            R(:,:,k,l,j) = functionRlocalscatteringApprox(M,angleBSj,ASDdeg,antennaSpacing,N); %Normalized matrix
        end
        
    end
    
end
