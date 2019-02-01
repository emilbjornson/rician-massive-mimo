function [Hhat_MMSE] = functionChannelEstimateMMSE(R,HMeanx,H,nbrOfRealizations,M,K,L,p,f,tau_p)
%Generating MMSE channel estimates
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


%Reshape the mean vectors to get same mean vector for all realizations
HMean=reshape(repmat(HMeanx,nbrOfRealizations,1),M,nbrOfRealizations,K,L,L);



%Generate pilot pattern
if f == 1
    
    pilotPattern = ones(L,1);
    
elseif f == 2 %Only works in the running example with its 16 BSs
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 2; 1; 2; 1]);
    
elseif f == 4 %Only works in the running example with its 16 BSs
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 3; 4; 3; 4]);
    
elseif f == 16 %Only works in the running example with its 16 BSs
    
    pilotPattern = (1:L)';
    
end


%Store identity matrix of size M x M
eyeM = eye(M);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(M,nbrOfRealizations,K,L,f) + 1i*randn(M,nbrOfRealizations,K,L,f));


%Prepare to store MMSE channel estimates
Hhat_MMSE = zeros(M,nbrOfRealizations,K,L,L);

% Go through all cells
for j = 1:L
    
    %Go through all f pilot groups
    for g = 1:f
        
        %Extract the cells that belong to pilot group g
        groupMembers = find(g==pilotPattern)';
        
        %Compute processed pilot signal for all UEs that use these pilots, according to (5)
        yp = sqrt(p)*tau_p*sum(H(:,:,:,g==pilotPattern,j),4) + sqrt(tau_p)*Np(:,:,:,j,g);
        yMean=sqrt(p)*tau_p*sum(HMean(:,:,:,g==pilotPattern,j),4);
        %Go through all UEs
        for k = 1:K
            
            %Compute the matrix that is inverted in the MMSE estimator
            PsiInv = (p*tau_p*sum(R(:,:,k,g==pilotPattern,j),4) + eyeM);
            
            
            %Go through the cells in pilot group g
            for l = groupMembers
                
                %Compute MMSE estimate of channel between BS l and UE k in
                %cell j using (6) in Lemma 1
                RPsi = R(:,:,k,l,j) / PsiInv;
                Hhat_MMSE(:,:,k,l,j) = HMean(:,:,k,l,j) + sqrt(p)*RPsi*(yp(:,:,k)-yMean(:,:,k));
                
                
            end
            
        end
        
    end
    
end

