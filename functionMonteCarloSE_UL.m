function [SE_MonteCarlo] = functionMonteCarloSE_UL(Hhat,H,prelogFactor,nbrOfRealizations,M,K,L,p)
%Computes the uplink SE by Monte Carlo simulations.
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


%Prepare to store simulation results
SE_MonteCarlo = zeros(K,L);

%Store the terms in SINR (20)
Term1x=zeros(K,L,nbrOfRealizations);
Term2x=zeros(K,L,nbrOfRealizations);
Term3x=zeros(K,L,nbrOfRealizations);

% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel realizations from all UEs to BS j
        Hallj = reshape(H(:,n,:,:,j),[M K*L]);
        
        
        %Go through all UEs in cell j
        for k = 1:K
            
            %Compute the terms in SINR (20)
            Term1x(k,j,n) =Hhat(:,n,k,j,j)'*H(:,n,k,j,j);
            Term2x(k,j,n)= sum(abs(Hhat(:,n,k,j,j)'*Hallj).^2) ;
            Term3x(k,j,n)= norm(Hhat(:,n,k,j,j))^2;
            
        end
        
    end
    
end

%Calculate the averages
Term1=mean(Term1x,3);
Term2=mean(Term2x,3);
Term3=mean(Term3x,3);

%Compute the SE in (19)for all UEs
for k=1:K
    for j=1:L
        SE_MonteCarlo(k,j)= prelogFactor*real(log2(1+ (p*abs(Term1(k,j))^2)/(p*Term2(k,j) - p*abs(Term1(k,j))^2 +    Term3(k,j) )    ))  ;
    end
end


end
