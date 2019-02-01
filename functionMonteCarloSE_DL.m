function [SE_MonteCarlo] = functionMonteCarloSE_DL(Hhat,H,prelogFactor,nbrOfRealizations,K,L,rho)
%Computes the downlink SE by Monte Carlo simulations.
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
SE_MonteCarlo=zeros(K,L);


%Prepare to store terms in (41)
Denominatorx=zeros(K,L,nbrOfRealizations);
Term1x=zeros(K,L,nbrOfRealizations);
Term2x=zeros(K,L,L,K,nbrOfRealizations);

% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Go through all UEs in cell j
        for k = 1:K
            %Compute denominator of the estimator
            Denominatorx(k,j,n) = norm(Hhat(:,n,k,j,j))^2;
        end
        
    end
    
end
%The denominator of the estimator
Denominator=sqrt(mean(Denominatorx,3));


for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Go through all UEs in cell j
        for k = 1:K
            
            %Compute the first term  in (41)
            Term1x(k,j,n) =(Hhat(:,n,k,j,j)/Denominator(k,j))'*H(:,n,k,j,j);
            
            for l=1:L
                for i=1:K
                    %Compute the second term in SINR (41)
                    Term2x(k,j,l,i,n)= abs((Hhat(:,n,i,l,l)/Denominator(i,l))'*H(:,n,k,j,l)).^2 ;
                end
            end
            
        end
        
    end
    
end

%Calculate the averages
Term1=abs(mean(Term1x,3)).^2;
Term2=mean(Term2x,5);

%Compute the SE in (40)for all UEs
for k=1:K
    for j=1:L
        
        Term2(k,j)=(sum(sum(Term2(k,j,:,:),3),4));
        SE_MonteCarlo(k,j)= prelogFactor*log2( 1 + (rho*Term1(k,j))/(rho*Term2(k,j)- rho*Term1(k,j)+ 1)) ;
    end
end



end
