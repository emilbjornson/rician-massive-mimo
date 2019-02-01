function SE_Theoretical = functionTheoreticalSE_UL_MMSE( R,HMean,M,K,L,p,tau_p,prelogFactor)
%Computes the UL SE with MMSE estimator in Theorem 1
%Note that the covariance matrices and the mean vectors are same for all
%channel realizations.
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


%Store identity matrix of size M x M
eyeM = eye(M);

%Prepare to store  (21) and (23)
CCterm1=zeros(K,L); %(21)
CCterm2=zeros(K,L); %(23)
CCterm2_p1=zeros(K,L,L,K);

%Prepare to store UL SE
SE_Theoretical=zeros(K,L);


% Go through all BSs
for j = 1:L
    
    %Go through all UEs
    for k = 1:K
        
        %Compute the matrix that is inverted in the MMSE estimator
        PsiInv = (p*tau_p*sum(R(:,:,k,:,j),4) + eyeM);
        %Compute (21) & (22) in Theorem 1
        CCterm1(k,j)=p*tau_p*trace(R(:,:,k,j,j)/PsiInv*R(:,:,k,j,j)) + norm(HMean(:,k,j,j))^2;
        
        %Go through all UEs
        for lx=1:L
            for i=1:K
                %Calculate (23) in Theorem 1
                %Compute the common terms for non-pilot contaminated and pilot contaminated UEs
                CCterm2_p1(k,j,lx,i)=p*tau_p*(trace(R(:,:,i,lx,j)*R(:,:,k,j,j)/PsiInv*R(:,:,k,j,j))) +...
                    (p*tau_p*HMean(:,i,lx,j)'*R(:,:,k,j,j)/PsiInv*R(:,:,k,j,j)*HMean(:,i,lx,j))+...
                    ( HMean(:,k,j,j)'*R(:,:,i,lx,j)*HMean(:,k,j,j))+...
                    abs(HMean(:,k,j,j)'*HMean(:,i,lx,j))^2;
                if k==i %Compute the additional terms for pilot contaminated UEs
                    CCterm2_p1(k,j,lx,i)= CCterm2_p1(k,j,lx,i)+ p*p*tau_p*tau_p*abs(trace(R(:,:,i,lx,j)/PsiInv*R(:,:,k,j,j)))^2+...
                        2*p*tau_p*real(trace(R(:,:,i,lx,j)/PsiInv*R(:,:,k,j,j))*HMean(:,k,j,j)'*HMean(:,i,lx,j));
                end
                
            end
        end
        
        %Sum all the interferences from all UEs
        CCterm2(k,j)=sum(sum(CCterm2_p1(k,j,:,:),3),4);
    end
end



%Compute the SE with (24) for each UE
for k=1:K
    for j=1:L
        SE_Theoretical(k,j)=prelogFactor*real(log2(1+ (p*abs(CCterm1(k,j))^2)/(p*CCterm2(k,j) -p*abs(CCterm1(k,j))^2 + CCterm1(k,j))));
    end
end

end
