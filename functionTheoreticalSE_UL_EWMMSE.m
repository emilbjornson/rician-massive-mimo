function SE_Theoretical = functionTheoreticalSE_UL_EWMMSE( R,HMean,M,K,L,p,tau_p,prelogFactor)   
%Computes the UL SE with EW-MMSE estimator in Theorem 2
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

%Prepare to store  (27),(28) and (29)
CCterm1=zeros(K,L); %(27)
CCterm2=zeros(K,L); %(29)
CCterm2_p1=zeros(K,L,L,K);
CCterm3=zeros(K,L); %(28)

%Prepare to store UL SE
SE_Theoretical=zeros(K,L);


% Go through all BSs
for j = 1:L
    
    %Go through all UEs
    for k = 1:K
       
        %Compute the matrix that is inverted in the MMSE estimator
        PsiInv = (p*tau_p*sum(R(:,:,k,:,j),4) + eyeM);
        %Compute the matrix that is inverted in the EWMMSE
        LambdaInv=diag(diag(PsiInv));
        %Taking the diagonal elements of the covariance matrices
        Djk=diag(diag(R(:,:,k,j,j)));
        %Compute (27) in Theorem 2
        CCterm1(k,j)=p*tau_p*trace(Djk/LambdaInv*Djk) + norm(HMean(:,k,j,j))^2;
        %Compute (28) in Theorem 2
        Sigmajk=p*tau_p*Djk/LambdaInv*PsiInv/LambdaInv*Djk;
        CCterm3(k,j)=trace(Sigmajk)+ norm(HMean(:,k,j,j))^2;
        %Go through all UEs
        for lx=1:L
            for i=1:K
                %Calculate (29) in Theorem 1
                %Compute the common terms for non-pilot contaminated and pilot contaminated UEs
                CCterm2_p1(k,j,lx,i)=trace(R(:,:,i,lx,j)*Sigmajk) +...
                    HMean(:,i,lx,j)'*Sigmajk*HMean(:,i,lx,j)+...
                    HMean(:,k,j,j)'*R(:,:,i,lx,j)*HMean(:,k,j,j)+...
                    abs(HMean(:,k,j,j)'*HMean(:,i,lx,j))^2;
                if k==i  %Compute the additional terms for pilot contaminated UEs
                    CCterm2_p1(k,j,lx,i)= CCterm2_p1(k,j,lx,i)+ p*p*tau_p*tau_p*abs(trace(R(:,:,i,lx,j)/LambdaInv*Djk))^2+...
                        2*p*tau_p*real(trace(R(:,:,i,lx,j)/LambdaInv*Djk)*HMean(:,k,j,j)'*HMean(:,i,lx,j));
                end
                
            end
        end
        
        %Sum all the interferences from all UEs
        CCterm2(k,j)=sum(sum(CCterm2_p1(k,j,:,:),3),4);
    end
end

    
    
    %Compute the SE with (30) for each UE
    for k=1:K
        for j=1:L
            SE_Theoretical(k,j)=prelogFactor*real(log2(1+ (p*abs(CCterm1(k,j))^2)/(p*CCterm2(k,j) -p*abs(CCterm1(k,j))^2 + CCterm3(k,j))));
        end
    end

end

