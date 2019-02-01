function SE_Theoretical = functionTheoreticalSE_UL_LS( R,HMean,M,K,L,p,tau_p,prelogFactor)
%Computes the UL SE with LS estimator in Theorem 3
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

%Prepare to store (31),(32) and (32)
Term1=zeros(K,L); %(31)
Term2=zeros(K,L); %(33)
Term2_p1=zeros(K,L,L,K);
Term2_p2=zeros(K,L,L,K);
Term1_p1=zeros(K,L,L,K);
Term3=zeros(K,L); %(32)

%Prepare to store UL SE
SE_Theoretical=zeros(K,L);


% Go through all BSs
for j = 1:L
    
    %Go through all UEs
    for k = 1:K
        
        PsiInv = (p*tau_p*sum(R(:,:,k,:,j),4) + eyeM);
        yMean=sqrt(p)*tau_p*sum(HMean(:,k,:,j),3);
        %Compute (32)  in Theorem 3
        Term3(k,j)=(1/(p*tau_p))*trace(PsiInv)+(1/(p*tau_p*tau_p))*norm(yMean)^2;
        
        %Go through all UEs
        for lx=1:L
            for i=1:K
                if k==i
                    
                    %Compute the second term in (31)
                    Term1_p1(k,j,lx,i)=HMean(:,i,lx,j)'*HMean(:,k,j,j);
                    %Calculate x_bar and Omega
                    x=yMean-sqrt(p)*tau_p*HMean(:,i,lx,j);
                    Omega=PsiInv-p*tau_p*R(:,:,i,lx,j);
                    
                    %Calculate (33) in Theorem 3 for pilot contaminated UEs
                    Term2_p1(k,j,lx,i)=(1/(p*tau_p*tau_p))*(tau_p*(trace(R(:,:,i,lx,j)*PsiInv)) ...
                        + abs(x'*HMean(:,i,lx,j))^2 +  p*tau_p*tau_p*norm(HMean(:,i,lx,j))^4+ p*tau_p*tau_p*abs(trace(R(:,:,i,lx,j)))^2 ...
                        + x'*R(:,:,i,lx,j)*x + tau_p*HMean(:,i,lx,j)'*Omega*HMean(:,i,lx,j) ....
                        + 2*sqrt(p)*tau_p*real(yMean'*HMean(:,i,lx,j)*trace(R(:,:,i,lx,j)) ...
                        + yMean'*R(:,:,i,lx,j)*HMean(:,i,lx,j) +x'*HMean(:,i,lx,j)*norm(HMean(:,i,lx,j))^2 ));
                    
                else
                    %Calculate (33) in Theorem 3 for non-pilot contaminated UEs
                    Term2_p2(k,j,lx,i)=(1/(p*tau_p*tau_p))*(tau_p*trace(R(:,:,i,lx,j)*PsiInv) ...
                        + abs(yMean'*HMean(:,i,lx,j))^2 ...
                        + yMean'*R(:,:,i,lx,j)*yMean+ tau_p*HMean(:,i,lx,j)'*PsiInv*HMean(:,i,lx,j) );
                end
            end
        end
        
        %Compute (31)
        Term1(k,j)=trace(R(:,:,k,j,j))+ sum(sum(Term1_p1(k,j,:,:),3),4);
        %Sum all the interferences from all UEs
        Term2(k,j)=sum(sum(Term2_p1(k,j,:,:),3),4)+ sum(sum(Term2_p2(k,j,:,:,:),3),4);
    end
end

%Compute the SE with (34) for each UE
for k=1:K
    for j=1:L
        SE_Theoretical(k,j)= prelogFactor*real(log2(1+ (p*abs(Term1(k,j))^2)/(p*Term2(k,j) -p*abs(Term1(k,j))^2 +Term3(k,j))));
    end
end


end
