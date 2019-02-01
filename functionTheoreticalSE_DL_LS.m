function SE_Theoretical = functionTheoreticalSE_DL_LS( R,HMean,M,K,L,p,tau_p,prelogFactor,rho)
%Computes the DL SE with LS estimator in Theorem 6
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

%Prepare to store  (50) and (51)
Term1=zeros(K,L);
Term1_p1=zeros(K,L,L,K);
Term2=zeros(K,L);
Term2_p1=zeros(K,L,L,K);
Term2_p2=zeros(K,L,L,K);
Term3=zeros(K,L);

%Prepare to store DL SE
SE_Theoretical=zeros(K,L);


% Go through all BSs
for j = 1:L
    
    %Go through all UEs
    for k = 1:K
        %Compute the mean norm square of the estimator
        PsiInv = (p*tau_p*sum(R(:,:,k,:,j),4) + eyeM);
        yMean=sqrt(p)*tau_p*sum(HMean(:,k,:,j),3);
        Term3(k,j)=(1/(p*tau_p))*trace(PsiInv)+(1/(p*tau_p*tau_p))*norm(yMean)^2;
        
        for lx=1:L
            for i=1:K
                
                %Compute the mean norm square of the estimator for all UEs
                PsiInv_li = (p*tau_p*sum(R(:,:,i,:,lx),4) + eyeM);
                yMean_li=sqrt(p)*tau_p*sum(HMean(:,i,:,lx),3);
                mu_li=(1/(p*tau_p))*trace(PsiInv_li)+(1/(p*tau_p*tau_p))*norm(yMean_li)^2;
                
                
                
                if k==i
                    %Compute the second term in the numerator of (50)
                    Term1_p1(k,j,lx,i)=HMean(:,i,lx,j)'*HMean(:,k,j,j);
                    
                    %Calculate x_bar and Omega
                    x=yMean_li-sqrt(p)*tau_p*HMean(:,k,j,lx); %Actually xMean
                    Omega=PsiInv_li-p*tau_p*R(:,:,k,j,lx);
                    
                    %Calculate (51) in Theorem 6 for pilot contaminated UEs
                    Term2_p1(k,j,lx,i)=((1/(p*tau_p*tau_p))*(tau_p*abs(trace(R(:,:,k,j,lx)*PsiInv_li)) + abs(x'*HMean(:,k,j,lx))^2 ...
                        + p*tau_p*tau_p*norm(HMean(:,k,j,lx))^4 + x'*R(:,:,k,j,lx)*x + tau_p*HMean(:,k,j,lx)'*Omega*HMean(:,k,j,lx) ...
                        + 2*sqrt(p)*tau_p*real(yMean'*HMean(:,k,j,lx)*trace(R(:,:,k,j,lx)) ...
                        + yMean'*R(:,:,k,j,lx)*HMean(:,k,j,lx) +x'*HMean(:,k,j,lx)*HMean(:,k,j,lx)'*HMean(:,k,j,lx))  )...
                        + abs(trace(R(:,:,k,j,lx)))^2 )/mu_li  ;
                    
                else
                    %Calculate (51) in Theorem 6 for non-pilot contaminated UEs
                    Term2_p2(k,j,lx,i)=((1/(p*tau_p*tau_p))*(tau_p*abs(trace(R(:,:,k,j,lx)*PsiInv_li))...
                        +yMean_li'*R(:,:,k,j,lx)*yMean_li+ tau_p*HMean(:,k,j,lx)'*PsiInv_li*HMean(:,k,j,lx) ...
                        + abs(yMean_li'*HMean(:,k,j,lx))^2 ))/mu_li;
                end
            end
        end
        
        %Compute the numerator of (50)
        Term1(k,j)=trace(R(:,:,k,j,j))+ sum(sum(Term1_p1(k,j,:,:),3),4);
        %Sum all the interferences from all UEs
        Term2(k,j)=sum(sum(Term2_p1(k,j,:,:),3),4)+ sum(sum(Term2_p2(k,j,:,:,:),3),4);
    end
end

%Compute the SE with (52) for each UE
for k=1:K
    for j=1:L
        
        SE_Theoretical(k,j)= prelogFactor*real(log2(1+ (((rho*abs(Term1(k,j))^2)/(Term3(k,j)))/(rho*Term2(k,j) -((rho*abs(Term1(k,j))^2)/(Term3(k,j))) + 1))));
    end
end


end
