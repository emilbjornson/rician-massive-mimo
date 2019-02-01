function R = functionRlocalscatteringApprox(M,theta,ASDdeg,antennaSpacing,N)
%Generate the spatial covariance matrix with multiple scattering clusters,
%each described by the local scattering model with Gaussian angular
%distribution.
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


%Set the antenna spacing if not specified by input
if  nargin < 4
    
    %Half a wavelength distance
    antennaSpacing = 1/2;
    
end


%Compute the ASD in radians based on input
ASD = ASDdeg*pi/180;


%The correlation matrix has a Toeplitz structure, so we only need to
%compute the first row of the matrix
firstRow = zeros(M,1);

%Go through all the N paths


for n=1:N
    
    thetaN=theta + (-0.691+1.3963*rand(1));  %U[-40,40] deviation around the nominal angle (in radians)
    
    %Go through all the columns of the first row
    for column = 1:M
        
        %Distance from the first antenna
        distance = column-1;
        
        
        %Compute the approximated integral as in (2.24)
        firstRow(column,n) = exp(1i*2*pi*antennaSpacing*sin(thetaN)*distance)*exp(-((ASD^2)/2) * ( 2*pi*antennaSpacing*cos(thetaN)*distance )^2);
        
    end
end
Rx=sum(firstRow,2)/N;

%Compute the spatial correlation matrix by utilizing the Toeplitz structure
R = toeplitz(Rx);
