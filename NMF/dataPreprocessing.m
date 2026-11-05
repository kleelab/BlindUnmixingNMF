function [Y,Y_sps, H_sps] = dataPreprocessing(A0,Y,Thrs)
% dataPreprocessing  Sparse pixel/timepoint selection for NMF-RI.
%
%   [Y, Y_sps, H_sps] = dataPreprocessing(A0, Y)
%   [Y, Y_sps, H_sps] = dataPreprocessing(A0, Y, Thrs)
%
%   Thrs  Sparseness threshold (default 0.7 = select top 30%).
%         Higher Thrs -> fewer but sparser samples.
%         For photometry, try 0.5~0.7 depending on signal overlap.
%
% Based on NMF-RI (Jimenez-Sanchez et al., Bioinformatics 2020)
% Synthetic experiments version, with optional Thrs parameter added.

if nargin < 3
    Thrs = 0.7;   % Original default: select 30% of pixels.
end

% Elimination of dark current.
Y = Y-min(Y);

% Matrix Initialization
H0 = max(100*eps,pinv([A0])*[Y]); 

% Sparseness selection method.
H_sps = H0;
Y_sps = Y;
% segmentation to obtain sparsest pixels
SPS = ( sqrt(size(H_sps,1))-(sum(abs(H_sps),1)./sqrt(sum(H_sps.^2,1))) )/(sqrt(size(H_sps,1))-1);                
% Cumulative distribution function to find sparsest pixels.
[f, x] = ecdf(SPS);
thresh = x(round(length(x)*Thrs));
Y_sps=Y(:,SPS>thresh);  
% Prepare H0 to unmix.-
H_sps = max(100*eps,pinv(A0)*Y_sps); 

end
