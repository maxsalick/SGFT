function [m_full_cov,im_noedge] = sft_cov_scan(im,scanjump,blocksize,umperpix,patestimate)

% This function scans the image and determines regions that have enough
% positive signal to be included in the analysis.

% The output, m_full_cov, is an array with 1s indicating regions with
% adequate signal, and 0s indicating regions that will be skipped to reduce
% computation time of Fourier analysis.

% v2 - changed to include patestimate./umperpix for dilation to ensure that
% dilated positive region will span negative regions

disp('Beginning preliminary COVERAGE scan loop...')
disp(' ')

disp('First scan loop will determine regions that contain cells.')
disp('Edges of positive regions will be reduced (within 0.75 blocksizes).')
disp('This reduces the false positives that result from edge effects during')
disp('Fast Fourier Transform image analysis.')



imagesize = size(im);

meanintense = ceil(mean(mean(im(im>0))));

dilated = imdilate(im,strel('disk',floor(patestimate./umperpix)));
meanintense = ceil(mean(mean(dilated(dilated>0))));
BW = roicolor(dilated,ceil(meanintense/2),255);
BW2 = imfill(BW,'holes');
roierode = imerode(BW2,strel('disk',floor(patestimate./umperpix)));
m_full_cov = roierode;

disp(' ')
disp('Positive regions of image detected.')
disp('Removing edges of image...')

D = bwdist(1-roierode);
D(D>blocksize*.2) = blocksize*.2;
D = (D/(blocksize*.2)).^2;


im_noedge = double(im).*D;


disp(' ')
disp('First scan loop (COVERAGE) completed!')

end

