function [imagesize, m_full_str, m_full_dir, m_full_sl, quiver] = ...
    sft_scan_looper(imvert, imhor, m_full_cov, im,  ... 
                    blocksize, scanjump, umperpix, overridebin, overridedir)

%   This file coordinates "sft_scan.m", giving the scan locations and
%   spacings, and constructing output arrays based on the sft_scan values.

%  Scans through the entire image to detect patterning.  First conducts
%  gradient analysis to determine directionality of the image.  Next, 1-D
%  Fourier transforms are taken in the determined directions, and peak
%  detection of averaged Fourier transforms is conducted to assess the
%  amount of organization at the scanned areas.  Spacing between patterns
%  is determined spatially using intensity profiles across the image in the
%  known direction of patterning.
                
                
disp('Beginning second scan loop.  This will utilize gradient methods')
disp('to determine the orientation of patterns within the image, and')
disp('will then apply 1D Fourier transforms in this direction to determine')
disp('pattern strength.')
disp(' ')
     

disp('Initializing main FFT scan loop...')
disp(' ')

fftw('planner','exhaustive')

[GX,GY] = gradient(im);
% [GX,GY] = gradient(max(im-median(im(im~=0)),0));
orientim = atan2(-GY,GX);
gradstrbig = (GY.^2+GX.^2).^.5;
if overridebin == 1
    orientim = overridedir.*(pi/180).*ones(size(orientim));
    gradstrbig = ones(size(orientim));
end
imagesize = [imvert imhor];
m_full_str = zeros(imagesize);
m_full_dir = zeros(imagesize);
m_full_sl = zeros(imagesize);
quiver = zeros(imagesize);
totaltoscan = sum(sum(m_full_cov));

tic
scanned=0;

mstart = zeros(scanjump);
mstart(floor(scanjump)/2,1:end)=1;

for xloc = ceil(blocksize/2)+1:scanjump:imagesize(2)-(ceil(blocksize/2)+1)
    
    timer = toc;
    timeremain = timer*totaltoscan/scanned - timer;
    minremain = floor(timeremain/60);
    secremain = floor(timeremain - 60*minremain);
    disp([num2str(ceil(100*scanned/totaltoscan)) '% complete. ' ...
        num2str(minremain) ':' num2str(secremain) ' remaining.'])

    for yloc = ceil(blocksize/2)+1:scanjump:imagesize(1)-(ceil(blocksize/2)+1)
        
      if m_full_cov(yloc, xloc) == 1
         scanned = scanned+scanjump^2;
         section = im(yloc-blocksize/2:yloc+blocksize/2-1,xloc-blocksize/2:xloc+blocksize/2-1);
         orientsection = orientim(yloc-blocksize/2:yloc+blocksize/2-1,xloc-blocksize/2:xloc+blocksize/2-1);
         gradstr = gradstrbig(yloc-blocksize/2:yloc+blocksize/2-1,xloc-blocksize/2:xloc+blocksize/2-1);

         section = section-min(min(section)).*255/max(max(1+section-min(min(section))));
        
 [fourstrengthval,sarclength,meandirect] = sft_scan3(section,orientsection,gradstr,blocksize,umperpix);


m_full_str(yloc-scanjump/2:yloc+scanjump/2,xloc-scanjump/2:xloc+scanjump/2)=fourstrengthval;
m_full_dir(yloc-scanjump/2:yloc+scanjump/2,xloc-scanjump/2:xloc+scanjump/2)=meandirect;
m_full_sl(yloc-scanjump/2:yloc+scanjump/2,xloc-scanjump/2:xloc+scanjump/2)=sarclength;



m = imrotate(mstart,double(meandirect*180/pi),'bilinear','crop');
quiver(yloc-scanjump/2:-1+yloc+scanjump/2,xloc-scanjump/2:-1+xloc+scanjump/2)=m;

        end
    end
   
end
                
      
disp(' ')
disp('Analysis complete!')        
                
       
                
end