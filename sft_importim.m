function [imagedata, imagesize, imhor, imvert] = ...
   sft_importim(PathName, FileName)
% Imports image and its parameters, given a path and file name.



    imagedata = imread([PathName FileName]);
    
    
    imagedata = max(max(imagedata(:,:,1),imagedata(:,:,2)),imagedata(:,:,3));
    
%     --------------------------------------------
% USER - for co-stain analysis, you may select out the R, G, or B channels.
%  See below for a quick example...

% imagedata = imagedata(:,:,1)   % this will only use the R channel for analysis
% imagedata = imagedata(:,:,2)   % this will only use the G channel for analysis
% imagedata = imagedata(:,:,3)   % this will only use the B channel for analysis

% Make sure to comment out line 10 if you do this.
%     --------------------------------------------
    


    if mean(mean(imagedata))<127
        imagedata = 255-imagedata;
    end
imagemax = max(max(imagedata));
imagecorrect= 255./imagemax;
imagedata = imagedata.*imagecorrect;
imagesize = size(imagedata);
imhor = imagesize(2);
imvert = imagesize(1);


end

