function im = sft_imfilter(imagedataorig,blocksize)

% Applies local contrast filter to image, and normalizes results.
disp(' ')
disp('Applying local contrast filter to image.')
disp('This involves producing a Gaussian blur of the image.')
disp('It may take a few minutes for large images.')

im_b = 255-imagedataorig;
filterbox = floor(blocksize/12);

close all
MD =imfilter(im_b,fspecial('disk',filterbox),'same');
im = im_b - MD;

% Normalize local contrasted image
disp('Normalizing filtered image.')
im = (im-min(min(im))).*(255/double(max(max(im-min(min(im))))));




end