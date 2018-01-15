function imagedataorig = sft_remscalebar(imagedata,imhor,imvert,vertres)

% Manually remove scalebar from image.  The targetted region will be
% blurred to prevent it from being analyzed and giving false results.  Also
% can be used on debris in the image that is not representative.

    disp(' ');
disp('Scalebar must now be removed to keep detection software from falsely including it in the data');
disp('Please outline the scalebar, then right click and select "FILL AREA" ');
    disp(' ')
    disp('PRESS ANY KEY TO BEGIN')
    close all
    pause
figure('position',[50,50,floor((vertres-100)*imhor/imvert),vertres-100]);
    colormap('gray')
    imagesc(imagedata);
    
    title('Select scalebar/debris. Right click and select "fill area" when finished.')

    imagedataorig=roifill;
    close
end