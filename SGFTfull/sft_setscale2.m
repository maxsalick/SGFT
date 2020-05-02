function umperpix = sft_setscale(imagedata,vertres,imhor,imvert)

% Provides user with many options of determining scale of the image.
% v2 - utilizing gscan2 to make zooming and selecting more user-friendly
  umperpix=1;
 disp('Scale must be set...')
 disp('      1 - Use scalebar in image to calibrate')
 disp('      2 - Direct input of pixel/um ratio')
 disp('      3 - Set no scale (FOR TROUBLESHOOTING PURPOSES ONLY!)')
 disp('      4 - Use last-used pixel/um ratio')

scalemode = input('Select Option: ');
if scalemode ==1

    disp('Please select both ends of the scalebar.')
    disp('   LEFT CLICK = ZOOM IN')
    disp('   DOUBLE LEFT CLICK = ZOOM OUT')
    disp('   RIGHT CLICK = SELECT')
    disp('   LEFT CLICK AND HOLD = PAN')
    disp(' ')
    disp('Press any key to begin')
    pause
    
    figure('position',[50,50,floor((vertres-100)*imhor/imvert),vertres-100])
         colormap('gray')
    imagesc(imagedata);
    
    [scalepointx,scalepointy] = ginput2(2,'.r');
    
    
        close all
        scaleinput = input('Input the known distance between points: ');
        pixeldist = ((scalepointx(1)-scalepointx(2))^2+(scalepointy(1)-scalepointy(2))^2)^.5;
        umperpix = scaleinput/pixeldist;
        disp(['Selected scale is ' num2str(umperpix) ' microns per pixel'])
            save('umperpix','umperpix')
elseif scalemode==2
    umperpix = input('Input known microns-per-pixel: ');
    save('umperpix','umperpix')
elseif scalemode==3
    sprintf('Setting arbitrary calibration of 1 micron per pixel')
    umperpix = 1;
elseif scalemode==4
    load umperpix
    disp(['Selected scale is ' num2str(umperpix) ' microns per pixel'])
    
    
else
    sprintf('Whoops! Unknown input detected!');
    stop
end
     colormap('gray')
end
