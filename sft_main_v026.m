function [results] = sft_main_v026(inputpath,inputfiles, ...
    scale,patternsize,scanresolution)
% Scanning Gradient Fourier Transform -- Max Salick, July 14, 2018 
% 
% Function structure:
% function [results] = sft_main_v026(inputpath,inputfiles,scale,patternsize,scanresolution)
% 
% (1)inputpath = directory containing image file(s)
%             NOTE: for batch analysis, all files should be in the same directory
% 
% (2)inputfiles = file name to be analyzed
%             NOTE: currently only TIF files are supported; please ensure
%             that all files are converted to TIF
% 
% (3)scale = scale of image, as MICRONS PER PIXEL. Set to 1 to keep units
%             as "pixel"
% 
% (4)patternsize = approximate estimation of pattern size
%             NOTE: patterns with variable wavelength are still acceptable,
%             but an accurate initial estimate will improve ability to
%             identify patterned regions
% 
% (5)scanresolution = number of pixels to skip between each analysis scan
%             NOTE: small numbers will lead to higher resolution output, at
%             cost of increased computation time. This is also dependent on
%             image resolution (recommended to be around 16 for first round
%             of analysis)
% 
% 
% 
% This code combines gradient analysis (for direction determination) with
% one-directional Fourier transforms to quantitatively assess regions of repeating patterns within heterogeneous
% 2-D image sets. This code has been designed for determination of
% sarcomere organization within cardiomyocyte cultures, but may be utilized
% for a variety of applications within and outside of the field of biology.

% This code outputs 5 figures.  Figures 1 and 2 plot direction data 
% measuring counter-clockwise from the horizontal (positive x) axis.  
% Figures 3 and 4 plot direction data measuring counter-clockwise
% from the vertical (positive y) axis.  Figure 5 is a visual representation
% of these reference axes.  Only Figures 1-4 are saved for each run of the software.
% Examples of these plots are provided in the SGFT user manual.

% Identify software
disp(' ')
disp('---------------------------------------------')
disp('Sarcomere Analysis, Max R. Salick 7-14-2018')
disp('This software is free to use and edit by third parties, as long as')
disp('acknowledgement is given to the original author.')
disp('---------------------------------------------')
disp(' ')
%% File Input
% Detect whether individual or multiple files have been input
inputfilesize = size(inputfiles);
if inputfilesize(1) == 1
%     batchmode = 0; Keep for batch implementation
    pathfile = [inputpath inputfiles];
    if exist(pathfile) < 2
        pathfile = [inputpath inputfiles '.tif'];
        if exist(pathfile) < 2
        pathfile = [inputpath '/' inputfiles];
            if exist(pathfile) < 2
            pathfile = [inputpath '/' inputfiles '.tif'];
            end
        end
    end
    
    if exist(pathfile) == 0
        disp(['Unable to find file ' inputpath inputfiles]);
        exit;
    end
    
    
    disp([pathfile ' detected...']);
    
    
    
else
%     batchmode = 1; Keep for batch implementation
    
end

     
%%  Main Analysis Block      
    disp('SINGLE IMAGE mode selected...')
    disp(' ')
    
% Import and normalize image data
    [imagedata, imagesize, imhor, imvert] = sft_importim(pathfile);

% Set Blocksize
    blocksize = sft_setblocksize(patternsize,scale);

% Save Parameters
save('rerunparams','scanresolution','patternsize','blocksize','inputpath','inputfiles','imagedata','scale');

disp(' ');
disp('------------------------------------');
disp('PARAMETERS SET!  BEGINNING ANALYSIS.')
disp('------------------------------------');
disp(['File Name: ' inputfiles]);
disp(['Path Name: ' inputpath]);
disp(['um Per Pixel: ' num2str(scale)]);
disp(['Pattern Width Estimate: ' num2str(patternsize)]);
disp(['Scanjump: ' num2str(scanresolution)]);
disp(['Blocksize: ' num2str(blocksize)]);
disp('------------------------------------');
disp(' ');
 
% Local intensity filtering (Gaussian) and normalization
im = sft_imfilter(imagedata,blocksize);

% COVERAGE loop to identify regions that may contain patterns
[m_full_cov,im_noedge] = sft_cov_scan2(im,blocksize,scale,patternsize);

% Run gradient and pattern strength analysis on filtered image
[m_full_str, m_full_dir, quiver] = ...
    sft_scan_looper(imvert, imhor, m_full_cov, im_noedge, blocksize, ...
    scanresolution, scale);

% Organize and compile data into accessible arrays
[data, m_im, m_orig, m_cov, m_dir, m_quiver, m_str, fillbins, superiorang,percsarc] = ...
    sft_compile(m_full_cov, imagesize, blocksize, m_full_dir, ...
    m_full_str, quiver, im, imagedata, scanresolution, scale, pathfile, inputpath, inputfiles);

% Exporting data into worksheets
[percsarc, p20, p15, p10] = ...
    sft_export(inputpath, inputfiles, pathfile, fillbins, scale, ...
    blocksize, scanresolution, timer, data, superiorang, percsarc);



    disp('---------------------------------------------');
    disp(['Analysis completed for ' inputfiles '!']);
    disp('---------------------------------------------');
    disp(' ')
    disp('View data as...')
    disp('  - Excel worksheet (in same folder as original file)')
    disp('  - .mat file exported to same folder')
    disp('  - "Results" structure containing data, plots, and parameters')
    disp(' ')
    disp('Example: Enter command "imagesc(results.plots.str)"')

    
results.data = [data];
results.params.filename = inputfiles;
results.params.pathname = inputpath;
results.params.scale = scale;
results.params.scanresolution = scanresolution;
results.params.blocksize = blocksize;
% results.summary.OI = OI;
% results.summary.AI = AI;
% results.summary.CMI = CMI;
results.summary.superiorang = superiorang;
results.summary.sarcarea = percsarc;
results.summary.p20 = p20;
results.summary.p15 = p15;
results.summary.p10 = p10;
results.plots.im = m_im;
results.plots.str = m_str;
results.plots.dir = m_dir;
results.plots.quiv = m_quiver;
results.plots.align = fillbins;

%% plot reference axes

figure
pos1 = [0.1 0.35 0.4 0.4];
subplot('Position',pos1)
hold on
plot([1 0; 0 0], [0 0; 0 1],'k');
scatter(1,0,'k','>','filled')
scatter(0,1,'k','^','filled')

th = linspace( 0, pi/2, 100);
R = 0.3;
x = R*cos(th) + 0.5;
y = R*sin(th);
plot(x,y,'--k');
scatter(0.5,R,'k','<','filled')
scatter(0.75,0.3,'k','+')

annotation('textbox',[0.41 0.38 0.1 0.1],'String','0','EdgeColor','none');
annotation('textbox',[0.175 0.62 0.1 0.1],'String','90','EdgeColor','none');
annotation('textbox',[0.2 0.34 0.25 0.1],'String','from horizontal','EdgeColor','none');

axis([-0.5 1.5 -0.5 1.5]);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]);
xticks([])
yticks([])
plot([1.5 1.5], [-0.5 1.5],'k')
plot([-0.5 1.5], [1.5 1.5],'k')
hold off

pos2 = [0.55 0.35 0.4 0.4];
subplot('Position',pos2)
hold on
plot([-1 0; 0 0], [0 0; 0 1],'k');
scatter(-1,0,'k','<','filled')
scatter(0,1,'k','^','filled')

annotation('textbox',[0.59 0.38 0.1 0.1],'String','90','EdgeColor','none');
annotation('textbox',[0.83 0.62 0.1 0.1],'String','0','EdgeColor','none');
annotation('textbox',[0.67 0.34 0.25 0.1],'String','from vertical','EdgeColor','none');

th = linspace( pi/2, pi, 100);
R = 0.3;
x = R*cos(th);
y = R*sin(th) + 0.5;
plot(x,y,'--k');
scatter(-R,0.5,'k','v','filled')
scatter(-0.3,0.75,'k','+')

axis([-1.5 0.5 -0.5 1.5]);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]);
xticks([])
yticks([])
plot([0.5 0.5], [-1.5 1.5],'k')
plot([-1.5 0.5], [1.5 1.5],'k')
hold off

end

function [imagedata, imagesize, imhor, imvert] = sft_importim(pathfile)
% Imports image and its parameters, given a path and file name.
%     --------------------------------------------
% USER - for co-stain analysis, you may select out the R, G, or B channels.
%  See below for a quick example...

% imagedata = imagedata(:,:,1)   % this will only use the R channel for analysis
% imagedata = imagedata(:,:,2)   % this will only use the G channel for analysis
% imagedata = imagedata(:,:,3)   % this will only use the B channel for analysis

% Make sure to comment out line 2 if you do this.
%     --------------------------------------------
    imagedata = imread(pathfile);
    imagedata = max(max(imagedata(:,:,1),imagedata(:,:,2)),imagedata(:,:,3));
    

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

function blocksize = sft_setblocksize(patternsize,scale)
% Sets blocksize to include ~8 patterns.  This was found to be an optimal
% number of patterns within each scan to give an adequate, but localized,
% representation of patterning in the region.

disp(' ')
disp('Setting blocksize to include ~8 pattern widths within each frame')

blocksize = 2*(floor((patternsize*10/scale)/2));

% Make blocksize non-prime.
while isprime(blocksize)==1
    blocksize=blocksize+1;
end
disp(['Setting BLOCKSIZE = ' num2str(blocksize)]);

end

function im = sft_imfilter(imagedata,blocksize)
% Applies local contrast filter to image, and normalizes results.
disp(' ')
disp('Applying local contrast filter to image.')
disp('This involves producing a Gaussian blur of the image.')
disp('It may take a few minutes for large images.')

im_b = 255-imagedata;
filterbox = floor(blocksize/12);

close all
MD =imfilter(im_b,fspecial('disk',filterbox),'same');
im = im_b - MD;

% Normalize local contrasted image
disp('Normalizing filtered image.')
im = (im-min(min(im))).*(255/double(max(max(im-min(min(im))))));




end

function [m_full_cov,im_noedge] = sft_cov_scan2(im,blocksize,scale,patternsize)
% This function scans the image and determines regions that have enough
% positive signal to be included in the analysis.
% 
% The output, m_full_cov, is an array with 1s indicating regions with
% adequate signal, and 0s indicating regions that will be skipped to reduce
% computation time of Fourier analysis.
% 
% v2 - changed to include patestimate./umperpix for dilation to ensure that
% dilated positive region will span negative regions

disp('Beginning preliminary COVERAGE scan loop...')
disp(' ')
disp('First scan loop will determine regions that contain cells.')
disp('Edges of positive regions will be reduced (within 0.75 blocksizes).')
disp('This reduces the false positives that result from edge effects during')
disp('Fast Fourier Transform image analysis.')

dilated = imdilate(im,strel('disk',floor(patternsize./scale)));
meanintense = ceil(mean(mean(dilated(dilated>0))));
BW = roicolor(dilated,ceil(meanintense/2),255);
BW2 = imfill(BW,'holes');
roierode = imerode(BW2,strel('disk',floor(patternsize./scale)));
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

function [m_full_str, m_full_dir, quiver] = ...
    sft_scan_looper(imvert, imhor, m_full_cov, im_noedge,  ... 
                    blocksize, scanresolution, scale)
%   This file coordinates "sft_scan.m", giving the scan locations and
%   spacings, and constructing output arrays based on the sft_scan values.
% 
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

[GX,GY] = gradient(im_noedge);

orientim = atan2(-GY,GX);
gradstrbig = (GY.^2+GX.^2).^.5;

imagesize = [imvert imhor];
m_full_str = zeros(imagesize);
m_full_dir = zeros(imagesize);
quiver = zeros(imagesize);
totaltoscan = sum(sum(m_full_cov));

tic
scanned=0;

mstart = zeros(scanresolution);
mstart(floor(scanresolution)/2,1:end)=1;

for xloc = ceil(blocksize/2)+1:scanresolution:imagesize(2)-(ceil(blocksize/2)+1)

    timer = toc;
    timeremain = timer*totaltoscan/scanned - timer;
    minremain = floor(timeremain/60);
    secremain = floor(timeremain - 60*minremain);
    disp([num2str(ceil(100*scanned/totaltoscan)) '% complete. ' ...
        num2str(minremain) ':' num2str(secremain) ' remaining.'])

    for yloc = ceil(blocksize/2)+1:scanresolution:imagesize(1)-(ceil(blocksize/2)+1)
 
      if m_full_cov(yloc, xloc) == 1
         scanned = scanned+scanresolution^2;
         section = im_noedge(yloc-blocksize/2:yloc+blocksize/2-1,xloc-blocksize/2:xloc+blocksize/2-1);
         orientsection = orientim(yloc-blocksize/2:yloc+blocksize/2-1,xloc-blocksize/2:xloc+blocksize/2-1);
         gradstr = gradstrbig(yloc-blocksize/2:yloc+blocksize/2-1,xloc-blocksize/2:xloc+blocksize/2-1);

         section = section-min(min(section)).*255/max(max(1+section-min(min(section))));
        
 [fourstrengthval,meandirect] = sft_scan3(section,orientsection,gradstr,blocksize);


m_full_str(yloc-scanresolution/2:yloc+scanresolution/2,xloc-scanresolution/2:xloc+scanresolution/2)=fourstrengthval;
m_full_dir(yloc-scanresolution/2:yloc+scanresolution/2,xloc-scanresolution/2:xloc+scanresolution/2)=meandirect;

m = imrotate(mstart,double(meandirect*180/pi),'bilinear','crop');
quiver(yloc-scanresolution/2:-1+yloc+scanresolution/2,xloc-scanresolution/2:-1+xloc+scanresolution/2)=m;

       end
    end
   
end
                
      
disp(' ')
disp('Analysis complete!')        

       
               
end

function [fourstrengthval,meandirect] = sft_scan3(section,orientsection,gradstr,blocksize)
%  Scans through the entire image to detect patterning.  First conducts
%  gradient analysis to determine directionality of the image.  Next, 1-D
%  Fourier transforms are taken in the determined directions, and peak
%  detection of averaged Fourier transforms is conducted to assess the
%  amount of organization at the scanned areas.  Spacing between patterns
%  is determined spatially using intensity profiles across the image in the
%  known direction of patterning.



C = sum(sum(gradstr.*cos(orientsection.*2))); % Center of mass X (circular form)
S = sum(sum(gradstr.*sin(orientsection.*2)));
meandirect = atan2(S,C);
if meandirect<0
    meandirect = meandirect+2*pi;
end
meandirect = meandirect/2;

%subplot(2,2,1)                                               %plots
%imagesc(orientsection)                                       %found in old code
%subplot(2,2,2)                                               %
%imagesc(gradstr)                                             %
%subplot(2,2,3)                                               %
%hist(orientsection(gradstr~=0),100)                          %
%subplot(2,2,4)                                               %
%gscatter(orientsection(gradstr~=0),gradstr(gradstr~=0))      %   
%pause(.0005)                                                 %                       

npixels = floor(blocksize/1.1);
clear fourlina
fourlina = {zeros(1,npixels) zeros(1,npixels) zeros(1,npixels) zeros(1,npixels) zeros(1,npixels)};
warning('off','signal:findpeaks:noPeaks');
k=0;
fourlinsum = 0;
intensetot = 0;


for Doff=-8:1:8
    k=k+1;
    xi = [blocksize/2-blocksize*cos(meandirect)/2.2-Doff*sin(meandirect); ...
          blocksize/2+blocksize*cos(meandirect)/2.2-Doff*sin(meandirect)];
    yi = [blocksize/2+blocksize*sin(meandirect)/2.2-Doff*cos(meandirect); ...
          blocksize/2-blocksize*sin(meandirect)/2.2-Doff*cos(meandirect)];
      

       fourlina{k} = improfile(section,double(xi),double(yi),npixels);
       FOUR{k} = abs(fft(fourlina{k}));

       fourlinsum = fourlinsum+FOUR{k};
       intensetot = intensetot+sum(fourlina{k});
end

linesmade = k;

fourlin_ave = smooth(fourlinsum ./ linesmade);

fourlin_ave = fourlin_ave(1:ceil(0.5*npixels));
for k=length(fourlin_ave):-1:1
    fourlin_ave(k) = fourlin_ave(k)-min(fourlin_ave(1:k));
end

[pks,locs] = findpeaks(fourlin_ave,'sortstr','descend');

    if isempty(pks) == 0

        fourstrengthval = pks(1);
        
    else
        fourstrengthval = 0;

        
    end

str_four = fourstrengthval/intensetot;
if isnan(str_four) == 1
    str_four = 0;
end

fourstrengthval = str_four*90;




end

function [data, m_im, m_orig, m_cov, m_dir, m_quiver, m_str, fillbins, superiorang, percsarc] = ...
    sft_compile(m_full_cov, imagesize, blocksize, m_full_dir, ...
    m_full_str, quiver, im, imagedata, scanresolution, scale, pathfile,inputpath,inputfiles)
% Takes output matrices and compiles them into data arrays.  Produces
% figures to report the resutls.

disp(' ')
disp('Producing report of results...')


m_cov = m_full_cov(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2));
ind = find(m_cov);

m_dir  = m_full_dir(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2));
m_quiver = quiver(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2));
m_str = (m_full_str(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2)).^1);
loweststrength = min(m_str(ind));
m_str = m_str + loweststrength.*(1-m_cov);

% plot m_str to scale
% m_str_max = max(max(m_str));
% m_str_scaled = 255.*(m_str./m_str_max);
% figure
% imshow(m_str_scaled,'DisplayRange',[0 255]);

m_im = (im(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2)));
m_orig = (255-imagedata(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2)));


m_full_dist = scale.*(bwdist(1-m_full_cov));
m_dist = m_full_dist(ceil(blocksize/2):imagesize(1)-1-ceil(blocksize/2),ceil(blocksize/2):imagesize(2)-1-ceil(blocksize/2));

% Make sparse matrices

m_sp_dist = m_dist(1:scanresolution:end,1:scanresolution:end);
m_sp_str = m_str(1:scanresolution:end,1:scanresolution:end);
m_sp_dir = m_dir(1:scanresolution:end,1:scanresolution:end);
m_sp_cov = m_cov(1:scanresolution:end,1:scanresolution:end);
[coveredY,coveredX] = size(m_sp_dist);
[Ymesh,Xmesh] = meshgrid(1:coveredY,1:coveredX);
Ymesh = Ymesh'; 
Ymesh = coveredY-Ymesh;
Xmesh = Xmesh';


[leng,wid]=size(m_full_cov);
for widthrow = 1:leng
    checkrow=[0 m_full_cov(widthrow,:) 0];
    rowends=strfind(checkrow,[1 0])-1;
	rowstrts=strfind(checkrow,[0 1]);
    if length(rowends)==1 && length(rowstrts)==1
    pixwidth(widthrow)=max(rowends-rowstrts+1);
    else
        pixwidth(widthrow)=0;
    end
end
lanewidth = mean(pixwidth(pixwidth>0))*scale;

data = [Xmesh(m_sp_cov>0.5) ...  % 1
                Ymesh(m_sp_cov>0.5)...  % 2
                m_sp_dist(m_sp_cov>0.5) ...  % 3
                m_sp_str(m_sp_cov>0.5) ...  %4
                m_sp_dir(m_sp_cov>0.5) ... %5
                m_sp_cov(m_sp_cov>0.5) .* lanewidth];  %6


%% Alignment Distribution

[n,bins] = histc(data(:,5).*180./pi,[0:1:180]);

%plot direction historgram
figure (1)                         
histogram('BinEdges',[0:1:181], 'BinCounts',[n]);
axis([0 180 0 inf]);
pbaspect([2 1 1])
title('Direction of pattern');
xlabel('Degrees from horizontal');
ylabel('Frequency');
print([inputfiles ' DirectionHistogram1.png'],'-dpng');
csvwrite([inputfiles ' DirectionData1.csv'],n)

% make all angles between 0 and 90
% n_plot = n;
% for i = 1:90
%     n_plot(i) = n_plot(i) + n_plot(182-i);
% end
% 
% n_plot = n_plot(1:91);

% make all angles between -90 and 90
n_plot = n;
n_plot(91:181) = n(1:91);
n_plot(1:90) = n(91:180);

%plot direction historgram
figure (2)                         
histogram('BinEdges',[-90:1:91], 'BinCounts',[n_plot],'Orientation','horizontal');
axis([0 inf -90 90]);
pbaspect([1 2 1])
yticks([-90 -60 -30 0 30 60 90]);
title('Direction of pattern');
ylabel('Degrees from horizontal');
xlabel('Frequency');
print([inputfiles ' DirectionHistogram2.png'],'-dpng');
csvwrite([inputfiles ' DirectionData2.csv'],n_plot)

figure (3)
histogram('BinEdges',[0:1:181], 'BinCounts',[n],'Orientation','horizontal');
set(gca,'ydir','reverse','xdir','reverse')
axis([0 inf 0 180]);
pbaspect([1 2 1])
title('Direction of pattern');
ylabel('Degrees from vertical');
xlabel('Frequency');
print([inputfiles ' DirectionHistogram3.png'],'-dpng');
csvwrite([inputfiles ' DirectionData3.csv'],n)

figure (4)                         
histogram('BinEdges',[-90:1:91], 'BinCounts',[n_plot]);
set(gca,'xdir','reverse')
axis([-90 90 0 inf]);
pbaspect([2 1 1])
xticks([-90 -60 -30 0 30 60 90]);
title('Direction of pattern');
xlabel('Degrees from vertical');
ylabel('Frequency');
print([inputfiles ' DirectionHistogram4.png'],'-dpng');
csvwrite([inputfiles ' DirectionData4.csv'],n_plot)

% plot polar histograms
figure (5)
polarhistogram(data(:,5),45)
title('Direction of pattern from horizontal');
print([inputfiles ' DirectionHistogram5.png'],'-dpng');
csvwrite([inputfiles ' DirectionData5.csv'],data(:,5))

figure (6)
polarhistogram(data(:,5),45)
title('Direction of pattern from vertical');
ax = gca;
ax.ThetaZeroLocation = 'top';
print([inputfiles ' DirectionHistogram6.png'],'-dpng');
csvwrite([inputfiles ' DirectionData6.csv'],data(:,5))

%%

fillbins = zeros(size(n));

for k =1:max(size(data))
    fillbins(bins(k)) = fillbins(bins(k)) + (data(k,4)); 
end
fillbins(end)=[];

[angpks,anglocs] = findpeaks(fillbins);
[toss,highestangpk] = max(angpks);
allangX = sum(fillbins'.*cos(pi.*[2:2:360]./180));
allangY = sum(fillbins'.*sin(pi.*[2:2:360]./180));
if allangX>0
    superiorang = atan(allangY/allangX)/2;
else
    superiorang = (atan(allangY/allangX) + pi)/2;
end

if superiorang<0
    superiorang = superiorang+pi;
end

superiorang = superiorang*180/pi;

% AI = ((sum(cos(data(:,5))).^2+sum(sin(data(:,5))).^2).^.5)/length(data(:,5));
%
%
%
% OI = mean(m_str(m_cov));
% CMI = OI*AI;

% data(:,7) = m_sp_cov(m_sp_cov>0.5) .* OI;
% data(:,8) = m_sp_cov(m_sp_cov>0.5) .* AI;
% data(:,9) = m_sp_cov(m_sp_cov>0.5) .* CMI;
percsarc = (nnz(m_str(m_cov)>.5)) / (nnz(m_cov)) * 100;

disp(' ')
disp('Data compiled and reported in "data" array:')
disp('  Col 1: X coordinate')
disp('  Col 2: Y coordinate')
disp('  Col 3: Coverage (1 = analyzed, 0 = bypassed)')
disp('  Col 4: Pattern Strength')
disp('  Col 5: Pattern Direction')
disp('  Col 6: Sample Lanewidth')
% disp('  Col 7: Sample Overall Strength Index')
% disp('  Col 8: Sample Overall Alignment Index')
% disp('  Col 9: Sample Combined Organization Index')

end

function [percsarc, p20, p15, p10] = ...
    sft_export(inputpath, inputfiles, pathfile, fillbins, scale, ...
    blocksize, scanresolution, timer, data, superiorang, percsarc)
% Exports data arrays into .m files that contain X, Y, distance from edge,
% pattern strength, and pattern spacing.  Also exports
% the data into an Excel workbook.
disp(['Exporting data to ' pathfile '.m'])
disp(['Exporting data to ' pathfile ' Excel Workbook'])

params.scale=scale;
params.blocksize=blocksize;
params.scanresolution=scanresolution;
params.timer=timer;



dataexport = {'x', 'y', 'edgedist(um)', 'strength', 'direction',...
    'lanewidth'};
dataexport2 = [data(:,1) data(:,2) data(:,3) ...
    data(:,4) data(:,5) data(:,6)];

paramexport = {'path' inputpath; 'file' inputfiles; 'umperpix' params.scale; ...
    'blocksize' params.blocksize; 'scanjump' params.scanresolution; ...
    'time' params.timer};

warning off MATLAB:xlswrite:AddSheet
xlswrite([inputfiles '.xls'], dataexport, 'Data Export');
xlswrite([inputfiles '.xls'], dataexport2, 'Data Export', 'A2');
xlswrite([inputfiles '.xls'], paramexport, 'Parameters');


superiorang = ceil(superiorang);
fillbins3 = [fillbins;fillbins;fillbins];
sumexport = {'% Area Containing Sarcomeres(>.5)' percsarc; ...
    'Superior Angle' superiorang;...
    '% Weighted Alignment Within 20 Degrees' 100*sum(fillbins3(160+superiorang:200+superiorang))/sum(fillbins);... 
'% Weighted Alignment Within 15 Degrees' 100*sum(fillbins3(165+superiorang:195+superiorang))/sum(fillbins);...
'% Weighted Alignment Within 10 Degrees' 100*sum(fillbins3(170+superiorang:190+superiorang))/sum(fillbins)};
xlswrite([inputfiles '.xls'], sumexport, 'Summary');

save([inputfiles '.mat'], 'inputpath','inputfiles','params','data','sumexport');


disp(' ');
% disp(['Average Sarcomere Strength: ' num2str(OI) ]);
disp(['% Area Containing Sarcomeres(>.5): ' num2str(percsarc) '%']);
disp(['Superior Angle: ' num2str(superiorang) ' degrees']);
disp(['% Weighted Alignment Within 20 Degrees: ' num2str(100*sum(fillbins3(160+superiorang:200+superiorang))/sum(fillbins)) '%']); 
disp(['% Weighted Alignment Within 15 Degrees: ' num2str(100*sum(fillbins3(165+superiorang:195+superiorang))/sum(fillbins)) '%']);
disp(['% Weighted Alignment Within 10 Degrees: ' num2str(100*sum(fillbins3(170+superiorang:190+superiorang))/sum(fillbins)) '%']);
% disp(' ')
% disp('-----------------------------------------')
% disp(['            ORGANIZATION INDEX: ' num2str(OI) ]);
% disp(['ALIGNMENT INDEX (MRL) (0 to 1): ' num2str(AI) ]);
% disp(['      COMBINED MYOFIBRIL INDEX: ' num2str(CMI) ]);
% disp('-----------------------------------------')
% disp(' ')


percsarc = percsarc/100;
p20 = sum(fillbins3(160+superiorang:200+superiorang))/sum(fillbins);
p15 = sum(fillbins3(165+superiorang:195+superiorang))/sum(fillbins);
p10 = sum(fillbins3(170+superiorang:190+superiorang))/sum(fillbins);
end
