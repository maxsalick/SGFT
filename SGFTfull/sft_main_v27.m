% Scanning Fourier Transform
% Max Salick, 7-7-2014
% 

% This code combines gradient analysis (for direction determination) with
% one-directional Fourier transforms (for pattern detection) to
% quantitatively assess regions of repeating patterns within heterogeneous
% 2-D image sets. This code has been designed for determination of
% sarcomere organization within cardiomyocyte cultures, but may be utilized
% for a variety of applications within and outside of the field of biology.

%% STARTUP
% Fresh start
clc
close all
clear all

% Set initial values
data = [];

% Acquire resolution data for current monitor
res = get(0,'ScreenSize');
horres = res(3);
vertres = res(4);

% Identify software
disp('Sarcomere Analysis, Max R. Salick 7-7-2014')
disp('This software is free to use and edit by third parties, as long as')
disp('acknowledgement is given to the original author.')
disp('------------------------------')

disp(' ')
disp('Is this a [RERUN]? (copy all parameters saved from previous run)')
rerun = input('>> [YES = 1], [NO = 0] :');

if rerun == 0

disp(' ')
disp('Is this being run on a [SINGLE IMAGE] or a [BATCH OF IMAGES]?');

batchmode = input('>> [SINGLE = 0], [BATCH = 1] : ');

else
    batchmode = 0;
end

switch batchmode
%% SINGLE START
  case 0
      if rerun == 0
    disp('SINGLE IMAGE mode selected.')
    disp(' ')
    disp('>> Please locate image.  Image must be in TIF format.')

    [FileName,PathName] = uigetfile('*.tif','Multiselect','off','Locate Image');
    
    disp('Image located.')
    disp(' ')
    disp('Please provide estimate of [PATTERN SIZE] (distance in um between neighboring stripes).');
    disp('This helps to optimize the scanning window size, improving signal-to-noise.');
    patestimate = input('>> PATTERN SIZE (default value = 2) : ');
    disp(' ')
    disp('Please provide a [SCANNING RESOLUTION].');
    disp('This is the distance between pixels of neighboring scans.');
    disp('Smaller values provide higher-resolution results, at the cost of computation time.');
    scanjump = input('>> SCANNING RESOLUTION (default value = 16) : ');
    disp(' ')
    disp('Override direction detection?')
    disp('You may input a direction (in degrees) that will force the software to determine')
    disp('sarcomere strength in that direction for all points.  Good for patterned substrates.')
    overridebin = input('>> [DETECT DIRECTION = 0], [OVERRIDE DIRECTION = 1] : ');
    if overridebin == 1
    disp('Input direction (0-180 degrees) (0 = left-right, 90 = up-down)')
    overridedir = input('>> OVERRIDE DIRECTION (default value = 90) : ');
    else
        overridedir = 0;
    end
      else
        load rerunparams
      end

% Import and normalize image data

[imagedata, imagesize, imhor, imvert] = sft_importim(PathName, FileName);

%% Set Scale        

if rerun == 0
umperpix = sft_setscale2(imagedata,vertres,imhor,imvert);

%% Set Blocksize

blocksize = sft_setblocksize(patestimate,umperpix);


%% Remove Scalebar

imagedataorig = sft_remscalebar(imagedata,imhor,imvert,vertres);

end

%% Save rerun parameters
save('rerunparams','scanjump','patestimate','blocksize','FileName','PathName','imagedataorig','umperpix','overridebin','overridedir');

disp(' ');
disp('------------------------------------');
disp('PARAMETERS SET!  BEGINNING ANALYSIS.')
disp('------------------------------------');
disp(['File Name: ' FileName]);
disp(['Path Name: ' PathName]);
disp(['um Per Pixel: ' num2str(umperpix)]);
disp(['Pattern Width Estimate: ' num2str(patestimate)]);
disp(['Scanjump: ' num2str(scanjump)]);
disp(['Blocksize: ' num2str(blocksize)]);
disp(['Override Direction Active: ' num2str(overridebin)]);
disp(['Override Direction: ' num2str(overridedir)]);
disp('------------------------------------');
disp(' ');

close all


%% Local Contrast Filtering and Normalizing

im = sft_imfilter(imagedataorig,blocksize);


%% COVERAGE loop


[m_full_cov,im_noedge] = sft_cov_scan2(im,scanjump,blocksize,umperpix,patestimate);
im = im_noedge;

imagesc(im)

%% GRADIENT + FFT Loop


[imagesize, m_full_str, m_full_dir, m_full_sl, quiver] = ...
    sft_scan_looper(imvert, imhor, m_full_cov, im, blocksize, scanjump, umperpix, overridebin, overridedir);

           

%% Reporting and Data Summarization


[data, m_im, m_orig, m_cov, m_sl, m_dir, m_quiver, m_str,fillbins, superiorang, OI, AI, CMI, percsarc] = ...
    sft_compile(m_full_cov, imagesize, blocksize, m_full_sl, m_full_dir, ...
    m_full_str, quiver, im, imagedataorig, scanjump, umperpix);


%% Save Data

[ave_sl, std_sl, p20, p15, p10] = sft_export(PathName, FileName, m_str, m_cov, fillbins, umperpix, ...
    blocksize, scanjump, timer, data, superiorang, OI, AI, CMI, percsarc);

    
    disp('---------------------------------------------');
    disp(['Analysis completed for ' FileName '!']);
    disp('---------------------------------------------');

    
results.data = [data];
results.params.filename = FileName;
results.params.pathname = PathName;
results.params.umperpix = umperpix;
results.params.scanjump = scanjump;
results.params.blocksize = blocksize;
results.summary.OI = OI;
results.summary.AI = AI;
results.summary.CMI = CMI;
results.summary.superiorang = superiorang;
results.summary.sarcarea = percsarc;
results.summary.sl_ave = ave_sl;
results.summary.sl_std = std_sl;
results.summary.p20 = p20;
results.summary.p15 = p15;
results.summary.p10 = p10;
results.plots.im = m_im;
results.plots.str = m_str;
results.plots.dir = m_dir;
results.plots.quiv = m_quiver;
results.plots.sl = m_sl;
results.plots.align = fillbins;


%% ----------------------------------------------------------    
%% BATCH MODE    
%% ----------------------------------------------------------    
%% Load Batch File Names
    case 1
     disp('BATCH mode selected.')
     disp(' ')
     
     disp('Please select from the following options:')
     disp('NOTE: Images must be in same folder, and must have same parameters)')
     disp(' ')
     disp('1 - Continue currently loaded batch from last uncompleted image')
if exist('currentjob.m')
     disp('    (currently loaded batch detected)');
else
     disp('    (currently loaded batch not detected. This option is not recommended)')
end
     disp(' ')
     disp('2 - Start new batch by selecting images manually')
   
     disp(' ')
     disp('3 - Load in batch info from text file containing image names')
     disp(' ')
     batchtype = input('>> [CONTINUE = 1], [NEW = 2], [NEW FROM TXT = 3] : ');
     disp(' ')
     switch batchtype
         case 1
             
          if exist('currentjob.m')
    batch = importdata('currentjob.m');
    PathName = importdata('currentpath.m');
    disp('Reading existing job status information...')
    [batch.importfilenames' batch.importfilestatus'];
    totalimagestorun = length(batch.importfilenames);
          else
             disp('No existing job status information.  Starting new job batch.')
             batchtype=2;
             disp(' ')
             disp('Please locate images.  Image must be in TIF format.')   
             disp('Images must all be in the same folder.')
             
          [FileName,PathName] = uigetfile('*.tif','Multiselect','on','Locate Image');
          totalimagestorun = length(FileName);
          batch.importfilenames = FileName;
          for ind1 = 1:totalimagestorun
           batch.importfilestatus{ind1} = 0;
          end 
          
          save('currentjob.m','batch');
          save('currentpath.m','PathName');
          [batch.importfilenames' batch.importfilestatus'];
       
          end
         
          case 2
          disp('Please locate images.  Image must be in TIF format.')   
          disp('Images must all be in the same folder.')
             
          [FileName,PathName] = uigetfile('*.tif','Multiselect','on','Locate Image');
          totalimagestorun = length(FileName);
          batch.importfilenames = FileName;
          for ind1 = 1:totalimagestorun
           batch.importfilestatus{ind1} = 0;
          end 
          
          save('currentjob.m','batch');
          save('currentpath.m','PathName');
          [batch.importfilenames' batch.importfilestatus'];
             
         case 3

        disp('Please select the job file (text file containing file names of images).')
        disp('The job file should be in the same directory as the images.')
        disp('The images should be listed with extensions, and should be ENTER-delimited.')
        [jobfilename,PathName] = uigetfile('*.txt','Multiselect','off','Locate Image');

        fid = fopen([PathName jobfilename],'r');
        batch.importfilenamesa = textscan(fid,'%s');
        batch.importfilenames = batch.importfilenamesa{1};
        totalimagestorun = length(batch.importfilenames);
          for ind1 = 1:totalimagestorun
           batch.importfilestatus{ind1} = 0;
          end 
          
          save('currentjob.m','batch');
          save('currentpath.m','PathName');
          [batch.importfilenames' batch.importfilestatus'];
     end

    save('currentjob.m','batch');
    save('currentpath.m','PathName');
    
    disp(' ')
    disp('Current job status (0=PENDING 1=COMPLETE 2=ERROR):')
    disp(' ')
    [batch.importfilenames' batch.importfilestatus']


%% Set Batch Parameters
    
    if batchtype == 1
        batchparams = importdata('batchparams.m');
        umperpix = batchparams.umperpix;
        scanjump = batchparams.scanjump;
        blocksize = batchparams.blocksize;
        patestimate = batchparams.patestimate;
        overridebin = batchparams.overridebin;
        overridedir = batchparams.overridedir;
        
    else
        disp(' ')
        disp('Batch parameters must be set.')
        disp('Once set, these parameters will be applied to ALL images in the batch.')
        disp(' ')
    
        disp(' ')
    disp('Please provide estimate of [PATTERN SIZE] (distance in um between neighboring stripes).');
    disp('This helps to optimize the scanning window size, improving signal-to-noise.');
    patestimate = input('>> PATTERN SIZE (default value = 2) : ');
    disp(' ')
    disp('Please provide a [SCANNING RESOLUTION].');
    disp('This is the distance between pixels of neighboring scans.');
    disp('Smaller values provide higher-resolution results, at the cost of computation time.');
    scanjump = input('>> SCANNING RESOLUTION (default value = 16) : ');
    disp(' ')
    disp('Override direction detection?')
    disp('You may input a direction (in degrees) that will force the software to determine')
    disp('sarcomere strength in that direction for all points.  Good for patterned substrates.')
    overridebin = input('>> [DETECT DIRECTION = 0], [OVERRIDE DIRECTION = 1] : ');
    if overridebin == 1
    disp('Input direction (0-180 degrees) (0 = left-right, 90 = up-down)')
    overridedir = input('>> OVERRIDE DIRECTION (default value = 90) : ');
    else
        overridedir = 0;
    end
    FileName = [batch.importfilenames{1}];    
    
    [imagedata, imagesize, imhor, imvert] = sft_importim(PathName, FileName);

% Set Scale        

umperpix = sft_setscale2(imagedata,vertres,imhor,imvert)

% Set Blocksize

blocksize = sft_setblocksize(patestimate,umperpix)


% Save rerun parameters

    end
    
save('batchparams.m','scanjump','patestimate','blocksize','umperpix','overridebin','overridedir');

disp(' ');
disp('------------------------------------');
disp('PARAMETERS SET!  BEGINNING ANALYSIS.')
disp('------------------------------------');
disp(['Path Name: ' PathName]);
disp(['um Per Pixel: ' num2str(umperpix)]);
disp(['Pattern Width Estimate: ' num2str(patestimate)]);
disp(['Scanjump: ' num2str(scanjump)]);
disp(['Blocksize: ' num2str(blocksize)]);
disp(['Override Direction Active: ' num2str(overridebin)]);
disp(['Override Direction: ' num2str(overridedir)]);
disp('------------------------------------');
disp(' ');


totalimagestorun = length(batch.importfilenames);    

while 1==1
    
for ind1 = 1:totalimagestorun+1
    if ind1>totalimagestorun
        disp('Reached end of batch!  Please find completed data in the folder that contains the images.');
        return
    end
    checkjobstatus = batch.importfilestatus{ind1};
    if checkjobstatus == 0
        currentjob = ind1;
        FileName = [batch.importfilenames{currentjob}];
        break
    end
 
end

disp(' ')
disp('--------------------------------------------------')
disp(['Next file selected: ' FileName])
disp(['This is file ' num2str(currentjob) ' / ' num2str(totalimagestorun)]);
disp(datestr(now))
disp('--------------------------------------------------')
disp(' ')

try
  
close all

[imagedata, imagesize, imhor, imvert] = sft_importim(PathName, FileName);
imagedataorig = imagedata;

% Local Contrast Filtering and Normalizing
im = sft_imfilter(imagedataorig,blocksize);

% COVERAGE loop
[m_full_cov,im_noedge] = sft_cov_scan2(im,scanjump,blocksize,umperpix,patestimate);
im = im_noedge;

imagesc(im)
colormap(gray)
% GRADIENT + FFT Loop
[imagesize, m_full_str, m_full_dir, m_full_sl, quiver] = ...
    sft_scan_looper(imvert, imhor, m_full_cov, im, blocksize, scanjump, umperpix, overridebin, overridedir);

% Reporting and Data Summarization
[data, m_im, m_orig, m_cov, m_sl, m_dir, m_quiver, m_str,fillbins, superiorang, OI, AI, CMI, percsarc] = ...
    sft_compile(m_full_cov, imagesize, blocksize, m_full_sl, m_full_dir, ...
    m_full_str, quiver, im, imagedataorig, scanjump, umperpix);

% Save Data
sft_export(PathName, FileName, m_str, m_cov, fillbins, umperpix, ...
    blocksize, scanjump, timer, data, superiorang, OI, AI, CMI, percsarc)


    disp('---------------------------------------------');
    disp(['Analysis completed for ' FileName '!']);
    batch.importfilestatus{currentjob} = 1;
    disp('---------------------------------------------');
    save('currentjob.m','batch')


catch
    disp(['WARNING: ERROR HAS OCCURED FOR SAMPLE ' FileName]);
    batch.importfilestatus{currentjob} = 2;
     save('currentjob.m','batch')   



end

end

end



