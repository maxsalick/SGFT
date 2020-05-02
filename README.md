# SGFT
Scanning gradient Fourier transform method for striated image quantification

The scanning gradient Fourier transform (SGFT) method incorporates gradient analysis along with fast Fourier transforms to determine regions of sarcomere organization within individual cells, and quantifies sarcomere organization and spacing at the sub-cellular scale.  Trends between sarcomere alignment, length, and location were observed in immature cardiomyocytes.  Utilization of the SGFT technique has been also demonstrated for additional applications, such as breast cancer collagen microstructure and neural rosette patterning.

Instructions for using Scanning Gradient Fourier Transform analysis code:

# Extract all files into the same directory
# Run sft_main_v026.m in MATLAB
This is the UI for running analysis on the striated images.  Use the "Add Image File(s)..." button to start loading grayscale TIFs for analysis.  You may set up analysis of one or several images at a time.

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
