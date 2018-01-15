# SGFT
Scanning gradient Fourier transform method for striated image quantification

The scanning gradient Fourier transform (SGFT) method incorporates gradient analysis along with fast Fourier transforms to determine regions of sarcomere organization within individual cells, and quantifies sarcomere organization and spacing at the sub-cellular scale.  Trends between sarcomere alignment, length, and location were observed in immature cardiomyocytes.  Utilization of the SGFT technique has been also demonstrated for additional applications, such as breast cancer collagen microstructure and neural rosette patterning.

Instructions for using Scanning Gradient Fourier Transform analysis code:

# Extract all files into the same directory
# Run sft_guifig.m in MATLAB
This is the GUI for running analysis on the striated images.  Use the "Add Image File(s)..." button to start loading grayscale TIFs for analysis.  You may set up analysis of one or several images at a time.
# Set the Microns/Pixel for all images
Select an image to display a preview.  If the scale is known, input it into the Microns/Pixel box.  If the image contains a scalebar, you may use the "Find From Scalebar..." button to open a dialog that utilizes the scalebar to determine the image scale.  Under the Scalebar dialog, place point 1 and 2 on opposite ends of the scalebar.  It is recommended to zoom in on the scalebar to make these selections as accurate as possible.  Once these points are selected, select "Compute" and type the length of the scalebar, in microns, in the input box.  
If all images were taken using the same objective/microscope settings, it is likely that this ratio is the same and you may copy the determined scale to all other images using the "Copy to All" button.
# Set the Scan Spacing for all images
This will set the SCANJUMP parameter within the code, determining the number of pixels that are skipped between scans.  This ultimately determines the resolution of the analyzed data sets (low SCANJUMP = higher resolution analysis).  This is restricted by computation time.  Select a SCANJUMP that leads to timely analysis, without making the resulting data too coarse.  For a typical 1000 X 1000 or 2000 X 2000 image, a value of ~16 works well, for example.
If analysis is taking too long, try increasing this number.  If analyzed data appears coarse, try reducing this number.
# Set the Block (Scan) Size for all images
This will set the BLOCKSIZE parameter within the code, determining the size of each individual scan as the code runs gradient/Fourier transform analysis.  The ideal BLOCKSIZE depends on the wavelength of the pattern being detected.  To get an initial estimate of a good BLOCKSIZE, use the "Estimate from Pattern Size" button, then input the approximate wavelength of the expected pattern (for example, ~2.0 µm for sarcomeres).
# Select all images to analyze, and select "Analyze Selected"
Analysis will run on the selected images.  This may take between a few seconds to several hours, depending on number of images, image resolution, and SCANJUMP size. 
Analyzed files will be added to the same directory as the originals, with "out_" appended before the filename.  For example, image_17d.tif will produce out_image_17d.tif in the same directory.
# Run sft_guidata.m in MATLAB
This is the data analysis portion of the code.  It will allow you to import a single or multiple analyzed data sets for assessment.
# Import "out_..." files using the Add Files... button
Select the produced datasets.  A slight delay may occur to allow the code to parse the data sets.  Once datasets are imported, you may select files and plot their various properties using the plotting tools along the bottom of the interface.  The Image Plot button shows the filtered image.  The Organization Plot displays the determined pattern strength (highly patterned regions should have a higher intensity).  Pattern Length Plot will show determined pattern lengths of the image regions.  Direction Plot will show values between 0 and pi, depending on the direction that has been calculated using gradient analysis.
# Combine data using "Combine Data" button
Once the appropriate file(s) have been selected, their values can be compiled into a single, large data set.  Once the data has been compiled, it may be saved and a filename may be given for the combined data.  The m-file that is produced this way can then be opened independently in MATLAB, allowing further investigation of the resulting data sets.
The data may also be plotted using the plotting tools on the right of the interface.  Simply select the variables to be compared, and choose PLOT.  If there is an extremely large number of data points, plotting may be slowed down; if this is the case, select a smaller number of data points to be shown by changing "Show All" to "Show 25%" or lower.
