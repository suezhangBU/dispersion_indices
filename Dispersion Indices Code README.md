**Dispersion Indices Code README**



CellProfiler, ImageJ and MATLAB are needed to run dispersion indices analyses. Save the following files to a desired location: ImageAnalysisPipeline\_081423\_v1.cpproj and statisticalDispersion\_DM\_SZ\_v2.m. 



**ImageJ**



Using ImageJ, create a maximum intensity projection of the image and split the channels. Channel 1 (C1) should be nuclear stain, channel 2 (C2) should be a stain for the cytoplasm, such as actin, and channel 3 (C3) is the subcellular structure to be quantified using dispersion indices. Save all the channels as separate .tif images.



**CellProfiler**



Use the .cpproj CellProfiler file to generate a cell mask for each image. Save the cell mask image as .tif and name the image with \[image name]\_CellMask. Ensure tha the cell mask image is in the same folder as the cell channel images. 





**MATLAB**



Open statisticalDispersion\_DM\_SZ\_v2.m in MATLAB. Change the path in Line 4 to the location where the .m file is saved. 



A prompt will ask for the bit depth of the image. Enter the bit depth of the image (usually 16 for widefield images and 12 for confocal images). A window pop-up will allow the user to select the directory where the .tif image files are saved. Select the C1 file of the image to be analyzed. The code will then analyze the dispersion indices of every cell in the image and will output a .mat file with an array of the values. 



The output will be an array as shown:



result = \[MLD\_mito, IMcov\_mito, Theil\_T\_mito, Gini\_mito,m\_mito, mActin, nucleiArea], where MLD\_mito is the mean log deviation, IMcov\_mito is the coefficient of variation, Theil\_T\_mito is the Theil's T value, Gini\_mito is the Gini coefficient, m\_mito is the mean of the pixel intensities of C3, mActin is the mean of the pixel intensities of C2, and nucleiArea is the area of all nuclei combined. 

