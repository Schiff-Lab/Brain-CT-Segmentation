# Brain-CT-Segmentation
Volumetric analysis for Brain, CSF, Calcium and Subdural spaces

Required Matlab packages: Image Processing Toolbox, Statistics and Machine Learning Toolbox, Parallel Computing Toolbox  

Within MATLAB_R2021a, navigate to the downloaded folder.  Open the PennstateCTsegmentation_calcium_App.mlapp file, and run the app to segment CT scans of the brain.  Within the App Designer, the code can be compiled into a standalone desktop application or a deployable web app.  

To use the standalone app provided here, download and install the application found in the MyAppInstaller file after downloading the matlab runtime compiler (found at https://www.mathworks.com/products/compiler/matlab-runtime.html - version R2021a for Mac).

Within earlier Matlab versions, open and run the PennstateCTsegmentation_calcium.m file.  

Once the app has been opened, press the Load CT Stack button.  A gui will open, and you can choose any dicom image out of the stack of a head CT scan set.  Once selected,the CT stack will be preprocessed and automatically segmented (into brain, CSF, and subdural is applicable).  Going through the slices, labels can be added or removed using the manual tools on the right side of the app.  The buttons toward the top right all use polygon selection (finished with a double click) to determine the area that is to be changed.  The buttons to the bottom right all use a single right-click in the center of a label region to change that region to the selected label.  

The threshold between brain and CSF can be adjusted using the Intensity Thresholding tool in the lower left.  Automated intra-cranial calcium segmentation can be accomlished using the Auto Segment Calcium button on the left (with more manual calcium tools on the right).  Be aware that shunts will be segmented as calcium when the automatic tool is used.

Once the scan set has been satisfactorily segmented, the Calculate Volumes button will determine the compartment volumes, and the Save Results button can be pressed. A .tiff file with the image segmentations and a .txt file with the volumes will be saved in the CT scan dicom directory.  A gui window will open, in which a name and directory can be chosen for the results .mat file.  This results .mat file can be opened using the Load Results button in order to make any changes to the segmentation.

Once the segmentation is done, the volume results from the .mat file can be used in this [Brain Volume Web Application](https://www.schiff-lab-webapps.esm.psu.edu/shinyappBV/) to determine age and sex-adjusted z-scores for the brain and CSF volumes, as well as to plot the volumes on normal growth curves.  If using this application for publishing, please cite this [paper from the Schiff lab](https://thejns.org/pediatrics/view/journals/j-neurosurg-pediatr/28/4/article-p458.xml).
