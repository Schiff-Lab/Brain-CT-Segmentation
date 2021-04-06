function varargout = PennstateCTsegmentation_calcium(varargin)
% clc;
% %clear variables;
% close all
% GUI_FLIS MATLAB code for GUI_FLIS.fig
%      GUI_FLIS, by itself, creates a new GUI_FLIS or raises the existing
%      singleton*.
%
%      H = GUI_FLIS returns the handle to a new GUI_FLIS or the handle to
%      the existing singleton*.
%
%      GUI_FLIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FLIS.M with the given input arguments.
%
%      GUI_FLIS('Property','Value',...) creates a new GUI_FLIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_FLIS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_FLIS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_FLIS

% Last Modified by GUIDE v2.5 27-Mar-2021 16:11:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PennstateCTsegmentation_calcium_OpeningFcn, ...
                   'gui_OutputFcn',  @PennstateCTsegmentation_calcium_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_FLIS is made visible.
function PennstateCTsegmentation_calcium_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FLIS (see VARARGIN)

% Choose default command line output for GUI_FLIS
clc
handles.output = hObject;
handles.flagStartSlice = 0;
handles.flagSubDuralSlice = 0;
handles.flagDisplaySeg = 0;
handles.flagSegFinish = 0;
handles.flagUndo = 0;
set(handles.Instruct,'Enable', 'Off')
set(handles.instructRegionGrow,'Enable', 'Off')
set(handles.regionInclude,'Enable','Off')
set(handles.regionBrain,'Enable','Off')
set(handles.regionCSF,'Enable','Off')
set(handles.regionSub,'Enable','Off')
set(handles.deleteRegion,'Enable','Off')
set(handles.addRegion,'Enable','Off')
set(handles.convertBrain,'Enable','Off')
set(handles.calculateVolumes,'Enable','Off')
set(handles.convertFluid,'Enable','Off')
set(handles.convertSub,'Enable','Off')
set(handles.saveResults,'Enable','Off')
set(handles.axes2,'visible','off') 
set(handles.axes3,'visible','off') 
set(handles.Undo,'Enable','Off')
set(handles.brainLeg,'visible', 'Off')
set(handles.csfLeg,'visible', 'Off')
set(handles.subLeg,'visible', 'Off')
set(handles.calcLeg,'visible', 'Off')
set(handles.calcium_seg,'Enable', 'Off')
set(handles.calcium_roi,'Enable', 'Off')
set(handles.threshold_button,'Enable', 'Off')
set(handles.threshold_slider,'Enable', 'Off')
set(handles.thresh_num,'visible', 'Off')
set(handles.shunt,'Enable', 'Off')
set(handles.axes5,'visible', 'Off')
set(handles.skull_calc_seg,'Enable', 'Off')
set(handles.calc_undo,'visible', 'Off')
set(handles.thresh_undo,'visible', 'Off')
%set(handles.volumeDisplay,'Enable','Off')
handles.code_dirname = cd;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_FLIS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PennstateCTsegmentation_calcium_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% clear variables
%clc
% close all
axes(handles.axes1)
cla
axes(handles.axes2)
cla
axes(handles.axes3)
cla
axes(handles.axes5)
cla
%% refreshing the code
clc
handles.output = hObject;
handles.flagStartSlice = 0;
handles.flagSubDuralSlice = 0;
handles.flagDisplaySeg = 0;
handles.flagSegFinish = 0;
%set(handles.FLIS,'Enable','On')
set(handles.Instruct,'Enable', 'Off')
set(handles.deleteRegion,'Enable','Off')
set(handles.addRegion,'Enable','Off')
set(handles.convertBrain,'Enable','Off')
set(handles.calculateVolumes,'Enable','Off')
set(handles.convertFluid,'Enable','Off')
set(handles.convertSub,'Enable','Off')
set(handles.saveResults,'Enable','Off')
set(handles.axes2,'visible','off') 
set(handles.axes3,'visible','off') 
set(handles.brainLeg,'visible', 'Off')
set(handles.csfLeg,'visible', 'Off')
set(handles.subLeg,'visible', 'Off')
set(handles.calcLeg,'visible', 'Off')
set(handles.instructRegionGrow,'Enable', 'Off')
set(handles.regionInclude,'Enable','Off')
set(handles.regionBrain,'Enable','Off')
set(handles.regionCSF,'Enable','Off')
set(handles.regionSub,'Enable','Off')
set(handles.calcium_seg,'Enable', 'Off')
set(handles.calcium_roi,'Enable', 'Off')
set(handles.skull_calc_seg,'Enable', 'Off')
set(handles.shunt,'Enable', 'Off')
set(handles.threshold_button,'Enable', 'Off')
set(handles.threshold_slider,'visible', 'Off')
set(handles.thresh_num,'visible', 'Off')
set(handles.axes5,'visible', 'Off')
set(handles.calc_undo,'visible', 'Off')
set(handles.thresh_undo,'visible', 'Off')
str1{1} = 'Select slices...';
set(handles.selectSlice,'Value',1);
%set(handles.startSlice,'Value',1);
%set(handles.segSliceNo,'Value',1);
%set(handles.subduralSlice,'Value',1);
%set(handles.displayCleanSlice,'Value', 0);
%set(handles.displaySeg,'Value', 0);
set(handles.volumeDisplay, 'String', 'Volumes not found'); 
%% start the code
[image_filename,image_dirname]=uigetfile(('*.ima;*.dcm'),'Select DICOM image from desired image stack');%selecting the directory
% Load image stack from selected image set:
handles.files = [dir(fullfile(image_dirname,'/*.ima'));dir(fullfile(image_dirname,'/*.dcm'))];
handles.info=dicominfo(fullfile(image_dirname, handles.files(1).name));
set(handles.File, 'String',['Patient ID:  ', handles.info.PatientID]); 
checkG = isfield(handles.info,'PatientSex');
if(checkG == 0)
    prompt = {'Enter Patient Gender (M or F)'};
    dlg_title = 'Gender';
    S = inputdlg(prompt,dlg_title);
    S = char(S);
else
   S = handles.info.PatientSex; 
end
checkA = isfield(handles.info,'PatientAge');
if(checkA == 0)
    prompt = {'Enter Patient Age in Months'};
    dlg_title = 'Age';
    S1 = inputdlg(prompt,dlg_title);
    S1 = [char(S1), 'M'];
else
   S1 = handles.info.PatientAge; 
end
strGender = ['Patient Gender: ',S];
strAge = ['Patient Age: ',S1];
ageVariable = S1(end);
ageLength = S1(1:end-1);
ageL = str2double(ageLength);
if(ageVariable == 'M')
    ageL = ageL*30;
elseif(ageVariable == 'W')
    ageL = ageL*7;
elseif(ageVariable == 'Y')
    ageL = ageL*365;
end
handles.Age = ageL;   
handles.Gender =  S;   
%set(handles.genderBox,'String',strGender);
%set(handles.ageBox,'String',strAge);
handles.imageDirName = image_dirname;
handles.rowNum = cast(handles.info.Rows,'double');
handles.colNum = cast(handles.info.Columns,'double');
handles.numSlices = length(handles.files);
handles.data = zeros(handles.rowNum,handles.colNum,handles.numSlices);
handles.dataOrig = zeros(handles.rowNum,handles.colNum,handles.numSlices);
handles.dataNormalized = zeros([handles.rowNum,handles.colNum,handles.numSlices]);
for p = 1:handles.numSlices
  handles.data(:,:,p) = dicomread(fullfile(image_dirname,handles.files(p).name));
  handles.dataOriginal(:,:,p) = handles.data(:,:,p);
  handles.dataOrig(:,:,p) = (handles.info.RescaleSlope.*handles.dataOriginal(:,:,p))+handles.info.RescaleIntercept;
  handles.data(:,:,p) = double((handles.data(:,:,p)*handles.info.RescaleSlope) + (handles.info.RescaleIntercept));%changing the scale of data
  Bin = (handles.data(:,:,p)>-50)&(handles.data(:,:,p)<=100);% excluding the non-brain content
  Bin1 = (handles.data(:,:,p)>100);
  Bin1 = Bin1*100;
  handles.data(:,:,p) = immultiply(double(handles.data(:,:,p)), double(Bin));
  handles.data(:,:,p) = double(handles.data(:,:,p)) + double(Bin1);
  handles.dataNormalized(:,:,p) = (handles.data(:,:,p) + 50)/150;
end
%axes(handles.axes1)
str1{1} = 'Select slices...';
for p = 1:handles.numSlices
    str1{p + 1} = ['Slice ',num2str(p)];
end
set(handles.selectSlice,'String',str1);
%set(handles.startSlice,'String',str1);
%set(handles.subduralSlice,'String',str1);
%set(handles.segSliceNo,'String',str1);
%slide show of images
% for p = 1 : handles.numSlices
%     handles.sliceNum = p;
%     imshow(handles.dataNormalized(:,:,handles.sliceNum),'Parent',handles.axes1);
%     str2 = ['Slice ',num2str(p)];
%     set(handles.MainBox,'String',str2);
%     set(handles.selectSlice,'Value',p+1);
%     pause(1);
% end
wb=waitbar(0,'Pre-Processing in progress...');
str2 = ['Slice ',num2str(1)];
set(handles.MainBox,'String',str2);
set(handles.selectSlice,'Value',2);
handles.threshk=0;
%% finding candidate region
handles.candidateRegion = zeros(512,512,handles.numSlices);
middleNum = round(handles.numSlices/2);
numberPartitions = ceil(handles.numSlices/6);
firstPartition = numberPartitions; % The first slices will be grouped together
firstTransition = firstPartition + floor((middleNum - firstPartition)/2); %slices of first transition
transMiddle = middleNum;
middleTrans = middleNum + (middleNum - firstTransition);
secondTransition = middleTrans + (firstTransition - firstPartition);
lastPartition = handles.numSlices;
%% alternative brain finding
X = handles.data(:);
meanX = mean(X(X > 0));
stdX = std(X(X > 0));
% for i = firstPartition + 1:numSlices
%     candidateRegion(:,:,i) = findBinaryNew(data(:,:,i), meanX);
%     figure(i);
%     imshow(candidateRegion(:,:,i), []);
% end
%% candidate region extraction
% first find it for the middle slice
wait = 0;
for slice = firstPartition + ceil((firstTransition - firstPartition)/2) + 1:middleTrans + ceil((secondTransition - middleTrans)/2)
    waitbar(wait/(handles.numSlices + 2),wb);
	wait = wait + 1;
    handles.candidateRegion(:,:,slice) = findBinaryMiddleSlice_1(handles.data(:,:,slice), meanX);
end
for slice = firstPartition + ceil((firstTransition - firstPartition)/2):-1:1
    waitbar(wait/(handles.numSlices),wb);
	wait = wait + 1;
    handles.candidateRegion(:,:,slice) = crossSliceBinary_Pre(handles.data(:,:,slice), handles.candidateRegion(:,:,slice + 1));
    %% commented for debugging
%     waitbar(wait/(handles.numSlices + 2),wb);
% 	wait = wait + 1;
%     handles.candidateRegion(:,:,slice) = crossSliceBinary(handles.data(:,:,slice), handles.candidateRegion(:,:,slice + 1));
end
for slice = middleTrans + ceil((secondTransition - middleTrans)/2) + 1:handles.numSlices
     waitbar(wait/(handles.numSlices),wb);
	wait = wait + 1;
    handles.candidateRegion(:,:,slice) = crossSliceBinaryNew_Pre(handles.data(:,:,slice), handles.candidateRegion(:,:,slice - 1), meanX);
    %% commented for debugging
%     waitbar(wait/(handles.numSlices + 2),wb);
% 	wait = wait + 1;
%     handles.candidateRegion(:,:,slice) = crossSliceBinaryNew(handles.data(:,:,slice), handles.candidateRegion(:,:,slice - 1), meanX);
end
%% removing Non Brain
for slice = firstPartition:-1:1
    waitbar(wait/(handles.numSlices + 2),wb);
	wait = wait + 1;
    handles.candidateRegion(:,:,slice) = removeNonBrain(handles.data(:,:,slice), handles.candidateRegion(:,:,slice));
end
%% Newly added code for automating beginning slice selection
checkBeg = ceil(firstPartition/2);
temp = handles.data;
maxVal = max(temp(:));
for slice = checkBeg-1:-1:1
    X1 = handles.data(:,:,slice);
    X2 = handles.data(:,:,slice + 1);
    len1 = length(find(X1 > 2/3*maxVal));
    len2 = length(find(X2 > 2/3*maxVal));
    checkVal = len1/len2;
    if(checkVal < .875)
        handles.segNo = slice+1;
        handles.flagStartSlice = 1;
        break;
    end
end
if(handles.flagStartSlice == 0)
   handles.segNo = 1;
   handles.flagStartSlice = 1; 
end
%% starting segmentation 
%% start segmenting the data
%start_spams
for slice = 1:handles.numSlices
    CR=handles.candidateRegion(:,:,slice);
    CR=imbinarize(CR);
    se = strel('cube',5);
    e = imerode(CR,se);
    ZZ=imfill(e,4,'holes'); 
    D = bwdist(ZZ);
    CR(D>0)=0;
    D=handles.candidateRegion(:,:,slice);
    D(CR==0)=0;
    handles.candidateRegion(:,:,slice)=D;
end
setup() ;
load('segNet_train_test');
netVerifySub = load('charscnn_large.mat') ;
handles.subNo = findSubduralNumber(handles.data, handles.candidateRegion, netVerifySub, net);
%load('meanVal');
%normalizing dictionaries
%% finding Kmean centers
temp = handles.data;
minInt = min(temp(:));
maxInt = max(temp(:));
%R = temp(handles.candidateRegion == 1);
thrsKm = min(temp(handles.candidateRegion == 1));
thrsKm  = max(thrsKm, 1/5*minInt);
%thrsMax = max(temp(handles.candidateRegion == 1));
thrsMax = 2/3*maxInt;%max(thrsMax, 2/3*maxInt);
L = (temp > thrsKm) & (temp < thrsMax) & (handles.candidateRegion == 1);
%L = (temp > thrsKm);
[~, centersKm] = kmeans(temp(L == 1), 2);
thrs2 = .5*centersKm(1) + .5*centersKm(2);
handles.thresh2=thrs2;
t = abs(centersKm(1) - centersKm(2));
[centersKm, ~] = sort(centersKm);
if(t < 12.5)
    centersKm(2) = centersKm(2) - 5;
    centersKm(1) = centersKm(1) - 12.5;
end
handles.centersKm = centersKm;
%% % parameters of sparse coefficents
param.mode = 0;
param.pos = true;
param.lambda = .15;
param.numThreads = -1;

patch = 10;
patchBy2 = patch/2;
lenTraining = 80;
handles.segmentedImage = zeros(512, 512, handles.numSlices);
%% extract the information from middle slice
tempImage = handles.data(:,:,middleNum);
candRegion = handles.candidateRegion(:,:,middleNum);
%addpath(genpath('SPAMS'))
[centers, subLabel] = findSubLabel(tempImage, candRegion, net);

handles.subLabel = subLabel;
%centers = [];
validData = handles.data(handles.candidateRegion == 1);
handles.thrs1 = mean(validData) + 3*std(validData);
handles.thrs = centersKm;
handles.sliceNum = 1;
close(wb);
imshow(handles.dataNormalized(:,:,handles.sliceNum),'Parent',handles.axes1);
%h = msgbox('To continue further, press "segment stack" !');

%% segmenting starts here 
set(handles.MainBox,'String','Segmentation Started!');
%% partitioning the data
%handles.candidateRegion = zeros(512,512,handles.numSlices);
wb=waitbar(0,'Segmentation in progress...');
middleNum = round(handles.numSlices/2);
numberPartitions = ceil(handles.numSlices/6);
firstPartition = numberPartitions; % The first slices will be grouped together
firstTransition = firstPartition + floor((middleNum - firstPartition)/2); %slices of first transition
transMiddle = middleNum;
middleTrans = middleNum + (middleNum - firstTransition);
secondTransition = middleTrans + (firstTransition - firstPartition);
lastPartition = handles.numSlices;
subLabelNum = handles.subNo;
segNo = handles.segNo;
handles.flagSegFinish = 1;
% %% alternative brain finding
% X = handles.data(:);
% meanX = mean(X(X > 0));
% stdX = std(X(X > 0));
% % for i = firstPartition + 1:numSlices
% %     candidateRegion(:,:,i) = findBinaryNew(data(:,:,i), meanX);
% %     figure(i);
% %     imshow(candidateRegion(:,:,i), []);
% % end
% %% candidate region extraction
% % first find it for the middle slice
% for slice = firstPartition + ceil((firstTransition - firstPartition)/2) + 1:middleTrans + ceil((secondTransition - middleTrans)/2)
%     handles.candidateRegion(:,:,slice) = findBinaryMiddleSlice_1(handles.data(:,:,slice), meanX);
% end
% for slice = firstPartition + ceil((firstTransition - firstPartition)/2):-1:1
%     handles.candidateRegion(:,:,slice) = crossSliceBinary(handles.data(:,:,slice), handles.candidateRegion(:,:,slice + 1));
% end
% for slice = middleTrans + ceil((secondTransition - middleTrans)/2) + 1:handles.numSlices
%     handles.candidateRegion(:,:,slice) = crossSliceBinaryNew(handles.data(:,:,slice), handles.candidateRegion(:,:,slice - 1), meanX);
% end
% %% removing Non Brain
% for slice = firstPartition:-1:1
%     handles.candidateRegion(:,:,slice) = removeNonBrain(handles.data(:,:,slice), handles.candidateRegion(:,:,slice));
% end
% 
% %% starting segmentation 
% %% start segmenting the data
% %start_spams
load('segNet_train_test');
%load('meanVal');
% %normalizing dictionaries
% %% finding Kmean centers
% temp = handles.data;
% minInt = min(temp(:));
% R = temp(handles.candidateRegion == 1);
% thrsKm = min(temp(handles.candidateRegion == 1));
% thrsKm  = max(thrsKm, 1/5*minInt);
% L = (temp > thrsKm);
% [~, centersKm] = kmeans(temp(handles.candidateRegion == 1), 2);
% [centersKm, ~] = sort(centersKm);
% handles.centersKm = centersKm;
% %% % parameters of sparse coefficents
% param.mode = 0;
% param.pos = true;
% param.lambda = .15;
% param.numThreads = -1;
% 
% patch = 10;
% patchBy2 = patch/2;
% lenTraining = 80;
% handles.segmentedImage = zeros(512, 512, handles.numSlices);
% %% extract the information from middle slice
% tempImage = handles.data(:,:,middleNum);
% candRegion = handles.candidateRegion(:,:,middleNum);
% dict = [dictBrain{6} dictFluid{6} dictSubDural{6}];
% w = [wBrain{6} wFluid{6} wSubDural{6}];
% meanDist = meanNormDist{6};
% meanInt = meanNormInt{6};
% %[centers, subLabel] = segmentationKsvdDFDL_findPar_final(tempImage, candRegion, dict,w,meanInt,meanDist);
% subLabelNum = handles.subNo;
% segNo = handles.segNo;
% subLabel = 2;
% centers = [];
% validData = handles.data(handles.candidateRegion == 1);
% thrs = mean(validData) + 3*std(validData);
% thrs = max(thrs, .667*max(validData));
%% segmentation for all images

%segments data partition wise
    for i = 1:ceil(firstPartition/2)
	    waitbar(i/(handles.numSlices),wb);
        if(i < segNo)
           handles.segmentedImage(:,:,i) = 0;
        else
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
       %handles.segmentedImage(:,:,i) = segFinal(tempImage, candRegion, handles.thrs1, handles.subLabel, net);
        handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        handles.segmentedImage(:,:,i) = segBegSlices(tempImage, candRegion, handles.thrs, handles.segmentedImage(:,:,i), handles.subLabel);
        end
        end
     end
    %creating second partition
    for i = ceil(firstPartition/2) + 1:firstPartition
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
            [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            %handles.segmentedImage(:,:,i) = segFinal(tempImage, candRegion, handles.thrs1, handles.subLabel, net);
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
            handles.segmentedImage(:,:,i) = segBegSlices(tempImage, candRegion, handles.thrs, handles.segmentedImage(:,:,i), handles.subLabel);
        end
        end
    end
    %creating third partition
    for i = firstPartition + 1:firstPartition + ceil((firstTransition - firstPartition)/2)
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
            [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating fourth partition
    for i = firstPartition + ceil((firstTransition - firstPartition)/2) + 1:firstTransition
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating fifth partition
    for i = firstTransition + 1:firstTransition + ceil((transMiddle - firstTransition)/2)
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating sixth partition
    for i = firstTransition + ceil((transMiddle - firstTransition)/2) + 1:transMiddle
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating seventh partition
    for i = transMiddle + 1:transMiddle + ceil((middleTrans - transMiddle)/2)
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating eigth partition
    for i = transMiddle + ceil((middleTrans - transMiddle)/2) + 1:middleTrans
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating nineth partition
    for i = middleTrans + 1:middleTrans + ceil((secondTransition - middleTrans)/2)
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating tenth partition
    for i = middleTrans + ceil((secondTransition - middleTrans)/2) + 1:secondTransition
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating eleventh partition
    for i = secondTransition + 1:secondTransition + ceil((lastPartition - secondTransition)/2)
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
       
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating twelveth partition
    for i = secondTransition + ceil((lastPartition - secondTransition)/2) + 1:lastPartition
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
%% display segmented region
close(wb);
cla(handles.axes1)
set(handles.axes1,'visible','off')
handles.sliceNum = 1;
set(handles.axes2,'visible','on') 
set(handles.axes3,'visible','on') 
%set(handles.displaySeg,'Value', 1);
axes(handles.axes2)
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.data(:,:,handles.sliceNum));
imshow(handles.dataNormalized(:,:,handles.sliceNum),'Parent',handles.axes3);
set(get(handles.axes3, 'xlabel'), 'string', 'Original Image')
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image') 
set(handles.MainBox,'String','Segmentation Complete!');
set(handles.deleteRegion,'Enable','On')
set(handles.Instruct,'Enable', 'On')
set(handles.addRegion,'Enable','On')
set(handles.convertBrain,'Enable','On')
set(handles.calculateVolumes,'Enable','On')
set(handles.convertFluid,'Enable','On')
set(handles.convertSub,'Enable','On')
set(handles.brainLeg,'visible', 'On')
set(handles.csfLeg,'visible', 'On')
set(handles.subLeg,'visible', 'On')
set(handles.calcLeg,'visible', 'On')
set(handles.regionBrain,'Enable','On')
set(handles.regionCSF,'Enable','On')
set(handles.regionSub,'Enable','On')
set(handles.regionInclude,'Enable','On')
set(handles.instructRegionGrow,'Enable', 'On')
set(handles.calcium_seg,'Enable', 'On')
set(handles.calcium_roi,'Enable', 'On')
set(handles.threshold_button,'Enable', 'On')
%set(handles.threshold_slider,'Enable', 'On')
set(handles.skull_calc_seg,'Enable', 'On')
%set(handles.shunt,'Enable', 'On')
set(handles.shunt,'Enable', 'On')
h = msgbox('Segmentation Complete! Now you can do manual corrections using the options towards the righ side of GUI. Once corrections are finished you can press calculate volumes');
%guidata(hObject, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object deletion, before destroying properties.
function axes1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in selectSlice.
function selectSlice_Callback(hObject, eventdata, handles)
% hObject    handle to selectSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selectSlice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectSlice
val1 = get(handles.selectSlice,'Value');% only brain is selected
if(handles.flagSegFinish == 0)
if(val1 == 1)
    set(handles.MainBox,'String','Select a Slice');
    guidata(hObject, handles);
elseif(get(handles.displaySeg,'Value')==0)
    handles.sliceNum = val1 - 1;
    str2 = ['Slice ',num2str(val1 - 1)];
    imshow(handles.dataNormalized(:,:,handles.sliceNum),'Parent',handles.axes1);
    set(handles.MainBox,'String',str2);
    guidata(hObject, handles);
else
    handles.sliceNum = val1 - 1;
    str2 = ['Slice ',num2str(val1 - 1)];
    dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.data(:,:,handles.sliceNum));
    set(handles.MainBox,'String',str2);
    guidata(hObject, handles);
end
else
    if(val1 == 1)
        cla(handles.axes1);
        set(handles.axes1,'visible','off')
        set(handles.MainBox,'String','Select a Slice');
        guidata(hObject, handles);
    else
        cla(handles.axes1);
        set(handles.axes1,'visible','off')
        handles.sliceNum = val1 - 1;
        str2 = ['Slice ',num2str(val1 - 1)];
        axes(handles.axes2)
        dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.data(:,:,handles.sliceNum));
        imshow(handles.dataNormalized(:,:,handles.sliceNum),'Parent',handles.axes3);
        set(get(handles.axes3, 'xlabel'), 'string', 'Original Image') 
        set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image') 
        set(handles.MainBox,'String',str2);
        guidata(hObject, handles);
    end
end
% --- Executes during object creation, after setting all properties.
function selectSlice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NextSlice.
function NextSlice_Callback(hObject, eventdata, handles)
% hObject    handle to NextSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = handles.sliceNum;
%val1 = get(handles.selectSlice,'Value');
str2 = ['Slice ',num2str(val + 1)];
if(handles.flagSegFinish == 0)
    if(val == handles.numSlices)
        set(handles.MainBox,'String','You have selected a slice that does not exist');
        guidata(hObject, handles);
    elseif(get(handles.displaySeg,'Value')==0)
        imshow(handles.dataNormalized(:,:,val + 1),'Parent',handles.axes1);
        % set(handles.selectSlice,'String', str2);
        set(handles.selectSlice, 'value', val+2);
        set(handles.MainBox,'String',str2);
        handles.sliceNum = val + 1;
        guidata(hObject, handles);
    else
        dataPlotAll(handles.segmentedImage(:,:,val + 1), handles.data(:,:,val + 1));
        % set(handles.selectSlice,'String', str2);
        set(handles.selectSlice, 'value', val+2);
        set(handles.MainBox,'String',str2);
        handles.sliceNum = val + 1;
        guidata(hObject, handles);
    end
else
    if(val == handles.numSlices)
        cla(handles.axes1);
        set(handles.axes1,'visible','off')
        set(handles.MainBox,'String','You have selected a slice that does not exist');
        guidata(hObject, handles);
    else
        cla(handles.axes1);
        set(handles.axes1,'visible','off')
        axes(handles.axes2)
        dataPlotAll(handles.segmentedImage(:,:,val + 1), handles.data(:,:,val + 1));
        imshow(handles.dataNormalized(:,:,val + 1),'Parent',handles.axes3);
        set(get(handles.axes3, 'xlabel'), 'string', 'Original Image')
        set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image')
        set(handles.selectSlice, 'value', val+2);
        set(handles.MainBox,'String',str2);
        handles.sliceNum = val + 1;
        guidata(hObject, handles);
    end
end

% --- Executes on button press in PreviousSlice.
function PreviousSlice_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = handles.sliceNum;
%val1 = get(handles.selectSlice,'Value');
str2 = ['Slice ',num2str(val - 1)];
if(handles.flagSegFinish == 0)
    if(val == 1)
        set(handles.MainBox,'String','You have selected a slice that does not exist');
        guidata(hObject, handles);
    elseif(get(handles.displaySeg,'Value')==0)
        imshow(handles.dataNormalized(:,:,val - 1),'Parent',handles.axes1);
        % set(handles.selectSlice,'String', str2);
        set(handles.selectSlice, 'value', val);
        set(handles.MainBox,'String',str2);
        handles.sliceNum = val - 1;
        guidata(hObject, handles);
    else
        dataPlotAll(handles.segmentedImage(:,:,val - 1), handles.data(:,:,val - 1));
        % set(handles.selectSlice,'String', str2);
        set(handles.selectSlice, 'value', val);
        set(handles.MainBox,'String',str2);
        handles.sliceNum = val - 1;
        guidata(hObject, handles);
    end
else
    if(val == 1)
        cla(handles.axes1);
        set(handles.axes1,'visible','off')
        set(handles.MainBox,'String','You have selected a slice that does not exist');
        guidata(hObject, handles);
    else
        guidata(hObject, handles);
        cla(handles.axes1);
        set(handles.axes1,'visible','off')
        axes(handles.axes2)
        dataPlotAll(handles.segmentedImage(:,:,val - 1), handles.data(:,:,val - 1));
        imshow(handles.dataNormalized(:,:,val - 1),'Parent',handles.axes3);
        set(get(handles.axes3, 'xlabel'), 'string', 'Original Image') 
        set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image')
        set(handles.selectSlice, 'value', val);
        set(handles.MainBox,'String',str2);
        handles.sliceNum = val - 1;
        guidata(hObject, handles);
    end
end


function MainBox_Callback(hObject, eventdata, handles)
% hObject    handle to MainBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MainBox as text
%        str2double(get(hObject,'String')) returns contents of MainBox as a double


% --- Executes during object creation, after setting all properties.
function MainBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MainBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FLIS.
function FLIS_Callback(hObject, eventdata, handles)
% hObject    handle to FLIS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.MainBox,'String','Segmentation Started!');
%% partitioning the data
%handles.candidateRegion = zeros(512,512,handles.numSlices);
wb=waitbar(0,'Segmentation in progress...');
middleNum = round(handles.numSlices/2);
numberPartitions = ceil(handles.numSlices/6);
firstPartition = numberPartitions; % The first slices will be grouped together
firstTransition = firstPartition + floor((middleNum - firstPartition)/2); %slices of first transition
transMiddle = middleNum;
middleTrans = middleNum + (middleNum - firstTransition);
secondTransition = middleTrans + (firstTransition - firstPartition);
lastPartition = handles.numSlices;
subLabelNum = handles.subNo;
segNo = handles.segNo;
handles.flagSegFinish = 1;
% %% alternative brain finding
% X = handles.data(:);
% meanX = mean(X(X > 0));
% stdX = std(X(X > 0));
% % for i = firstPartition + 1:numSlices
% %     candidateRegion(:,:,i) = findBinaryNew(data(:,:,i), meanX);
% %     figure(i);
% %     imshow(candidateRegion(:,:,i), []);
% % end
% %% candidate region extraction
% % first find it for the middle slice
% for slice = firstPartition + ceil((firstTransition - firstPartition)/2) + 1:middleTrans + ceil((secondTransition - middleTrans)/2)
%     handles.candidateRegion(:,:,slice) = findBinaryMiddleSlice_1(handles.data(:,:,slice), meanX);
% end
% for slice = firstPartition + ceil((firstTransition - firstPartition)/2):-1:1
%     handles.candidateRegion(:,:,slice) = crossSliceBinary(handles.data(:,:,slice), handles.candidateRegion(:,:,slice + 1));
% end
% for slice = middleTrans + ceil((secondTransition - middleTrans)/2) + 1:handles.numSlices
%     handles.candidateRegion(:,:,slice) = crossSliceBinaryNew(handles.data(:,:,slice), handles.candidateRegion(:,:,slice - 1), meanX);
% end
% %% removing Non Brain
% for slice = firstPartition:-1:1
%     handles.candidateRegion(:,:,slice) = removeNonBrain(handles.data(:,:,slice), handles.candidateRegion(:,:,slice));
% end
% 
% %% starting segmentation 
% %% start segmenting the data
% %start_spams
load('segNet_train_test');
%load('meanVal');
% %normalizing dictionaries
% %% finding Kmean centers
% temp = handles.data;
% minInt = min(temp(:));
% R = temp(handles.candidateRegion == 1);
% thrsKm = min(temp(handles.candidateRegion == 1));
% thrsKm  = max(thrsKm, 1/5*minInt);
% L = (temp > thrsKm);
% [~, centersKm] = kmeans(temp(handles.candidateRegion == 1), 2);
% [centersKm, ~] = sort(centersKm);
% handles.centersKm = centersKm;
% %% % parameters of sparse coefficents
% param.mode = 0;
% param.pos = true;
% param.lambda = .15;
% param.numThreads = -1;
% 
% patch = 10;
% patchBy2 = patch/2;
% lenTraining = 80;
% handles.segmentedImage = zeros(512, 512, handles.numSlices);
% %% extract the information from middle slice
% tempImage = handles.data(:,:,middleNum);
% candRegion = handles.candidateRegion(:,:,middleNum);
% dict = [dictBrain{6} dictFluid{6} dictSubDural{6}];
% w = [wBrain{6} wFluid{6} wSubDural{6}];
% meanDist = meanNormDist{6};
% meanInt = meanNormInt{6};
% %[centers, subLabel] = segmentationKsvdDFDL_findPar_final(tempImage, candRegion, dict,w,meanInt,meanDist);
% subLabelNum = handles.subNo;
% segNo = handles.segNo;
% subLabel = 2;
% centers = [];
% validData = handles.data(handles.candidateRegion == 1);
% thrs = mean(validData) + 3*std(validData);
% thrs = max(thrs, .667*max(validData));
%% segmentation for all images
%segments data partition wise
    for i = 1:ceil(firstPartition/2)
	    waitbar(i/(handles.numSlices),wb);
        if(i < segNo)
           handles.segmentedImage(:,:,i) = 0;
        else
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
       %handles.segmentedImage(:,:,i) = segFinal(tempImage, candRegion, handles.thrs1, handles.subLabel, net);
        handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        handles.segmentedImage(:,:,i) = segBegSlices(tempImage, candRegion, handles.thrs, handles.segmentedImage(:,:,i), handles.subLabel);
        end
        end
     end
    %creating second partition
    for i = ceil(firstPartition/2) + 1:firstPartition
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
            [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            %handles.segmentedImage(:,:,i) = segFinal(tempImage, candRegion, handles.thrs1, handles.subLabel, net);
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
            handles.segmentedImage(:,:,i) = segBegSlices(tempImage, candRegion, handles.thrs, handles.segmentedImage(:,:,i), handles.subLabel);
        end
        end
    end
    %creating third partition
    for i = firstPartition + 1:firstPartition + ceil((firstTransition - firstPartition)/2)
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
            [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating fourth partition
    for i = firstPartition + ceil((firstTransition - firstPartition)/2) + 1:firstTransition
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating fifth partition
    for i = firstTransition + 1:firstTransition + ceil((transMiddle - firstTransition)/2)
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating sixth partition
    for i = firstTransition + ceil((transMiddle - firstTransition)/2) + 1:transMiddle
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating seventh partition
    for i = transMiddle + 1:transMiddle + ceil((middleTrans - transMiddle)/2)
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating eigth partition
    for i = transMiddle + ceil((middleTrans - transMiddle)/2) + 1:middleTrans
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating nineth partition
    for i = middleTrans + 1:middleTrans + ceil((secondTransition - middleTrans)/2)
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating tenth partition
    for i = middleTrans + ceil((secondTransition - middleTrans)/2) + 1:secondTransition
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating eleventh partition
    for i = secondTransition + 1:secondTransition + ceil((lastPartition - secondTransition)/2)
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
       
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    %creating twelveth partition
    for i = secondTransition + ceil((lastPartition - secondTransition)/2) + 1:lastPartition
	    waitbar(i/(handles.numSlices),wb);
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, handles.centersKm);
        else
        
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            handles.segmentedImage(:,:,i) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
        end
    end
    
close(wb);
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
guidata(hObject, handles);
    
%% display segmented region
close(wb);
cla(handles.axes1)
set(handles.axes1,'visible','off')
handles.sliceNum = 1;
set(handles.axes2,'visible','on') 
set(handles.axes3,'visible','on') 
%set(handles.displaySeg,'Value', 1);
axes(handles.axes2)
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.data(:,:,handles.sliceNum));
imshow(handles.dataNormalized(:,:,handles.sliceNum),'Parent',handles.axes3);
set(get(handles.axes3, 'xlabel'), 'string', 'Original Image') 
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image')
set(handles.MainBox,'String','Segmentation Complete!');
set(handles.deleteRegion,'Enable','On')
set(handles.Instruct,'Enable', 'On')
set(handles.addRegion,'Enable','On')
set(handles.convertBrain,'Enable','On')
set(handles.calculateVolumes,'Enable','On')
set(handles.convertFluid,'Enable','On')
set(handles.convertSub,'Enable','On')
set(handles.saveResults,'Enable','On')
set(handles.calcium_seg,'Enable', 'On')
set(handles.calcium_roi,'Enable', 'On')
set(handles.threshold_button,'Enable', 'On')
set(handles.threshold_slider,'Enable', 'On')
set(handles.thresh_num,'visible', 'Off')
set(handles.axes5,'visible', 'Off')
set(handles.skull_calc_seg,'Enable', 'On')
set(handles.shunt,'Enable', 'On')
h = msgbox('Segmentation Complete! Now you can do manual corrections using the options towards the righ side of GUI. Once corrections are finished you can press calculate volumes');
guidata(hObject, handles);
% str2 = ['Slice ',num2str(1)];
% set(handles.MainBox,'String',str2);
% set(handles.selectSlice,'Value',2);
% guidata(hObject, handles);
% val1 = get(handles.selectSlice,'Value');% only brain is selected
% set(handles.displaySeg,'Value', 1);
% if(val1 == 1)
%     set(handles.MainBox,'String','Select a Slice');
%     guidata(hObject, handles);
% else
%     handles.sliceNum = val1 - 1;
%     str2 = ['Slice ',num2str(val1 - 1)];
%     dataPlotAll_axes(handles.segmentedImage(:,:,handles.sliceNum), handles.data(:,:,handles.sliceNum),handles.axes1);
%     set(handles.MainBox,'String',str2);
%     guidata(hObject, handles);
% end

% --- Executes on selection change in startSlice.
function startSlice_Callback(hObject, eventdata, handles)
% hObject    handle to startSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns startSlice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from startSlice
val1 = get(handles.startSlice,'Value');% only brain is selected
if(val1 == 1)
    set(handles.MainBox,'String','Select a Slice');
    guidata(hObject, handles);
else
    handles.segNo = val1 - 1;
    handles.flagStartSlice = 1;
    if(handles.flagStartSlice == 1 && handles.flagSubDuralSlice == 1)
        set(handles.FLIS,'Enable','On')
    end
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function startSlice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in subduralSlice.
function subduralSlice_Callback(hObject, eventdata, handles)
% hObject    handle to subduralSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns subduralSlice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from subduralSlice
val1 = get(handles.subduralSlice,'Value');% only brain is selected
if(val1 == 1)
    set(handles.MainBox,'String','Select a Slice');
    guidata(hObject, handles);
else
    handles.subNo = val1 - 1;
    handles.flagSubDuralSlice = 1;
    if(handles.flagStartSlice == 1 && handles.flagSubDuralSlice == 1)
        set(handles.FLIS,'Enable','On')
    end
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function subduralSlice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subduralSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in displayCleanSlice.
function displayCleanSlice_Callback(hObject, eventdata, handles)
% hObject    handle to displayCleanSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayCleanSlice
set(handles.displaySeg,'Value', 0);
imshow(handles.dataNormalized(:,:,handles.sliceNum),'Parent',handles.axes1);
guidata(hObject, handles);

% --- Executes on button press in displaySeg.
function displaySeg_Callback(hObject, eventdata, handles)
% hObject    handle to displaySeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displaySeg
set(handles.displayCleanSlice,'Value', 0);
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.data(:,:,handles.sliceNum));
guidata(hObject, handles);


% --- Executes on button press in deleteRegion.
function deleteRegion_Callback(hObject, eventdata, handles)
% hObject    handle to deleteRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.BW = roipoly();
X = handles.segmentedImage(:,:,handles.sliceNum);
set(handles.Undo,'Enable','On')
handles.undoSliceNum = handles.sliceNum;
handles.undoSegmented = X;
X(handles.BW) = 0;
handles.segmentedImage(:,:,handles.sliceNum) = X;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);


% --- Executes on button press in addRegion.
function addRegion_Callback(hObject, eventdata, handles)
% hObject    handle to addRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.BW = roipoly();
X = handles.segmentedImage(:,:,handles.sliceNum);
set(handles.Undo,'Enable','On')
handles.undoSliceNum = handles.sliceNum;
handles.undoSegmented = X;
temp = handles.data(:,:, handles.sliceNum);
maxVal = max(temp(:));
check = temp < .9*maxVal;
binVal = handles.BW & check;
handles.candidateRegion(:,:,handles.sliceNum)  = handles.candidateRegion(:,:,handles.sliceNum) | binVal;
[xCoord, yCoord] = find(handles.BW & check );
if(isempty(xCoord))
   X(handles.BW) = 0;
else
    data = handles.data(:,:,handles.sliceNum);
    for k = 1:length(xCoord)
        intVal = data(xCoord(k), yCoord(k));
        distC1 = abs(intVal - handles.centersKm(1));
        distC2 = abs(intVal - handles.centersKm(2));
        [~, lab] = min([distC2, distC1]);
        X(xCoord(k), yCoord(k)) = lab;
    end
end
handles.segmentedImage(:,:,handles.sliceNum) = X;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);

% --- Executes on button press in convertFluid.
function convertFluid_Callback(hObject, eventdata, handles)
% hObject    handle to convertFluid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.BW = roipoly();
X = handles.segmentedImage(:,:,handles.sliceNum);
C = handles.candidateRegion(:,:,handles.sliceNum);
set(handles.Undo,'Enable','On')
handles.undoSliceNum = handles.sliceNum;
handles.undoSegmented = X;
X(handles.BW & C == 1) = 2;
handles.segmentedImage(:,:,handles.sliceNum) = X;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);

% --- Executes on button press in convertBrain.
function convertBrain_Callback(hObject, eventdata, handles)
% hObject    handle to convertBrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.BW = roipoly();
X = handles.segmentedImage(:,:,handles.sliceNum);
C = handles.candidateRegion(:,:,handles.sliceNum);
set(handles.Undo,'Enable','On')
handles.undoSliceNum = handles.sliceNum;
handles.undoSegmented = X;
X(handles.BW & C == 1) = 1;
handles.segmentedImage(:,:,handles.sliceNum) = X;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);


% --- Executes on button press in calculateVolumes.
function calculateVolumes_Callback(hObject, eventdata, handles)
% hObject    handle to calculateVolumes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.saveResults,'Enable','On')
X = handles.segmentedImage;
brainTissue = find(X == 1);
fluidCSF = find(X == 2);
subdural = find(X == 3);
calcification = find(X == 4);
voxelVol = (handles.info.PixelSpacing(1)*handles.info.PixelSpacing(2)*handles.info.SliceThickness);
handles.brainVol = length(brainTissue)*voxelVol/10^3;
handles.fluidVol = length(fluidCSF)*voxelVol/10^3;
handles.subVol =  length(subdural)*voxelVol/10^3;
handles.calc =  length(calcification)*voxelVol/10^3;
set(handles.volumeDisplay, 'Max',2);
set(handles.volumeDisplay,'String',sprintf(['Tissue volume = ',num2str(handles.brainVol),' cm^3   ',  '   Fluid volume = ',num2str(handles.fluidVol),' cm^3  ',  ' \n Subdural volume = ',num2str(handles.subVol),' cm^3   ','Calcification volume = ',num2str(handles.calc),' cm^3']));
h = msgbox('Volumes are calculated! Now you can save results and load a new CT stack');
guidata(hObject, handles);

 %%%%REMOVE FOR FINAL%%%%%
% X = handles.segmentedImage35;
% brainTissue35 = find(X == 1);
% fluidCSF35 = find(X == 2);
% subdural35 = find(X == 3);
% calcification35 = find(X == 4);
% voxelVol = (handles.info.PixelSpacing(1)*handles.info.PixelSpacing(2)*handles.info.SliceThickness);
% handles.brainVol35 = length(brainTissue35)*voxelVol/10^3;
% handles.fluidVol35 = length(fluidCSF35)*voxelVol/10^3;
% handles.subVol35 =  length(subdural35)*voxelVol/10^3;
% handles.calc35 =  length(calcification35)*voxelVol/10^3;
% guidata(hObject, handles);
% 
% X = handles.segmentedImage45;
% brainTissue45 = find(X == 1);
% fluidCSF45 = find(X == 2);
% subdural45 = find(X == 3);
% calcification45 = find(X == 4);
% voxelVol = (handles.info.PixelSpacing(1)*handles.info.PixelSpacing(2)*handles.info.SliceThickness);
% handles.brainVol45 = length(brainTissue45)*voxelVol/10^3;
% handles.fluidVol45 = length(fluidCSF45)*voxelVol/10^3;
% handles.subVol45 =  length(subdural45)*voxelVol/10^3;
% handles.calc45 =  length(calcification45)*voxelVol/10^3;
% guidata(hObject, handles);
% 
% X = handles.segmentedImage50;
% brainTissue50 = find(X == 1);
% fluidCSF50 = find(X == 2);
% subdural50 = find(X == 3);
% calcification50 = find(X == 4);
% voxelVol = (handles.info.PixelSpacing(1)*handles.info.PixelSpacing(2)*handles.info.SliceThickness);
% handles.brainVol50 = length(brainTissue50)*voxelVol/10^3;
% handles.fluidVol50 = length(fluidCSF50)*voxelVol/10^3;
% handles.subVol50 =  length(subdural50)*voxelVol/10^3;
% handles.calc50 =  length(calcification50)*voxelVol/10^3;
% guidata(hObject, handles);
 %%%%REMOVE FOR FINAL%%%%%

% --- Executes on button press in saveResults.
function saveResults_Callback(hObject, eventdata, handles)
% hObject    handle to saveResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
segmentedImage = handles.segmentedImage;
dataOriginal=handles.dataOrig;
data = handles.data;
fluidVol = handles.fluidVol;
brainVol = handles.brainVol;
subVol = handles.subVol;
calcVol=handles.calc;
threshold=handles.threshk;
sliceNo = handles.numSlices;
info = handles.info;
str_dir = string(handles.imageDirName);
str_dir = [handles.imageDirName, 'result_image.tif'];
str_vol = [handles.imageDirName, 'result_vol.txt'];
%str_dir = [handles.imageDirName, ['handles.info.PatientID,'result_image.tif'']];
%str_vol = [handles.imageDirName, [handles.info.PatientID,'result_vol.txt']];
str_growth = [handles.imageDirName, 'growthPlot1.jpg'];
fileID = fopen(str_vol,'w');
fprintf(fileID,'Brain Volume : %3.2f cm^3\r\n',brainVol);
fprintf(fileID,'CSF Volume : %3.2f cm^3\r\n',fluidVol);
fprintf(fileID,'Subdural Volume : %3.2f cm^3\r\n',subVol);
fprintf(fileID,'Calcium Volume : %3.2f cm^3\r\n',calcVol);
fprintf(fileID,'Threshold : %3.4f \r\n',threshold);
fclose(fileID);
maxD = max(handles.data(:));
minD = min(handles.data(:));
dataN = (data - minD)/(maxD - minD);
wb=waitbar(0,'Saving results...');
wait = 0;
for i = 1:sliceNo
    storeFiles = zeros(512, 1050, 3);
    waitbar(wait/(sliceNo + 2),wb);
    wait = wait + 1;
    D = zeros(512,512,3);
    D1 = zeros(512,512,3);
    seg = segmentedImage(:,:,i);
    tempData = dataN(:,:,i);
    BinTissue = seg == 1;
    BinFluid = seg == 2;
    BinSubDural = seg == 3;
    BinCalc = seg == 4;
    rgbImage1 = cat(3, tempData, tempData, tempData);
    storeFiles(1:512, 1:512,1:3) = rgbImage1;
    D(:,:,2)=BinTissue;
    D(:,:,1)=BinFluid;
    D(:,:,3)=BinSubDural;
    %D(:,:,1) = BinCalc;
    D1 = .25*D + rgbImage1;
    D22=D1(:,:,2);
    D11=D1(:,:,1);
    D33=D1(:,:,3);
    D22(BinCalc==1)=1;
    D11(BinCalc==1)=1;
    D33(BinCalc==1)=0;
    %D1(:,:,3)=.25*D(:,:,3) + rgbImage1(:,:,3);
    D1(:,:,1)=D11;
    D1(:,:,2)=D22;
    D1(:,:,3)=D33;
    storeFiles(1:512, 539:1050,1:3) = D1;
    imwrite(storeFiles, str_dir, 'WriteMode', 'append');
end
     %%%%REMOVE FOR FINAL%%%%%
% segmentedImage = handles.segmentedImage35;
% data = handles.data;
% fluidVol = handles.fluidVol35;
% brainVol = handles.brainVol35;
% subVol = handles.subVol35;
% calcVol=handles.calc35;
% threshold=handles.threshk;
% sliceNo = handles.numSlices;
% info = handles.info;
% str_dir = string(handles.imageDirName);
% str_dir = [handles.imageDirName, 'result_image35.tif'];
% str_vol = [handles.imageDirName, 'result_vol35.txt'];
% fileID = fopen(str_vol,'w');
% fprintf(fileID,'Brain Volume : %3.2f cm^3\r\n',brainVol);
% fprintf(fileID,'CSF Volume : %3.2f cm^3\r\n',fluidVol);
% fprintf(fileID,'Subdural Volume : %3.2f cm^3\r\n',subVol);
% fprintf(fileID,'Calcium Volume : %3.2f cm^3\r\n',calcVol);
% fprintf(fileID,'Threshold : %3.4f \r\n',threshold);
% fclose(fileID);
% maxD = max(handles.data(:));
% minD = min(handles.data(:));
% dataN = (data - minD)/(maxD - minD);
% for i = 1:sliceNo
%     storeFiles = zeros(512, 1050, 3);
%     D = zeros(512,512,3);
%     D1 = zeros(512,512,3);
%     seg = segmentedImage(:,:,i);
%     tempData = dataN(:,:,i);
%     BinTissue = seg == 1;
%     BinFluid = seg == 2;
%     BinSubDural = seg == 3;
%     BinCalc = seg == 4;
%     rgbImage1 = cat(3, tempData, tempData, tempData);
%     storeFiles(1:512, 1:512,1:3) = rgbImage1;
%     D(:,:,2)=BinTissue;
%     D(:,:,1)=BinFluid;
%     D(:,:,3)=BinSubDural;
%     %D(:,:,1) = BinCalc;
%     D1 = .25*D + rgbImage1;
%     D22=D1(:,:,2);
%     D11=D1(:,:,1);
%     D33=D1(:,:,3);
%     D22(BinCalc==1)=1;
%     D11(BinCalc==1)=1;
%     D33(BinCalc==1)=0;
%     %D1(:,:,3)=.25*D(:,:,3) + rgbImage1(:,:,3);
%     D1(:,:,1)=D11;
%     D1(:,:,2)=D22;
%     D1(:,:,3)=D33;
%     storeFiles(1:512, 539:1050,1:3) = D1;
%     imwrite(storeFiles, str_dir, 'WriteMode', 'append');
% end
%      
% segmentedImage = handles.segmentedImage45;
% data = handles.data;
% fluidVol45 = handles.fluidVol45;
% brainVol45 = handles.brainVol45;
% subVol45 = handles.subVol45;
% calcVol45=handles.calc45;
% threshold=handles.threshk;
% sliceNo = handles.numSlices;
% info = handles.info;
% str_dir = string(handles.imageDirName);
% str_dir = [handles.imageDirName, 'result_image45.tif'];
% str_vol = [handles.imageDirName, 'result_vol45.txt'];
% str_growth = [handles.imageDirName, 'growthPlot1.jpg'];
% fileID = fopen(str_vol,'w');
% fprintf(fileID,'Brain Volume : %3.2f cm^3\r\n',brainVol45);
% fprintf(fileID,'CSF Volume : %3.2f cm^3\r\n',fluidVol45);
% fprintf(fileID,'Subdural Volume : %3.2f cm^3\r\n',subVol45);
% fprintf(fileID,'Calcium Volume : %3.2f cm^3\r\n',calcVol45);
% fprintf(fileID,'Threshold : %3.4f \r\n',threshold);
% fclose(fileID);
% maxD = max(handles.data(:));
% minD = min(handles.data(:));
% dataN = (data - minD)/(maxD - minD);
% for i = 1:sliceNo
%     storeFiles = zeros(512, 1050, 3);
%     D = zeros(512,512,3);
%     D1 = zeros(512,512,3);
%     seg = segmentedImage(:,:,i);
%     tempData = dataN(:,:,i);
%     BinTissue = seg == 1;
%     BinFluid = seg == 2;
%     BinSubDural = seg == 3;
%     BinCalc = seg == 4;
%     rgbImage1 = cat(3, tempData, tempData, tempData);
%     storeFiles(1:512, 1:512,1:3) = rgbImage1;
%     D(:,:,2)=BinTissue;
%     D(:,:,1)=BinFluid;
%     D(:,:,3)=BinSubDural;
%     %D(:,:,1) = BinCalc;
%     D1 = .25*D + rgbImage1;
%     D22=D1(:,:,2);
%     D11=D1(:,:,1);
%     D33=D1(:,:,3);
%     D22(BinCalc==1)=1;
%     D11(BinCalc==1)=1;
%     D33(BinCalc==1)=0;
%     %D1(:,:,3)=.25*D(:,:,3) + rgbImage1(:,:,3);
%     D1(:,:,1)=D11;
%     D1(:,:,2)=D22;
%     D1(:,:,3)=D33;
%     storeFiles(1:512, 539:1050,1:3) = D1;
%     imwrite(storeFiles, str_dir, 'WriteMode', 'append');
% end
% 
% segmentedImage = handles.segmentedImage50;
% data = handles.data;
% fluidVol50 = handles.fluidVol50;
% brainVol50 = handles.brainVol50;
% subVol50 = handles.subVol50;
% calcVol50=handles.calc50;
% threshold=handles.threshk;
% sliceNo = handles.numSlices;
% info = handles.info;
% str_dir = string(handles.imageDirName);
% str_dir = [handles.imageDirName, 'result_image50.tif'];
% str_vol = [handles.imageDirName, 'result_vol50.txt'];
% fileID = fopen(str_vol,'w');
% fprintf(fileID,'Brain Volume : %3.2f cm^3\r\n',brainVol50);
% fprintf(fileID,'CSF Volume : %3.2f cm^3\r\n',fluidVol50);
% fprintf(fileID,'Subdural Volume : %3.2f cm^3\r\n',subVol50);
% fprintf(fileID,'Calcium Volume : %3.2f cm^3\r\n',calcVol50);
% fprintf(fileID,'Threshold : %3.4f \r\n',threshold);
% fclose(fileID);
% maxD = max(handles.data(:));
% minD = min(handles.data(:));
% dataN = (data - minD)/(maxD - minD);
% for i = 1:sliceNo
%     storeFiles = zeros(512, 1050, 3);
%     D = zeros(512,512,3);
%     D1 = zeros(512,512,3);
%     seg = segmentedImage(:,:,i);
%     tempData = dataN(:,:,i);
%     BinTissue = seg == 1;
%     BinFluid = seg == 2;
%     BinSubDural = seg == 3;
%     BinCalc = seg == 4;
%     rgbImage1 = cat(3, tempData, tempData, tempData);
%     storeFiles(1:512, 1:512,1:3) = rgbImage1;
%     D(:,:,2)=BinTissue;
%     D(:,:,1)=BinFluid;
%     D(:,:,3)=BinSubDural;
%     %D(:,:,1) = BinCalc;
%     D1 = .25*D + rgbImage1;
%     D22=D1(:,:,2);
%     D11=D1(:,:,1);
%     D33=D1(:,:,3);
%     D22(BinCalc==1)=1;
%     D11(BinCalc==1)=1;
%     D33(BinCalc==1)=0;
%     %D1(:,:,3)=.25*D(:,:,3) + rgbImage1(:,:,3);
%     D1(:,:,1)=D11;
%     D1(:,:,2)=D22;
%     D1(:,:,3)=D33;
%     storeFiles(1:512, 539:1050,1:3) = D1;
%     imwrite(storeFiles, str_dir, 'WriteMode', 'append');
% end
     %%%%REMOVE FOR FINAL%%%%%
     
%% saving volumes
% if(handles.Gender == 'M')
%     plotGrowthCurveMale(handles.brainVol, handles.Age);
%     saveas(gcf,str_growth)
% else
%     plotGrowthCurveFemale(handles.brainVol, handles.Age);
%     saveas(gcf,str_growth)
% end

%append_pdfs('final.pdf', input_slice{:});
close(wb);

resultVar = {'segmentedImage', 'data','dataOriginal', 'fluidVol', 'brainVol','info','subVol','calcVol','threshold'};
cd(handles.imageDirName);
% Use uisave to save results in user-specified directory:
uisave(resultVar,[handles.info.PatientID,'segmentation_results'])

cd(handles.code_dirname);

guidata(hObject,handles)


function volumeDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to volumeDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of volumeDisplay as text
%        str2double(get(hObject,'String')) returns contents of volumeDisplay as a double


% --- Executes during object creation, after setting all properties.
function volumeDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volumeDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in segSliceNo.
function segSliceNo_Callback(hObject, eventdata, handles)
% hObject    handle to segSliceNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns segSliceNo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from segSliceNo
load('segNet_train_test');
set(handles.axes2,'visible','on') 
set(handles.axes3,'visible','on') 
cla(handles.axes1)
set(handles.axes1,'visible','off') 
val1 = get(handles.segSliceNo,'Value');% only brain is selected
if(val1 == 1)
    set(handles.MainBox,'String','Select a Slice');
    guidata(hObject, handles);
else
    handles.sliceNum = val1 - 1;
    str2 = ['Slice ',num2str(val1 - 1)];
    tempImage = handles.data(:,:,handles.sliceNum);
    candRegion = handles.candidateRegion(:,:,handles.sliceNum);
    if(handles.sliceNum < handles.subNo)
        handles.segmentedImage(:,:,handles.sliceNum) = kmSeg(tempImage, candRegion, handles.centersKm);
    else
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,handles.sliceNum));
        if(isempty(xCoord))
            handles.segmentedImage(:,:,handles.sliceNum) = zeros(512,512);
        else
            handles.segmentedImage(:,:,handles.sliceNum) = segFinal1(tempImage, candRegion, handles.thrs, handles.subLabel, net);
        end
    end
    axes(handles.axes2);
    dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.data(:,:,handles.sliceNum));
    imshow(handles.dataNormalized(:,:,handles.sliceNum),'Parent',handles.axes3);
    set(get(handles.axes3, 'xlabel'), 'string', 'Original Image') 
    set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image') 
    set(handles.MainBox,'String',str2);
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function segSliceNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segSliceNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.



% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3


% --- Executes on button press in convertSub.
% function convertSub_Callback(hObject, eventdata, handles)
% % hObject    handle to convertSub (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% cla('reset')
% axes(handles.axes2) 
% dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
% %handles.h = imfreehand(gca);
% %handles.BW = createMask(handles.h);
% handles.BW = roipoly();
% X = handles.segmentedImage(:,:,handles.sliceNum);
% X(handles.BW) = 1;
% handles.segmentedImage(:,:,handles.sliceNum) = X;
% dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
% guidata(hObject, handles);


% --- Executes on button press in convertSub.
function convertSub_Callback(hObject, eventdata, handles)
% hObject    handle to convertSub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.BW = roipoly();
X = handles.segmentedImage(:,:,handles.sliceNum);
C = handles.candidateRegion(:,:,handles.sliceNum);
set(handles.Undo,'Enable','On')
handles.undoSliceNum = handles.sliceNum;
handles.undoSegmented = X;
X(handles.BW & C == 1) = 3;
handles.segmentedImage(:,:,handles.sliceNum) = X;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);


% --- Executes on button press in Undo.
function Undo_Callback(hObject, eventdata, handles)
% hObject    handle to Undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
sliceNo = handles.undoSliceNum;
dataPlotAll(handles.segmentedImage(:,:,sliceNo), handles.dataNormalized(:,:,sliceNo));
handles.segmentedImage(:,:,sliceNo) = handles.undoSegmented;
dataPlotAll(handles.segmentedImage(:,:,sliceNo), handles.dataNormalized(:,:,sliceNo));
set(handles.Undo,'Enable','Off')
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);


% --- Executes on button press in Instruct.
function Instruct_Callback(hObject, eventdata, handles)
% hObject    handle to Instruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% load('buttonIcon');
% % x = createButtonLines(Strings,...
% %                 {'fontsize',8,'fontweight','b'});
%set(handles.Instruct , 'Cdata' , off);            
instructions;
guidata(hObject, handles);


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function text6_CreateFcn(hObject, eventdata, handles)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text6.
function text6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
instructions;
guidata(hObject, handles);



%function genderBox_Callback(hObject, eventdata, handles)
% hObject    handle to genderBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of genderBox as text
%        str2double(get(hObject,'String')) returns contents of genderBox as a double


% --- Executes during object creation, after setting all properties.
%function genderBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to genderBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    %set(hObject,'BackgroundColor','white');
%end



function ageBox_Callback(hObject, eventdata, handles)
% hObject    handle to ageBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ageBox as text
%        str2double(get(hObject,'String')) returns contents of ageBox as a double


% --- Executes during object creation, after setting all properties.
function ageBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ageBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in regionBrain.
function regionBrain_Callback(hObject, eventdata, handles)
% hObject    handle to regionBrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
[sizeX, sizeY] = size(handles.data(:,:,handles.sliceNum));
[x, y] = getpts;
x = round(x);
y = round(y);
wb=waitbar(0,'Manual operation in progress...');
%reg = [x-2:x+2, y-2:y+2];
handles.BW  = zeros(sizeX, sizeY);
handles.BW(y-2:y+2, x-2:x+2) = 1;
X = handles.segmentedImage(:,:,handles.sliceNum);

set(handles.Undo,'Enable','On')
handles.undoSliceNum = handles.sliceNum;
handles.undoSegmented = X;
%X(handles.BW) = 1;
D = handles.data(:,:,handles.sliceNum);
C = handles.candidateRegion(:,:,handles.sliceNum);
handles.segmentedImage(:,:,handles.sliceNum) = regionGrowBrain_1(D, X, C, handles.BW );
close(wb);
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);


% --- Executes on button press in regionCSF.
function regionCSF_Callback(hObject, eventdata, handles)
% hObject    handle to regionCSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
[sizeX, sizeY] = size(handles.data(:,:,handles.sliceNum));
[x, y] = getpts;
x = round(x);
y = round(y);
wb=waitbar(0,'Manual operation in progress...');
%reg = [x-2:x+2, y-2:y+2];
handles.BW  = zeros(sizeX, sizeY);
handles.BW(y-2:y+2, x-2:x+2) = 1;
X = handles.segmentedImage(:,:,handles.sliceNum);

set(handles.Undo,'Enable','On')
handles.undoSliceNum = handles.sliceNum;
handles.undoSegmented = X;
%X(handles.BW) = 1;
D = handles.data(:,:,handles.sliceNum);
C = handles.candidateRegion(:,:,handles.sliceNum);
handles.segmentedImage(:,:,handles.sliceNum) = regionGrowFluid_1(D, X, C, handles.BW );
close(wb);
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);

% --- Executes on button press in regionSub.
function regionSub_Callback(hObject, eventdata, handles)
% hObject    handle to regionSub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
[sizeX, sizeY] = size(handles.data(:,:,handles.sliceNum));
[x, y] = getpts;
x = round(x);
y = round(y);
wb=waitbar(0,'Manual operation in progress...');
%reg = [x-2:x+2, y-2:y+2];
handles.BW  = zeros(sizeX, sizeY);
handles.BW(y-2:y+2, x-2:x+2) = 1;
X = handles.segmentedImage(:,:,handles.sliceNum);

set(handles.Undo,'Enable','On')
handles.undoSliceNum = handles.sliceNum;
handles.undoSegmented = X;
%X(handles.BW) = 1;
D = handles.data(:,:,handles.sliceNum);
C = handles.candidateRegion(:,:,handles.sliceNum);
handles.segmentedImage(:,:,handles.sliceNum) = regionGrowSub_1(D, X, C, handles.BW );
close(wb);
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);


% --- Executes on button press in regionInclude.
function regionInclude_Callback(hObject, eventdata, handles)
% hObject    handle to regionInclude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
[sizeX, sizeY] = size(handles.data(:,:,handles.sliceNum));
[x, y] = getpts;
x = round(x);
y = round(y);
wb=waitbar(0,'Manual operation in progress...');
%reg = [x-2:x+2, y-2:y+2];
handles.BW  = zeros(sizeX, sizeY);
handles.BW(y-2:y+2, x-2:x+2) = 1;
X = handles.segmentedImage(:,:,handles.sliceNum);

set(handles.Undo,'Enable','On')
handles.undoSliceNum = handles.sliceNum;
handles.undoSegmented = X;
%X(handles.BW) = 1;
D = handles.data(:,:,handles.sliceNum);
C = handles.candidateRegion(:,:,handles.sliceNum);
handles.segmentedImage(:,:,handles.sliceNum) = regionGrowInclude(D, X, C, handles.BW, handles.centersKm);
close(wb);
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);


% --- Executes on button press in instructRegionGrow.
function instructRegionGrow_Callback(hObject, eventdata, handles)
% hObject    handle to instructRegionGrow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
instructions_regionGrow;
guidata(hObject, handles);


% --- Executes on button press in calcium_seg.
function calcium_seg_Callback(hObject, eventdata, handles)
% hObject    handle to calcium_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
% 
% X = handles.segmentedImage(:,:,handles.sliceNum);
% set(handles.Undo,'Enable','On')
% handles.undoSliceNum = handles.sliceNum;
% handles.undoSegmented = X;

%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.undoSegmented2=handles.segmentedImage;
handles.BW2 = zeros(handles.rowNum,handles.colNum,handles.numSlices);
handles.Calc= zeros(handles.rowNum,handles.colNum,handles.numSlices);
wait=0;
wb=waitbar(0,'Calcium segmentation in process...');
for slice = 1:handles.numSlices
    CR=zeros(512,512);
    CR2=zeros(512,512);
    BB=zeros(512,512);
    waitbar(wait/(handles.numSlices + 2),wb);
	wait = wait + 1;
    handles.BR(:,:,slice)=handles.segmentedImage(:,:,slice);
    CR(handles.BR(:,:,slice)~=0)=1;
%     se = strel('cube',5);
%     e = imerode(CR,se);
%     ZZ=imfill(e,4,'holes'); 
%     D = bwdist(ZZ);
%     CR(D>0)=0;
    %se = strel('disk',2);
    %handles.BW2(:,:,slice) = imerode(CR,se);
    handles.BW2(:,:,slice) = imfill(CR,4,'holes'); 
    BB(handles.BW2(:,:,slice)~=0)=1;
    se = strel('cube',8);
    E = imerode(BB,se);
%   handles.Calc(:,:,slice)= handles.BW2(:,:,slice)~=CR;
    Original=handles.dataOrig(:,:,slice);
    %C=handles.Calc(:,:,slice);
    D=handles.segmentedImage(:,:,slice);
%     D(C==1)=4;
%     handles.segmentedImage(:,:,slice)=D;
    CR2(Original>=50 & E==1)=1;
    D(CR2==1)=4;
    handles.segmentedImage(:,:,slice)=D;
    %%%%REMOVE FOR FINAL%%%%%
%     CR3=zeros(512,512);
%     D=handles.BR(:,:,slice);
%     CR3(Original>=35 & E==1)=1;
%     D(CR3==1)=4;
%     handles.segmentedImage35(:,:,slice)=D;
%     CR4=zeros(512,512);
%     D=handles.BR(:,:,slice);
%     CR4(Original>=45 & E==1)=1;
%     D(CR4==1)=4;
%     handles.segmentedImage45(:,:,slice)=D;
%     CR5=zeros(512,512);
%     D=handles.BR(:,:,slice);
%     CR5(Original>=50 & E==1)=1;
%     D(CR5==1)=4;
%     handles.segmentedImage50(:,:,slice)=D;
     %%%%REMOVE FOR FINAL%%%%%
end
close(wb);
set(handles.calc_undo,'visible', 'On')
set(handles.calc_undo,'Enable', 'On')
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);


% --- Executes on button press in calcium_roi.
function calcium_roi_Callback(hObject, eventdata, handles)
% hObject    handle to calcium_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.BW = roipoly();
X = handles.segmentedImage(:,:,handles.sliceNum);
C = handles.candidateRegion(:,:,handles.sliceNum);
set(handles.Undo,'Enable','On')
handles.undoSliceNum = handles.sliceNum;
handles.undoSegmented = X;
X(handles.BW & C == 1) = 4;
handles.segmentedImage(:,:,handles.sliceNum) = X;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);


% --- Executes on button press in threshold_button.
function threshold_button_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.thresh_num,'visible', 'On')
set(handles.threshold_slider,'visible', 'On')
set(handles.threshold_slider,'Enable', 'On')
set(handles.axes5,'visible', 'On')
set(handles.thresh_undo,'visible', 'On')
set(handles.thresh_undo,'Enable','On')
sliceNo = handles.numSlices;
sliceNum = handles.numSlices;
handles.newSeg=zeros(512,512,handles.numSlices);
handles.undoSegmented3=handles.segmentedImage;
for i=1:sliceNo
    %NS=zeros(512,512);
    %NSn=zeros(512,512);
    %NS=handles.data(:,:,i);
    origD=handles.dataOrig(:,:,i);
    origD(origD>50 | origD<=0)=0;
    %NSn(origD<=50 & origD>=0)=origD;
    %NSn=handles.dataNormalized(:,:,i);
    CR=handles.candidateRegion(:,:,i);
    %SRI=handles.segmentedImage(:,:,i);
    %NS(CR==0)=0;
    %NS(SRI==3)=0;
    %NSn(CR==0)=0;
    origD(CR==0)=0;
    %handles.newSeg(:,:,i)=NS;
    %handles.newSeg2(:,:,i)=NSn;
    handles.newSeg2(:,:,i)=origD;
end
%handles.A=reshape(handles.newSeg,[],1);
%[handles.Ka,handles.cka]=kmeans(handles.A,2);
handles.A2=reshape(handles.newSeg2,[],1);
handles.A2=nonzeros(handles.A2);
[handles.Ka2,handles.cka2]=kmeans(handles.A2,2);
h=histogram(handles.newSeg2);
[maxcount, whichbin] = max(h.Values);
handles.C_thresh=(h.BinWidth*whichbin);
%handles.newSeg=uint8(handles.newSeg);
handles.threshk=(handles.cka2(1)+handles.cka2(2))/2;
handles.L_adjusted=zeros(512,512,sliceNo);
for i=1:sliceNo
    NSnew=zeros(512,512);
    fillI=zeros(512,512);
%     CR=handles.candidateRegion(:,:,i);
%     CR=imbinarize(CR);
%     se = strel('cube',1);
%     e = imerode(CR,se);
%     ZZ=imfill(e,4,'holes'); 
%     D = bwdist(ZZ);
%     CR(D>0)=0;
%     D=handles.candidateRegion(:,:,i);
%     D(CR==0)=0;
%     handles.candidateRegion(:,:,i)=D;
    %SubKDN=handles.dataNormalized(:,:,i);
    %SubKDN=handles.dataOrig(:,:,i);
    SubKDN=handles.newSeg2(:,:,i);
    SubK=handles.segmentedImage(:,:,i);
    subCR=handles.candidateRegion(:,:,i);
    %Kseg=handles.Lnew(:,:,i);
    %NSnew(SubKDN>handles.threshk & SubKDN>handles.C_thresh & SubKDN<1 & subCR==1)=1;
    %NSnew(SubKDN<=handles.threshk & SubKDN>handles.C_thresh & SubKDN<1 & subCR==1)=2;
    NSnew(SubKDN>handles.threshk & subCR==1 & SubK~=0)=1;
    NSnew(SubKDN<=handles.threshk & subCR==1 & SubK~=0)=2;
    fillI(NSnew==2)=1;
    fillI=imfill(fillI);
    NSnew(fillI==1)=2;
    %NSnew(SubK==0)=0;
    %NSnew(Kseg==3)=1;
    NSnew(SubK==3)=3;
    NSnew(SubK==4)=4;
    handles.L_adjusted(:,:,i)=NSnew;
end
%[handles.Lnew,handles.centersnew] = imsegkmeans3(handles.newSeg,4);
% handles.L_adjusted=zeros(512,512,sliceNo);
% for i=1:sliceNo
%     NSnew=zeros(512,512);
%     SubK=handles.segmentedImage(:,:,i);
%     Kseg=handles.Lnew(:,:,i);
%     NSnew(Kseg==1)=0;
%     NSnew(Kseg==2)=2;
%     NSnew(Kseg==3)=0;
%     NSnew(Kseg==4)=1;
%    %NSnew(Kseg==3)=1;
%     NSnew(SubK==3)=3;
%     NSnew(SubK==4)=4;
%     handles.L_adjusted(:,:,i)=NSnew;
% end
axes(handles.axes5); % Switch current axes to axes11.
cla
handles.Hist=histogram(handles.newSeg2);
ylim([0 200000]);
%xlim([0 1]);
xlim([1 50]);
hold on;
%line([handles.centersKm(1)/100, handles.centersKm(1)/100], ylim, 'LineWidth', 2, 'Color', 'r');
line([handles.threshk, handles.threshk], ylim, 'LineWidth', 2, 'Color', 'r');
axes(handles.axes2);
handles.segmentedImage=handles.L_adjusted;
set(handles.threshold_slider, 'Max',50);
set(handles.threshold_slider, 'Min', 0);
set(handles.threshold_slider, 'Value', handles.threshk);
set(handles.thresh_num,'String',sprintf(num2str(handles.threshk)));
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function threshold_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
 set(handles.threshold_slider, 'Max',50);
 set(handles.threshold_slider, 'Min', -15);
 set(handles.threshold_slider, 'SliderStep' , [2,2] );
 %handles.V1=round(handles.threshk*100)/100;
 handles.V1=round(handles.threshk);
 set(handles.threshold_slider, 'Value', handles.V1);
 guidata(hObject, handles);
 
% --- Executes on slider movement.
function threshold_slider_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sliderValue = get(handles.threshold_slider,'Value');
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliceNo = handles.numSlices;
sliceNum = handles.numSlices;
% handles.newSeg=zeros(512,512,handles.numSlices);
% for i=1:sliceNo
%     %NS=zeros(512,512);
%     NSn=zeros(512,512);
%     %NS=handles.data(:,:,i);
%     NSn=handles.dataNormalized(:,:,i);
%     CR=handles.candidateRegion(:,:,i);
%     SRI=handles.segmentedImage(:,:,i);
%     %NS(CR==0)=0;
%     %NS(SRI==3)=0;
%     NSn(CR==0)=0;
%     NSn(CR==0)=0;
%     %handles.newSeg(:,:,i)=NS;
%     handles.newSeg2(:,:,i)=NSn;
% end
% %handles.A=reshape(handles.newSeg,[],1);
% %[handles.Ka,handles.cka]=kmeans(handles.A,2);
% handles.A2=reshape(handles.newSeg2,[],1);
% handles.A2=nonzeros(handles.A2);
% [handles.Ka2,handles.cka2]=kmeans(handles.A2,2);
%handles.newSeg=uint8(handles.newSeg);
handles.threshk=handles.sliderValue;
handles.L_adjusted=zeros(512,512,sliceNo);
for i=1:sliceNo
    NSnew=zeros(512,512);
    fillI=zeros(512,512);
%     CR=handles.candidateRegion(:,:,i);
%     CR=imbinarize(CR);
%     se = strel('cube',1);
%     e = imerode(CR,se);
%     ZZ=imfill(e,4,'holes'); 
%     D = bwdist(ZZ);
%     CR(D>0)=0;
%     D=handles.candidateRegion(:,:,i);
%     D(CR==0)=0;
%     handles.candidateRegion(:,:,i)=D;
        %SubKDN=handles.dataOrig(:,:,i);
    SubKDN=handles.newSeg2(:,:,i);
    SubK=handles.segmentedImage(:,:,i);
    subCR=handles.candidateRegion(:,:,i);
    %Kseg=handles.Lnew(:,:,i);
    %NSnew(SubKDN>handles.threshk & SubKDN>handles.C_thresh & SubKDN<1 & subCR==1)=1;
    %NSnew(SubKDN<=handles.threshk & SubKDN>handles.C_thresh & SubKDN<1 & subCR==1)=2;
    NSnew(SubKDN>handles.threshk & subCR==1 & SubK~=0)=1;
    NSnew(SubKDN<=handles.threshk & subCR==1 & SubK~=0)=2;
    fillI(NSnew==2)=1;
    fillI=imfill(fillI);
    NSnew(fillI==1)=2;
    %NSnew(SubK==0)=0;
    %NSnew(Kseg==3)=1;
    NSnew(SubK==3)=3;
    NSnew(SubK==4)=4;
    handles.L_adjusted(:,:,i)=NSnew;
end
%[handles.Lnew,handles.centersnew] = imsegkmeans3(handles.newSeg,4);
% handles.L_adjusted=zeros(512,512,sliceNo);
% for i=1:sliceNo
%     NSnew=zeros(512,512);
%     SubK=handles.segmentedImage(:,:,i);
%     Kseg=handles.Lnew(:,:,i);
%     NSnew(Kseg==1)=0;
%     NSnew(Kseg==2)=2;
%     NSnew(Kseg==3)=0;
%     NSnew(Kseg==4)=1;
%    %NSnew(Kseg==3)=1;
%     NSnew(SubK==3)=3;
%     NSnew(SubK==4)=4;
%     handles.L_adjusted(:,:,i)=NSnew;
% end
axes(handles.axes5); % Switch current axes to axes11.
cla
handles.Hist=histogram(handles.newSeg2);
ylim([0 200000]);
%xlim([0 1]);
xlim([1 50]);
hold on;
%line([handles.centersKm(1)/100, handles.centersKm(1)/100], ylim, 'LineWidth', 2, 'Color', 'r');
line([handles.threshk, handles.threshk], ylim, 'LineWidth', 2, 'Color', 'r');
axes(handles.axes2);
handles.segmentedImage=handles.L_adjusted;
set(handles.thresh_num,'String',sprintf(num2str(handles.threshk)));
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
guidata(hObject, handles);

function thresh_num_Callback(hObject, eventdata, handles)
% hObject    handle to thresh_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh_num as text
%        str2double(get(hObject,'String')) returns contents of thresh_num as a double


% --- Executes during object creation, after setting all properties.
function thresh_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in skull_calc_seg.
function skull_calc_seg_Callback(hObject, eventdata, handles)
% hObject    handle to skull_calc_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.BW = roipoly();
X = handles.segmentedImage(:,:,handles.sliceNum);
C = handles.candidateRegion(:,:,handles.sliceNum);
Data = handles.dataNormalized(:,:,handles.sliceNum);
Original =handles.dataOrig(:,:,handles.sliceNum);
set(handles.Undo,'Enable','On')
handles.undoSliceNum = handles.sliceNum;
handles.undoSegmented = X;
%X(handles.BW ==1 & Data >= 0.7) = 4;
X(handles.BW ==1 & Original >= 45) = 4;
handles.segmentedImage(:,:,handles.sliceNum) = X;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);


% --- Executes on button press in shunt.
function shunt_Callback(hObject, eventdata, handles)
% hObject    handle to shunt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla('reset')
axes(handles.axes2) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.BW = roipoly();
X = handles.segmentedImage(:,:,handles.sliceNum);
set(handles.Undo,'Enable','On')
handles.undoSliceNum = handles.sliceNum;
handles.undoSegmented = X;
X(handles.BW ==1 & X == 4) = 0;
handles.segmentedImage(:,:,handles.sliceNum) = X;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image');
guidata(hObject, handles);


% --- Executes on button press in editor.
function editor_Callback(hObject, eventdata, handles)
% hObject    handle to editor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes5)
cla
set(handles.volumeDisplay, 'String', ''); 
set(handles.thresh_num,'visible', 'Off')
set(handles.threshold_slider,'visible', 'Off')
set(handles.threshold_slider,'Enable', 'Off')
set(handles.saveResults,'Enable', 'Off')
%% start the code
[image_filename,image_dirname]=uigetfile(('*.mat'),'Select DICOM image from desired image stack');%selecting the directory
% Load image stack from selected image set:
load(fullfile(image_dirname, image_filename));
info = info;
handles.info = info;
set(handles.File, 'String',['Patient ID:  ', handles.info.PatientID]); 
handles.segmentedImage = segmentedImage;
handles.data = data;
handles.dataOrig=dataOriginal;
cand = segmentedImage > 0;
handles.candidateRegion = cand;
[handles.rowNum, handles.colNum, handles.numSlices] = size(handles.segmentedImage);
handles.imageDirName = image_dirname;
handles.dataNormalized = (handles.data + 50)/150;
%axes(handles.axes1)
str1{1} = 'Select slices...';
for p = 1:handles.numSlices
    str1{p + 1} = ['Slice ',num2str(p)];
end
set(handles.selectSlice,'String',str1);
%set(handles.startSlice,'String',str1);
%set(handles.subduralSlice,'String',str1);
%set(handles.segSliceNo,'String',str1);
%slide show of images
% for p = 1 : handles.numSlices
%     handles.sliceNum = p;
%     imshow(handles.dataNormalized(:,:,handles.sliceNum),'Parent',handles.axes1);
%     str2 = ['Slice ',num2str(p)];
%     set(handles.MainBox,'String',str2);
%     set(handles.selectSlice,'Value',p+1);
%     pause(1);
% end
str2 = ['Slice ',num2str(1)];
set(handles.MainBox,'String',str2);
set(handles.selectSlice,'Value',2);
%% finding candidate region

%% finding Kmean centers
temp = handles.data;
L = handles.candidateRegion;
%L = (temp > thrsKm);
[~, centersKm] = kmeans(temp(L == 1), 2);
handles.centersKm = centersKm;
handles.flagSegFinish = 1;
cla(handles.axes1)
set(handles.axes1,'visible','off')
handles.sliceNum = 1;
set(handles.axes2,'visible','on') 
set(handles.axes3,'visible','on') 

cla(handles.axes1)
set(handles.axes1,'visible','off')
handles.sliceNum = 1;
set(handles.axes2,'visible','on') 
set(handles.axes3,'visible','on') 
%set(handles.displaySeg,'Value', 1);
axes(handles.axes2)
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.data(:,:,handles.sliceNum));
imshow(handles.dataNormalized(:,:,handles.sliceNum),'Parent',handles.axes3);
set(get(handles.axes3, 'xlabel'), 'string', 'Original Image') 
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image') 
set(handles.MainBox,'String','Segmentation Complete!');
set(handles.deleteRegion,'Enable','On')
set(handles.Instruct,'Enable', 'On')
set(handles.addRegion,'Enable','On')
set(handles.convertBrain,'Enable','On')
set(handles.calculateVolumes,'Enable','On')
set(handles.convertFluid,'Enable','On')
set(handles.convertSub,'Enable','On')
set(handles.saveResults,'Enable','On')
set(handles.calcium_seg,'Enable', 'On')
set(handles.calcium_roi,'Enable', 'On')
set(handles.threshold_button,'Enable', 'On')
set(handles.threshold_slider,'Enable', 'On')
set(handles.thresh_num,'visible', 'Off')
set(handles.axes5,'visible', 'Off')
set(handles.skull_calc_seg,'Enable', 'On')
set(handles.shunt,'Enable', 'On')
set(handles.shunt,'visible', 'On')
set(handles.regionBrain,'Enable','On')
set(handles.regionCSF,'Enable','On')
set(handles.regionSub,'Enable','On')
set(handles.regionInclude,'Enable','On')
set(handles.instructRegionGrow,'Enable', 'On')
set(handles.brainLeg,'visible', 'On')
set(handles.csfLeg,'visible', 'On')
set(handles.subLeg,'visible', 'On')
set(handles.calcLeg,'visible', 'On')
set(handles.calc_undo,'visible', 'Off')
set(handles.thresh_undo,'visible', 'Off')
guidata(hObject, handles);


% --- Executes on button press in calc_undo.
function calc_undo_Callback(hObject, eventdata, handles)
% hObject    handle to calc_undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
handles.segmentedImage = handles.undoSegmented2;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(handles.calc_undo,'Enable','Off')
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image')
guidata(hObject, handles);


% --- Executes on button press in thresh_undo.
function thresh_undo_Callback(hObject, eventdata, handles)
% hObject    handle to thresh_undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes2) 
handles.segmentedImage = handles.undoSegmented3;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
set(handles.thresh_undo,'Enable','Off')
set(handles.thresh_num,'visible', 'Off')
axes(handles.axes5); % Switch current axes to axes11.
cla
set(handles.axes5,'visible', 'Off')
set(handles.threshold_slider,'visible', 'Off')
set(get(handles.axes2, 'xlabel'), 'string', 'Segmented Image')
guidata(hObject, handles);



function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of File as text
%        str2double(get(hObject,'String')) returns contents of File as a double


% --- Executes during object creation, after setting all properties.
function File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
