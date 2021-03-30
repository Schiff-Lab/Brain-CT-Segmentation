function varargout = GUI_FLIS(varargin)
% clc;
% clear global variables;
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

% Last Modified by GUIDE v2.5 04-Apr-2017 10:08:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_FLIS_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_FLIS_OutputFcn, ...
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
function GUI_FLIS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FLIS (see VARARGIN)

% Choose default command line output for GUI_FLIS
handles.output = hObject;
handles.flagStartSlice = 0;
handles.flagSubDuralSlice = 0;
handles.flagDisplaySeg = 0;
set(handles.FLIS,'Enable','Off')
set(handles.deleteRegion,'Enable','Off')
set(handles.addRegion,'Enable','Off')
set(handles.convertBrain,'Enable','Off')
set(handles.calculateVolumes,'Enable','Off')
set(handles.convertFluid,'Enable','Off')
set(handles.saveResults,'Enable','Off')
%set(handles.volumeDisplay,'Enable','Off')
handles.code_dirname = cd;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_FLIS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_FLIS_OutputFcn(hObject, eventdata, handles) 
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
[image_filename,image_dirname]=uigetfile(('*.ima;*.dcm'),'Select DICOM image from desired image stack');%selecting the directory
% Load image stack from selected image set:
handles.files = [dir(fullfile(image_dirname,'/*.ima'));dir(fullfile(image_dirname,'/*.dcm'))];
handles.info=dicominfo(fullfile(image_dirname, handles.files(1).name));
handles.imageDirName = image_dirname;
handles.rowNum = cast(handles.info.Rows,'double');
handles.colNum = cast(handles.info.Columns,'double');
handles.numSlices = length(handles.files);
handles.data = zeros(handles.rowNum,handles.colNum,handles.numSlices);
handles.dataOrig = zeros(handles.rowNum,handles.colNum,handles.numSlices);
handles.dataNormalized = zeros([handles.rowNum,handles.colNum,handles.numSlices]);
for p = 1:handles.numSlices
  handles.data(:,:,p) = dicomread(fullfile(image_dirname,handles.files(p).name));
  handles.dataOrig(:,:,p) = handles.data(:,:,p);
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
set(handles.startSlice,'String',str1);
set(handles.subduralSlice,'String',str1);
% 
handles.sliceNum = 1;
imshow(handles.dataNormalized(:,:,handles.sliceNum),'Parent',handles.axes1);
str2 = ['Slice ',num2str(1)];
set(handles.MainBox,'String',str2);
set(handles.selectSlice,'Value',2);
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

% --- Executes on button press in PreviousSlice.
function PreviousSlice_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = handles.sliceNum;
%val1 = get(handles.selectSlice,'Value');
str2 = ['Slice ',num2str(val - 1)];
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
for slice = firstPartition + ceil((firstTransition - firstPartition)/2) + 1:middleTrans + ceil((secondTransition - middleTrans)/2)
    handles.candidateRegion(:,:,slice) = findBinaryMiddleSlice_1(handles.data(:,:,slice), meanX);
end
for slice = firstPartition + ceil((firstTransition - firstPartition)/2):-1:1
    handles.candidateRegion(:,:,slice) = crossSliceBinary(handles.data(:,:,slice), handles.candidateRegion(:,:,slice + 1));
end
for slice = middleTrans + ceil((secondTransition - middleTrans)/2) + 1:handles.numSlices
    handles.candidateRegion(:,:,slice) = crossSliceBinaryNew(handles.data(:,:,slice), handles.candidateRegion(:,:,slice - 1), meanX);
end
%% removing Non Brain
for slice = firstPartition:-1:1
    handles.candidateRegion(:,:,slice) = removeNonBrain(handles.data(:,:,slice), handles.candidateRegion(:,:,slice));
end

%% starting segmentation 
%% start segmenting the data
start_spams
load('testDictksvdDFDL_1');
load('meanVal');
%normalizing dictionaries
%% finding Kmean centers
temp = handles.data;
minInt = min(temp(:));
R = temp(handles.candidateRegion == 1);
thrsKm = min(temp(handles.candidateRegion == 1));
thrsKm  = max(thrsKm, 1/5*minInt);
L = (temp > thrsKm);
[~, centersKm] = kmeans(temp(handles.candidateRegion == 1), 2);
[centersKm, ~] = sort(centersKm);
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
dict = [dictBrain{6} dictFluid{6} dictSubDural{6}];
w = [wBrain{6} wFluid{6} wSubDural{6}];
meanDist = meanNormDist{6};
meanInt = meanNormInt{6};
%[centers, subLabel] = segmentationKsvdDFDL_findPar_final(tempImage, candRegion, dict,w,meanInt,meanDist);
subLabelNum = handles.subNo;
segNo = handles.segNo;
subLabel = 2;
centers = [];
%% segmentation for all images
%segments data partition wise
    for i = 1:ceil(firstPartition/2)
        if(i < segNo)
           handles.segmentedImage(:,:,i) = 0;
        else
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, centersKm);
        else
        dict = [dictBrain{2} dictFluid{2} dictSubDural{2}];
        w = [wBrain{2} wFluid{2} wSubDural{2}];
        meanDist = meanNormDist{1};
        meanInt = meanNormInt{1};
        handles.segmentedImage(:,:,i) = segmentationFLIS(tempImage, candRegion, dict,w,meanInt,meanDist, centers, subLabel);
        end
        end
     end
    %creating second partition
    for i = ceil(firstPartition/2) + 1:firstPartition
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, centersKm);
        else
        dict = [dictBrain{2} dictFluid{2} dictSubDural{2}];
        w = [wBrain{2} wFluid{2} wSubDural{2}];
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            meanDist = meanNormDist{2};
            meanInt = meanNormInt{2};
            handles.segmentedImage(:,:,i) = segmentationFLIS(tempImage, candRegion, dict,w,meanInt,meanDist, centers, subLabel);
        end
        end
    end
    %creating third partition
    for i = firstPartition + 1:firstPartition + ceil((firstTransition - firstPartition)/2)
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, centersKm);
        else
        dict = [dictBrain{3} dictFluid{3} dictSubDural{3}];
        w = [wBrain{3} wFluid{3} wSubDural{3}];
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            meanDist = meanNormDist{3};
            meanInt = meanNormInt{3};
            handles.segmentedImage(:,:,i) = segmentationFLIS(tempImage, candRegion, dict,w,meanInt,meanDist, centers, subLabel);
        end
        end
    end
    %creating fourth partition
    for i = firstPartition + ceil((firstTransition - firstPartition)/2) + 1:firstTransition
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, centersKm);
        else
        dict = [dictBrain{4} dictFluid{4} dictSubDural{4}];
        w = [wBrain{4} wFluid{4} wSubDural{4}];
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            meanDist = meanNormDist{4};
            meanInt = meanNormInt{4};
            handles.segmentedImage(:,:,i) = segmentationFLIS(tempImage, candRegion, dict,w,meanInt,meanDist, centers, subLabel);
        end
        end
    end
    %creating fifth partition
    for i = firstTransition + 1:firstTransition + ceil((transMiddle - firstTransition)/2)
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, centersKm);
        else
        dict = [dictBrain{5} dictFluid{5} dictSubDural{5}];
        w = [wBrain{5} wFluid{5} wSubDural{5}];
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            meanDist = meanNormDist{5};
            meanInt = meanNormInt{5};
            handles.segmentedImage(:,:,i) = segmentationFLIS(tempImage, candRegion, dict,w,meanInt,meanDist, centers, subLabel);
        end
        end
    end
    %creating sixth partition
    for i = firstTransition + ceil((transMiddle - firstTransition)/2) + 1:transMiddle
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, centersKm);
        else
        dict = [dictBrain{6} dictFluid{6} dictSubDural{6}];
        w = [wBrain{6} wFluid{6} wSubDural{6}];
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            meanDist = meanNormDist{6};
            meanInt = meanNormInt{6};
            handles.segmentedImage(:,:,i) = segmentationFLIS(tempImage, candRegion, dict,w,meanInt,meanDist, centers, subLabel);
        end
        end
    end
    %creating seventh partition
    for i = transMiddle + 1:transMiddle + ceil((middleTrans - transMiddle)/2)
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, centersKm);
        else
        dict = [dictBrain{6} dictFluid{6} dictSubDural{6}];
        w = [wBrain{6} wFluid{6} wSubDural{6}];
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            meanDist = meanNormDist{6};
            meanInt = meanNormInt{6};
            handles.segmentedImage(:,:,i) = segmentationFLIS(tempImage, candRegion, dict,w,meanInt,meanDist, centers, subLabel);
        end
        end
    end
    %creating eigth partition
    for i = transMiddle + ceil((middleTrans - transMiddle)/2) + 1:middleTrans
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, centersKm);
        else
        dict = [dictBrain{8} dictFluid{8} dictSubDural{8}];
        w = [wBrain{8} wFluid{8} wSubDural{8}];
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            meanDist = meanNormDist{8};
            meanInt = meanNormInt{8};
            handles.segmentedImage(:,:,i) = segmentationFLIS(tempImage, candRegion, dict,w,meanInt,meanDist, centers, subLabel);
        end
        end
    end
    %creating nineth partition
    for i = middleTrans + 1:middleTrans + ceil((secondTransition - middleTrans)/2)
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, centersKm);
        else
        dict = [dictBrain{9} dictFluid{9} dictSubDural{9}];
        w = [wBrain{9} wFluid{9} wSubDural{9}];
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            meanDist = meanNormDist{9};
            meanInt = meanNormInt{9};
            handles.segmentedImage(:,:,i) = segmentationFLIS(tempImage, candRegion, dict,w,meanInt,meanDist, centers, subLabel);
        end
        end
    end
    %creating tenth partition
    for i = middleTrans + ceil((secondTransition - middleTrans)/2) + 1:secondTransition
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, centersKm);
        else
        dict = [dictBrain{9} dictFluid{9} dictSubDural{9}];
        w = [wBrain{9} wFluid{9} wSubDural{9}];
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            meanDist = meanNormDist{9};
            meanInt = meanNormInt{9};
            handles.segmentedImage(:,:,i) = segmentationFLIS(tempImage, candRegion, dict,w,meanInt,meanDist, centers, subLabel);
        end
        end
    end
    %creating eleventh partition
    for i = secondTransition + 1:secondTransition + ceil((lastPartition - secondTransition)/2)
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, centersKm);
        else
        dict = [dictBrain{11} dictFluid{11} dictSubDural{11}];
        w = [wBrain{11} wFluid{11} wSubDural{11}];
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            meanDist = meanNormDist{11};
            meanInt = meanNormInt{11};
            handles.segmentedImage(:,:,i) = segmentationFLIS(tempImage, candRegion, dict,w,meanInt,meanDist, centers, subLabel);
        end
        end
    end
    %creating twelveth partition
    for i = secondTransition + ceil((lastPartition - secondTransition)/2) + 1:lastPartition
        tempImage = handles.data(:,:,i);
        candRegion = handles.candidateRegion(:,:,i);
        if(i < subLabelNum)
            handles.segmentedImage(:,:,i) = kmSeg(tempImage, candRegion, centersKm);
        else
        dict = [dictBrain{12} dictFluid{12} dictSubDural{12}];
        w = [wBrain{12} wFluid{12} wSubDural{12}];
        [xCoord, yCoord] = find(handles.candidateRegion(:,:,i));
        if(isempty(xCoord))
            continue;
        else
            meanDist = meanNormDist{12};
            meanInt = meanNormInt{12};
            handles.segmentedImage(:,:,i) = segmentationFLIS(tempImage, candRegion, dict,w,meanInt,meanDist, centers, subLabel);
        end
        end
    end
%% display segmented region
handles.sliceNum = 1;
%% display segmented slice
% greenMap = zeros(512,512,3);
% greenMap(:,:,2) = 1;
% redMap = zeros(512,512,3);
% redMap(:,:,1) = 1;
% blueMap = zeros(512,512,3);
% blueMap(:,:,3) = 1;
% BinTissue = handles.segmentedImage(:,:,handles.sliceNum) == 1;
% BinFluid = handles.segmentedImage(:,:,handles.sliceNum) == 2;
% BinSubDural = handles.segmentedImage(:,:,handles.sliceNum) >= 3;
% data = handles.data(:,:,handles.sliceNum);
% minD = min(data(:));
% maxD = max(data(:));
% imagesc(data, [minD maxD]);
% colormap(gray);
% hold on
% X = imshow(greenMap,'Parent',handles.axes1);
% set(X, 'AlphaData', BinTissue.*.3);
% hold on
% Y = imshow(redMap,'Parent',handles.axes1);
% set(Y, 'AlphaData', BinFluid.*.3);
% Y = imshow(blueMap,'Parent',handles.axes1);
% set(Y, 'AlphaData', BinSubDural.*.3);
%% 
set(handles.displaySeg,'Value', 1);
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.data(:,:,handles.sliceNum));
set(handles.MainBox,'String','Segmentation Complete!');
set(handles.deleteRegion,'Enable','On')
set(handles.addRegion,'Enable','On')
set(handles.convertBrain,'Enable','On')
set(handles.calculateVolumes,'Enable','On')
set(handles.convertFluid,'Enable','On')
set(handles.saveResults,'Enable','On')
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
axes(handles.axes1) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.BW = roipoly();
X = handles.segmentedImage(:,:,handles.sliceNum);
X(handles.BW) = 0;
handles.segmentedImage(:,:,handles.sliceNum) = X;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
guidata(hObject, handles);


% --- Executes on button press in addRegion.
function addRegion_Callback(hObject, eventdata, handles)
% hObject    handle to addRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes1) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.BW = roipoly();
X = handles.segmentedImage(:,:,handles.sliceNum);
[xCoord, yCoord] = find(handles.BW);
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
guidata(hObject, handles);

% --- Executes on button press in convertFluid.
function convertFluid_Callback(hObject, eventdata, handles)
% hObject    handle to convertFluid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes1) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.BW = roipoly();
X = handles.segmentedImage(:,:,handles.sliceNum);
X(handles.BW) = 2;
handles.segmentedImage(:,:,handles.sliceNum) = X;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
guidata(hObject, handles);

% --- Executes on button press in convertBrain.
function convertBrain_Callback(hObject, eventdata, handles)
% hObject    handle to convertBrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla('reset')
axes(handles.axes1) 
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
%handles.h = imfreehand(gca);
%handles.BW = createMask(handles.h);
handles.BW = roipoly();
X = handles.segmentedImage(:,:,handles.sliceNum);
X(handles.BW) = 1;
handles.segmentedImage(:,:,handles.sliceNum) = X;
dataPlotAll(handles.segmentedImage(:,:,handles.sliceNum), handles.dataNormalized(:,:,handles.sliceNum));
guidata(hObject, handles);

% --- Executes on button press in calculateVolumes.
function calculateVolumes_Callback(hObject, eventdata, handles)
% hObject    handle to calculateVolumes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X = handles.segmentedImage;
brainTissue = find(X == 1);
fluidCSF = find(X == 2);
voxelVol = (handles.info.PixelSpacing(1)*handles.info.PixelSpacing(2)*handles.info.SliceThickness);
handles.brainVol = length(brainTissue)*voxelVol/10^3;
handles.fluidVol = length(fluidCSF)*voxelVol/10^3;
set(handles.volumeDisplay,'String',['Tissue volume = ',num2str(handles.brainVol),' cm^3   ', '   Fluid volume = ',num2str(handles.fluidVol),' cm^3']);
guidata(hObject, handles);

% --- Executes on button press in saveResults.
function saveResults_Callback(hObject, eventdata, handles)
% hObject    handle to saveResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
segmentedImage = handles.segmentedImage;
data = handles.data;
fluidVol = handles.fluidVol;
brainVol = handles.brainVol;
resultVar = {'segmentedImage', 'data', 'fluidVol', 'brainVol'};
cd(handles.imageDirName);
% Use uisave to save results in user-specified directory:
uisave(resultVar,'segmentation_results')

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
