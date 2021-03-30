function [X, Y] = dataPlotAll(imBinary, data)
[indx, indy] = size(data);
greenMap = zeros(indx,indy,3);
greenMap(:,:,2) = 1;
redMap = zeros(indx,indy,3);
redMap(:,:,1) = 1;
blueMap = zeros(indx,indy,3);
blueMap(:,:,3) = 1;
yellowMap = zeros(indx,indy,3);
yellowMap(:,:,3) = 0;
yellowMap(:,:,1) = 1;
yellowMap(:,:,2) = 1;
BinTissue = imBinary == 1;
BinFluid = imBinary == 2;
BinSubDural = imBinary == 3;
BinCalc = imBinary == 4;
minD = min(data(:));
maxD = max(data(:));
imagesc(data, [minD maxD]);
colormap(gray);
hold on
X = imshow(greenMap);
set(X, 'AlphaData', BinTissue.*.3);
hold on
Y = imshow(redMap);
set(Y, 'AlphaData', BinFluid.*.3);
Y = imshow(blueMap);
set(Y, 'AlphaData', BinSubDural.*.3);
Y = imshow(yellowMap);
set(Y, 'AlphaData', BinCalc.*.3);