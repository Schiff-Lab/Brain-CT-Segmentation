

function FLISImage = segFinal1(data, cand, thrs, subLabel, net)

tempImage = data;
testCand = cand;
binImage = testCand; 
distBinImage = bwdist(~binImage); % distance transform to calculate distance values
%labelImage = zeros(512,512);
maxInt = max(tempImage(:));
R = tempImage(testCand == 1);
thrsSeg = mean(R) + 3*std(R); % making sure only candidate region is being segmented
thrsSeg = min(thrsSeg, 2/3*maxInt);
%% Apply CNN filter to find subdurals
image = imresize(tempImage, 1/8);
image = uint8(255 * mat2gray(image));
image = single(image)/255;
res = vl_simplenn(net, image);
subDural = (res(end).x*255) > 0;
subDural = imresize(subDural, 8);
subDural = finalizeSub(subDural, tempImage, distBinImage);
%% find brain and fluid
labelImage = totalSegStack(tempImage, distBinImage, subDural, thrs, thrsSeg);
FLISImage = labelImage;
