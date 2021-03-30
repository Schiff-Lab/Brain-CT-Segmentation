

function [centers, subLabel] = findSubLabel(data, cand, net)

%start_spams
tempImage = data;
distBinImage = bwdist(~cand);
image = imresize(data, 1/8);
image = uint8(255 * mat2gray(image));
image = single(image)/255;
res = vl_simplenn(net, image);
subDural = (res(end).x*255) > 0;
subDural = imresize(subDural, 8);
subDural = finalizeSub(subDural, data, distBinImage);

maxInt = max(tempImage(:));
R = tempImage(cand == 1);
thrsSeg = mean(R) + 3*std(R); % making sure only candidate region is being segmented
thrsSeg = min(thrsSeg, 2/3*maxInt);

X = distBinImage > 0;
cand = X & (~subDural) & (tempImage < thrsSeg);
%% use kmeans for segmenting
data1 = tempImage(cand == 1);
[~, centers] = kmeans(data1, 2);

centers = sort(centers);
thrs = .5*centers(1) + .5*centers(2);
fluid = (tempImage <= thrs) & (cand == 1);
brain = (tempImage > thrs) & (cand == 1);
finalImage = zeros(512,512);
finalImage(fluid == 1) = 2;
finalImage(brain == 1) = 1;
finalImage(subDural == 1) = 3;

meanB = mean(tempImage(finalImage == 1));
meanF = mean(tempImage(finalImage == 2));
meanS = mean(tempImage(finalImage == 3));
distB = abs(meanB - meanS);
distF = abs(meanF - meanS);

if(distB < distF)
    subLabel = 1;
else
    subLabel = 2;
end








