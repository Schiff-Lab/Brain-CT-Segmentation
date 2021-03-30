function finalBinary = crossSliceBinary_Pre(imageSliceCurr, binaryImageSliceOld)

meanX = mean(imageSliceCurr(imageSliceCurr > 0));
maxVal = max(imageSliceCurr(:));
thrs = max(meanX, .5*maxVal);
thrs = min(thrs, .56*maxVal);
%stdX = std(imageSliceCurr(imageSliceCurr > 0));
B1 = imageSliceCurr > thrs;
B1 = (~B1)&(imageSliceCurr ~= 0);
B2 = ~B1;
B2 = bwareaopen(B2, 20);
B1 = ~B2;
se = strel('disk', 10);
fillCurrent = imopen(B1, se);
finalBinary = imfill(fillCurrent, 'holes');
%% to remove small components
[L, num] = bwlabel(finalBinary);
maxVal = 0;
for i = 1:num
    F = find(L==i);
    if(length(F) > maxVal)
        maxVal = length(F);
    end
end

for i = 1: num
    [x, y] = find(L == i);
    ratio = maxVal/(length(x));
    X = L == i;
    Y = X & binaryImageSliceOld;
    ratioCross = length(find(Y == 1))/length(find(X == 1));
    maxX = max(x);
    if(maxX > 255)
        thrsRatio = .60;
    else
        thrsRatio = .75;
    end
    %ratio > 20 || 
    if(ratioCross < thrsRatio)
        finalBinary(L == i) = 0;
    end
end

