function finalBinary = crossSliceBinaryNew_Pre(imageSliceCurr, binaryImageSliceOld, meanX)
%% initial classification
check = imageSliceCurr > .55*max(imageSliceCurr(:)) ;
test = (~check) & (imageSliceCurr ~= 0);
test1 = bwareaopen(test, 500);
len = length(find(test1));
if(len < 7500)
    finalBinary = zeros(512,512);
else
    Z = ~test;
    B2 = bwareaopen(Z, 20);
    B1 = ~B2;
    se = strel('disk', 20);
    fillCurrent = imopen(B1, se);
    finalBinary = fillCurrent & binaryImageSliceOld;
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
        if(length(x) < maxVal)
            finalBinary(L == i) = 0;
        end
    end
end
finalBinary = imfill(finalBinary, 'holes');
% se = strel('disk', 30);
% finalBinary = imopen(finalBinary, se);
