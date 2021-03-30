function finalBinary = findBinaryMiddleSlice_1(imageSlice, meanX)

thrsVal = .6*max(imageSlice(:));
B1 = imageSlice > 0;
B2 = imageSlice > thrsVal;
binaryImage = B1&(~B2);
[L, num] = bwlabel(binaryImage);
maxL = 0;
for i = 1:num
    lenL = length(find(L == i));
    if(lenL > maxL)
        maxL = lenL;
    end
end

for i = 1:num
    lenL = length(find(L == i));
    if(lenL < maxL)
        binaryImage(L ==i) = 0;
    end
end
finalBinary = imfill(binaryImage, 'holes');








