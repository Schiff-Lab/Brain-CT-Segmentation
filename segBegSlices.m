

function FLISImage = segBegSlices(data, cand, thrs, segImage, subLabel)
%maxVal = max(data(:));
thrsB = (thrs(1) + thrs(2))/2;
if(subLabel == 2)
    segImage(data > thrsB & cand == 1) = 1;
    dist = bwdist(~cand);
    X = segImage == 3;
    Y = imfill(X,'holes');
    %% finding fluid
    [L, num] = bwlabel(Y);
    for i = 1:num
        binIm = L == i;
        edgeIm = edge(binIm);
        Dmin = min(dist(edgeIm));
        totLen = length(find(edgeIm));
        totLenB = length(find(binIm));
        thrsLen = length(find(dist(edgeIm) < Dmin + 6));
        ratio = thrsLen/totLen;
        if(ratio < .33 || Dmin > 5 || totLenB < 300)
            segImage(L == i) = 2;
        else
            segImage(L == i) = 3;
        end
    end
    X = segImage == 2;
    Y = bwareaopen(X, 50);
    Z = X&(~Y);
    segImage(Z == 1) = 1;
else
    segImage(data < thrsB & cand == 1) = 2;
    dist = bwdist(~cand);
    X = segImage == 3;
    Y = imfill(X,'holes');
    %% finding fluid
    [L, num] = bwlabel(Y);
    for i = 1:num
        binIm = L == i;
        edgeIm = edge(binIm);
        Dmin = min(dist(edgeIm));
        totLen = length(find(edgeIm));
        totLenB = length(find(binIm));
        thrsLen = length(find(dist(edgeIm) < Dmin + 6));
        ratio = thrsLen/totLen;
        if(ratio < .33 || Dmin > 5 || totLenB < 300)
            segImage(L == i) = 1;
        else
            segImage(L == i) = 3;
        end
    end
    X = segImage == 1;
    Y = bwareaopen(X, 50);
    Z = X&(~Y);
    segImage(Z == 1) = 2;
end

FLISImage = segImage;
