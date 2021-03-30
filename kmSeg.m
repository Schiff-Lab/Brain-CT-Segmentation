

function segImage = kmSeg(tempImage, testCand, centersKm)


segTemp = zeros(512,512);

%param.iter = 50;
%% segmentation only for test Image
maxInt = max(tempImage(:));
minInt = min(tempImage(:));
%R = tempImage(testCand == 1);
%thrsSeg = mean(R) + 3*std(R);
thrsSeg = 2/3*maxInt;%min(thrsSeg, 2/3*maxInt);
thrsKm = min(tempImage(testCand == 1));
thrsKm  = max(thrsKm, 1/5*minInt);

%% kMeans Clustering
cand = testCand;
[xCoord, yCoord] = find(cand);
if(isempty(xCoord))
    segTemp = 0;
else
    data = tempImage;
    for k = 1:length(xCoord)
        intVal = data(xCoord(k), yCoord(k));
        if(intVal > thrsSeg || intVal < thrsKm)
            continue;
        else
            distC1 = abs(intVal - centersKm(1));
            distC2 = abs(intVal - centersKm(2));
            [~, lab] = min([distC2, distC1]);
            segTemp(xCoord(k), yCoord(k)) = lab;
        end
    end
end
%% added newly
[L, num] = bwlabel(segTemp == 2);
for i = 1:num
    X = L == i;
    len = length(find(X == 1));
    if(len < 200)
        segTemp(L == i) = 1;
    end
end
   % segTemp = smallCorr(segTemp);
    segImage = segTemp;
    
    
    
    
    
    
    
    
    
