function segmentedImage = regionGrowInclude(data, seg, cand, reg, meanVal)

[imwidth, imheight] = size(data);
temp = zeros(imwidth, imheight);
temp(2:imwidth-1, 2:imheight-1) = 1;
segmentedImage = seg;

%meaSub = mean(data(seg == 3));

maxThrs = 2/3*(max(data(:)));
minThrs = 1/5*(min(data(:)));

tempReg = reg;
iter = 0;
check = 1;

    meanReg = mean(data(reg == 1));
while(iter < 10 && check > 0)
    testReg = (reg == 0) & (temp == 1);
    [indx, indy] = find(testReg == 1);
    lenInd = length(indx);

    for ind = 1:lenInd
        X = reg(indx(ind)-1:indx(ind)+1,indy(ind)-1:indy(ind)+1);
        x = indx(ind);
        y = indy(ind);
        checkVal = sum(X(:));
        segVal = seg(x,y);
        if(checkVal > 0 && segVal ==0)
            pixelVal = data(x,y);
            checkThrs1 = (pixelVal > minThrs) && (pixelVal < maxThrs);
%             checkVal1 = abs(meanReg - pixelVal);
%             checkVal2 = abs(meaSub - pixelVal);
            if(checkThrs1)
                reg(x,y) = 1;
            end
        end
        
    end
    Diff = abs(tempReg - reg);
    tempReg = reg;
    check = sum(Diff(:));
    iter = iter + 1;
end

[xCoord, yCoord] = find(reg == 1);
if(isempty(xCoord))
   segmentedImage(reg == 1) = 0;
else
    for k = 1:length(xCoord)
        intVal = data(xCoord(k), yCoord(k));
        distC1 = abs(intVal - meanVal(1));
        distC2 = abs(intVal - meanVal(2));
        [~, lab] = min([distC2, distC1]);
        segmentedImage(xCoord(k), yCoord(k)) = lab;
    end
end