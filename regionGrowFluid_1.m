function segmentedImage = regionGrowFluid_1(data, seg, cand, reg)

%[imwidth, imheight] = size(data);
segmentedImage = seg;

%meaSub = mean(data(seg == 3));

maxThrs = 2/3*(max(data(:)));
minThrs = 1/5*(min(data(:)));

tempReg = reg;
segRegVal = mode(seg(reg == 1));
iter = 0;
check = 1;
meanFluid = mean(data(seg == 2));
meanReg = mean(data(reg == 1));
    %meanReg = mean(data(reg == 1));
while(iter < 25 && check > 0)
    testReg = (reg == 0) & (cand == 1);
    [indx, indy] = find(testReg == 1);
    lenInd = length(indx);

    for ind = 1:lenInd
        X = reg(indx(ind)-1:indx(ind)+1,indy(ind)-1:indy(ind)+1);
        x = indx(ind);
        y = indy(ind);
        checkVal = sum(X(:));
        if(checkVal > 0)
            pixelVal = data(x,y);
            pixelSeg = seg(x,y);
            checkThrs1 = (pixelVal > minThrs) && (pixelVal < maxThrs);
            if(segRegVal == 3)
                if(pixelSeg == segRegVal && checkThrs1)
                    reg(x,y) = 1;
                end
            else
                checkVal1 = abs(meanReg - pixelVal);
                checkVal2 = abs(meanFluid - pixelVal);
                if(pixelSeg == segRegVal && checkVal1 <= 10 && checkThrs1 && checkVal2 <= 10 && checkThrs1)
                    reg(x,y) = 1;
                end
            end
        end
        
    end
    Diff = abs(tempReg - reg);
    tempReg = reg;
    check = sum(Diff(:));
    iter = iter + 1;
end
reg = imfill(reg, 'holes');
segmentedImage(reg == 1) = 2;