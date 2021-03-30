function segmentedImage = regionGrowBrain_1(data, seg, cand, reg)

%[imwidth, imheight] = size(data);
segmentedImage = seg;

%meaSub = mean(data(seg == 3));

maxThrs = 2/3*(max(data(:)));
minThrs = 1/5*(min(data(:)));

tempReg = reg;
segRegVal = mode(seg(reg == 1));
iter = 0;
check = 1;

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
            
            if(pixelSeg == segRegVal && checkThrs1)
                reg(x,y) = 1;
            end
        end
        
    end
    Diff = abs(tempReg - reg);
    tempReg = reg;
    check = sum(Diff(:));
    iter = iter + 1;
end

segmentedImage(reg == 1) = 1;