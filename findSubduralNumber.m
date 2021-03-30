function subNo = findSubduralNumber(data, cand, net, netSeg)

sliceNumber = length(data(1,1,:));
middleNumber = round(sliceNumber/2);
predSub = zeros(1,5);
k = 1;
for i = middleNumber:middleNumber + 4
    image = data(:,:,i);
    image = imresize(image, 1/8);
    image = uint8(255 * mat2gray(image));
    image = single(image)/255;
    image = 256 * (image - net.imageMean) ;
    res = vl_simplenn_class(net, image) ;
    [score, predSub(k)] = max(squeeze(res(end).x(1,1,:))) ;
    k = k+1;
end
l = zeros(1,5);
k = 1;
lenScore = length(find(predSub == 1));
if(lenScore < 3)
    subNo = sliceNumber + 2;
else
    for i = middleNumber:middleNumber + 4
        X = data(:,:,i);
        C = cand(:,:,i);
        l(k) = lenSub(X, C, netSeg);
        k = k+1;
    end
    lengthL = length(find(l < 10000));
    if(lengthL > 1)
        subNo = sliceNumber + 2;
    else
        for i = 1:sliceNumber
            image = data(:,:,i);
            image = imresize(image, 1/8);
            image = uint8(255 * mat2gray(image));
            image = single(image)/255;
            image = 256 * (image - net.imageMean) ;
            res = vl_simplenn_class(net, image) ;
            [score, pred] = max(squeeze(res(end).x(1,1,:)));
            if(pred == 1)
                subNo = i;
                break;
            end
        end
    end
end