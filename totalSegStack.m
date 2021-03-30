function finalImage = totalSegStack(image, dist, sub, centers, thrShunt)

%% find cand reg for segmenting into brain and CSF
X = dist > 0;
cand = X & (~sub) & (image < thrShunt);
finalImage = zeros(512,512);
%% use kmeans for segmenting
data = image(cand == 1);
[~, centers1] = kmeans(data, 2);
centers1 = sort(centers1);
thrs = .5*centers(1) + .5*centers(2);
fluid = (image <= thrs) & (cand == 1);
brain = (image > thrs) & (cand == 1);
finalImage(fluid == 1) = 2;
finalImage(brain == 1) = 1;
finalImage(sub == 1) = 3;
%% do final post processing
meanB = mean(image(finalImage == 1));
meanF = mean(image(finalImage == 2));
meanS = mean(image(finalImage == 3));
distB = abs(meanB - meanS);
distF = abs(meanF - meanS);
if(distB < distF)
    notSub = (finalImage == 3) & (image <= centers(1));
    finalImage(notSub == 1) = 2;
    [L, num] = bwlabel(finalImage == 2);
    for i = 1:num
        X = L == i;
        len = length(find(X == 1));
        if(len < 200)
            finalImage(L == i) = 1;
        end
    end
end
notSub = (finalImage == 3) & (image >= centers1(2)) & (dist >= 7.5);
finalImage(notSub == 1) = 1;

%% remove small subs 
[L, num] = bwlabel(finalImage == 3);
for i = 1:num
    X = L == i;
    len = length(find(X == 1));
    if(len < 200)
        m = mean(image(X == 1));
        if(m < thrs)
          finalImage(L == i) = 2;
        else
           finalImage(L == i) = 1; 
        end
    end
end

[L, num] = bwlabel(finalImage == 2);
for i = 1:num
    X = L == i;
    len = length(find(X == 1));
    if(len < 200)
        finalImage(L == i) = 1;
    end
end
%% remove shunt

finalImage(image >= thrShunt) = 0;
%% remove subs from large distances
[L, num] = bwlabel(finalImage == 3);
for i = 1:num
    X = L == i;
    minD = min(dist(X));
    if(minD > 10)
        m = mean(image(X == 1));
        if(m < thrs)
           finalImage(L == i) = 2;
        else
           finalImage(L == i) = 1; 
        end
    end
end
%% close small brains inside sub
sub = bwlabel(finalImage == 3);
negSub = ~sub;
negSub = bwareaopen(negSub, 200);
sub = ~negSub;
finalImage(sub == 1) = 3;


