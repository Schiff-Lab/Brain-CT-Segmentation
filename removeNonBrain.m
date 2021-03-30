function brain = removeNonBrain(data, cand)

L = cand == 1;
[L, num] = bwlabel(L);
brain = cand;
maxL = 1;
%% removing small portions
for i = 1:num
    [xSub, ~] = find(L == i);
    if(length(xSub) > maxL)
        maxL = length(xSub);
    end
end
for i = 1:num
    [xSub, ~] = find(L == i);
    if(length(xSub) <= .1*maxL && length(xSub) < 1000)
        brain(L == i) = 0;
    end
end

[L, num] = bwlabel(brain);
for i = 1:num
    [x, y] = find(L == i);
    xMin = min(x);
    yMin = min(y);
    if(xMin < 128 && yMin < 192 && length(x)< .33*maxL)
        brain(L == i) = 0;
    end
end

[L, num] = bwlabel(brain);
for i = 1: num
    [x, y] = find(L == i);
    xMin = min(x);
    ratio = maxL/(length(x));
    if(ratio > 1 && xMin < 150 && length(x) < 2000)
        brain(L == i) = 0;
    end
end
