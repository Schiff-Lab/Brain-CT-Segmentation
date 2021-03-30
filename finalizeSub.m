function finalImage = finalizeSub(sub, image, dist)

%% eliminate which are not in cand region
X = dist > 0;
sub = sub & X;
%% next eliminate which are far from boundary
[L, num] = bwlabel(sub);
for i = 1:num
    X = L == i;
    minD = min(dist(X));
    if(minD > 10)
        sub(L == i) = 0;
    end
end
finalImage = sub;
