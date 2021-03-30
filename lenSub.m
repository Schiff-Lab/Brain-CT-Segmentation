

function l = lenSub(data, cand, net)

%start_spams
tempImage = data;
distBinImage = bwdist(~cand);
image = imresize(data, 1/8);
image = uint8(255 * mat2gray(image));
image = single(image)/255;
res = vl_simplenn(net, image);
subDural = (res(end).x*255) > 0;
subDural = imresize(subDural, 8);
subDural = finalizeSub(subDural, data, distBinImage);

thrs = 2/3*max(data(:));
X = subDural & (cand == 1) & (data < thrs);
l = length(find(X == 1));








