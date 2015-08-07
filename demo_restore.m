clear all;
close all;

global disc_k;
global data_k;
global lambda;
global numLabels;

disc_k = 200;           % truncation of discontinuity cost
data_k = 10000;         % truncation of data cost
lambda = 0.05;          % weighting of data cost
numLabels = 256;        %

numIter = 5;        % number of BP iterations at each scale
numLevels = 5;      % number of scales


img = importdata('penguin.png');
img = img.cdata;
figure(1);
imshow(img);
img = imnoise(img,'gaussian',0,0.005);
figure(2);
imshow(img);

[height, width] = size(img);
img = double(img);
%%%%%%%%%%%%%%%%%%%%%%%%% compute data cost %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataCost = zeros(height, width, numLabels);

for i=1:numLabels
    gg = min((img-(i-1)).^2,data_k);
    dataCost(:,:,i) = lambda*gg;
end

%%%%%%%%%%%%%%%%%%%%%%%%% construct the data pyramid %%%%%%%%%%%%%%%%%%%%%

dataCostPyramid = cell(numLevels,1);
dataCostPyramid{1} = dataCost;
for j=2:numLevels
    [w,h] = size(imresize(dataCostPyramid{j-1}(:,:,1),0.5));
    dataCostPyramid{j} = zeros(w,h,numLabels);
    for i=1:numLabels
        d1 = downsample(downsample(dataCostPyramid{j-1}(:,:,i),2)',2)';
        d2 = downsample(downsample(dataCostPyramid{j-1}(:,:,i),2,1)',2,1)';
        d3 = downsample(downsample(dataCostPyramid{j-1}(:,:,i),2,0)',2,1)';
        d4 = downsample(downsample(dataCostPyramid{j-1}(:,:,i),2,1)',2,0)';
        dataCostPyramid{j}(:,:,i) = d1 + padarray(d2,size(d1)-size(d2),0,'post') + padarray(d3,size(d1)-size(d3),0,'post')+ padarray(d4,size(d1)-size(d4),0,'post');
    end
end
u = dataCostPyramid;
d = dataCostPyramid;
l = dataCostPyramid;
r = dataCostPyramid;
%%%%%%%%%%%%%%%%%%%%%%%%% run bp from coarse to fine %%%%%%%%%%%%%%%%%%%%%
for j=numLevels:-1:1
    if(j==numLevels)
         u{j} = u{j}*0;
         d{j} = d{j}*0;
         l{j} = l{j}*0;
         r{j} = r{j}*0;
    end
    if(j<numLevels)
        [h,w] = size(u{j}(:,:,i));
       for i=1:numLabels
            uu = imresize(u{j+1}(:,:,i),2,'nearest');
            dd = imresize(d{j+1}(:,:,i),2,'nearest');
            ll = imresize(l{j+1}(:,:,i),2,'nearest');
            rr = imresize(r{j+1}(:,:,i),2,'nearest');
           if(size(imresize(u{j+1}(:,:,i),2,'nearest'),1)>size(u{j}(:,:,i),1))
               uu(end,:) = [];
               dd(end,:) = [];
               ll(end,:) = [];
               rr(end,:) = [];
           end
           if(size(imresize(u{j+1}(:,:,i),2,'nearest'),2)>size(u{j}(:,:,i),2))
               uu(:,end) = [];
               dd(:,end) = [];
               ll(:,end) = [];
               rr(:,end) = [];
           end           
            u{j}(:,:,i) = uu;
            d{j}(:,:,i) = dd;
            l{j}(:,:,i) = ll;
            r{j}(:,:,i) = rr;          
        end 
    end
    [uu,dd,ll,rr] = bp_cb1(u{j},d{j},l{j},r{j},dataCostPyramid{j},numIter);
    u{j} = uu;
    d{j} = dd;
    l{j} = ll;
    r{j} = rr;
end

out = zeros(size(img));
%%%%%%%%%%%%%%%%%% generate output from current messages %%%%%%%%%%%%%%%%%%
for y = 2:height-1
    for x=2:width-1
        
        %%% Keep track of best balue for current pixel

        v = u{1}(y+1,x,:)+d{1}(y-1,x,:)+l{1}(y,x+1,:)+r{1}(y,x-1,:)+dataCostPyramid{1}(y,x);
        [vv,ii] = min(v);
        out(y,x) =  ii-1;
    end
end

out = uint8(out);
figure(3);
imshow(out);