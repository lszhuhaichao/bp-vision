function [u,d,l,r] = bp_cb2(u,d,l,r,data,numIter)

%%% belief propagation using checkerboard update scheme
[height, width,~] = size(data);

for iter=1:numIter
    
    disp(['Iteration ' int2str(iter)]);
    for y=2:height-1
        for x=mod(y+iter,2)+2:2:width-1
            u(y,x,:) = msg(u(y+1,x,:), l(y,x+1,:), r(y,x-1,:), data(y,x,:));
            d(y,x,:) = msg(d(y-1,x,:), l(y,x+1,:), r(y,x-1,:), data(y,x,:));
            r(y,x,:) = msg(u(y+1,x,:), d(y-1,x,:), r(y,x-1,:), data(y,x,:));
            l(y,x,:) = msg(u(y+1,x,:), d(y-1,x,:), l(y,x+1,:), data(y,x,:));
        end
    end
    
end

end

%%% compute message
function dst = msg(s1,s2,s3,s4)

global disc_k;
global data_k;
global lambda;
global numLabels;

%%%% aggregate and find min
dst = s1 + s2 + s3 + s4;
miniC = min(dst(:));

dst = dt(dst);

% truncate and store in destination vector
miniC = miniC + disc_k;

dst(dst>miniC) = miniC;

% normalize
val = sum(dst(:));
val = val/numLabels;
dst = dst - val;

end

%%%%%%% dt of 1d function
function h = dt(h)
    for i=2:length(h)
        prev = h(i-1)+1;
        if(prev<h(i))
            h(i) = prev;
        end
    end
    
    for i=length(h)-1:-1:1
        prev = h(i+1)+1;
        if(prev<h(i))
            h(i) = prev;
        end
    end
end