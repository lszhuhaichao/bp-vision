function [u,d,l,r] = bp_cb1(u,d,l,r,data,numIter)

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

tmp = dt(dst);

% truncate and store in destination vector
miniC = miniC + disc_k;

dst = min(tmp,miniC);

% normalize
val = sum(dst(:));
val = val/numLabels;
dst = dst - val;

end

%%%%%%% dt of 1d function
function d = dt(h)
    n = length(h);
    
    d = zeros(n,1);
    v = zeros(n,1);     % Locations of parabolas in lower envelope
    z = zeros(n+1,1);   % Locations of boundaries between parabolas
    
    v(1)=1;
    z(1) = -inf;
    z(2) = inf;
    j = 1;
    %% Compute lower envelope
    for q=2:n
        s = ((h(q)+(q-1)^2)-(h(v(j))+(v(j)-1)^2))/(2*(q-(v(j))));
        while(s<=z(j))
            j = j-1;
            s = ((h(q)+(q-1)^2)-(h(v(j))+(v(j)-1)^2))/(2*(q-(v(j))));
        end
        j = j + 1;
        v(j) = q;
        z(j) = s;
        z(j+1) = inf;
    end
    %% Fill in values of min convolution
    j=1;
    for q=1:n
        while(z(j+1)<q)
            j = j + 1;
        end
        d(q) = (q-(v(j)))^2+h(v(j));
    end
end