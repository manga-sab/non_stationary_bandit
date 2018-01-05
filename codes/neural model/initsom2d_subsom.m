function wt = initsom2d_subsom(a,b,m,n, X, mode)

% mode = 1: regular initialization from train data
% 		 = 2: random order initialization from train data
%		 = 3: initialization with random data

[N,dim] = size(X);
		% N = # of input vectors
		% dim = dimension of input vector

wt = zeros(a,b,m,n,dim);

% Regular initialization
if(mode == 1)
   M = m*n;
   diff = floor(N/M);
   ind = 1;
   for i = 1:m,
      for j = 1:n,
         wt(i,j,:) = X(ind,:);
         ind = ind + diff;
      end
   end
end

   
% Random order initialization from train data

if(mode == 2)   
   
   ind = ceil(N*rand(1,1));
    for i = 1:a,
        for j = 1:b,
            for k = 1:m,
                for l = 1:n,
                    wt(i,j,k,l,:) = X(ind,:) + rand(size(X(ind,:)))/100;%(rand(size(X(ind,:)))-0.0)/1;%;0.0+0.01*(rand(size(X(ind,:)))-0.0);
                end
            end
        end
    end
end

% Random initialization 
if(mode == 3)   
    ind = ceil(N*rand(1,1));
    for i = 1:a,
        for j = 1:b,
            for k = 1:m,
                for l = 1:n,
                    wt(i,j,k,l,:) = (rand(size(X(ind,:)))-0.0)/1;%;0.0+0.01*(rand(size(X(ind,:)))-0.0);
                end
            end
        end
    end
end
