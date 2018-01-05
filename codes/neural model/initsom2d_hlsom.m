function wt = initsom2d_hlsom(m,n, X, mode)

% mode = 1: regular initialization from train data
% 		 = 2: random order initialization from train data
%		 = 3: initialization with random data

[N,dim] = size(X);
		% N = # of input vectors
		% dim = dimension of input vector

wt = zeros(m,n,dim);

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
   
   for i = 1:m,
      for j = 1:n,
          %ind = (randi(m-1)-1)*n+randi(n-1);
             ind=(i-1)*n+j;
             wt(i,j,:) = X(ind,:) + (2*rand(size(X(ind,:)))-1)/5;
      end
   end
end

% Random initialization 
if(mode == 3)   
   ind = ceil(N*rand(1,1));
   for i = 1:m,
      for j = 1:n,
         wt(i,j,:) = (rand(size(X(ind,:)))-0.0)/1;%;0.0+0.01*(rand(size(X(ind,:)))-0.0);
%          if rand()<0.5
%              wt(i,j,1)=0.25;
%          else
%              wt(i,j,1)=0.75;
%          end
      end
   end
end
