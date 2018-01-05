function y=respsom2d(x,wt,sig)

[m,n,dim] = size(wt);
y = zeros(m,n);
if(dim ~= length(x))
   fprintf('Invalid input size in respsom2d()\n')
   return
end

[mm, nn] = size(x);
if(nn == dim)
   x = x';
end
%sig = 0.2;
for i = 1:m,
   for j = 1:n,
      v = reshape(wt(i,j,:), dim,1);
      dif = v-x;
      dis = (sum(dif.*dif));
      y(i,j) = exp(-dis/(sig*sig));
   end
end

