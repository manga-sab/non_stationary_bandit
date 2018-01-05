function ind = nstnbrind(wt,v)

[m,n,dim] = size(wt);
wtvec = reshape(wt(1,1,:), dim,1);
dis = (wtvec - v)'*(wtvec - v);
ind = [1 1];

for i = 1:m,
   for j = 1:n,
      wv = reshape(wt(i,j,:), dim,1);
      tmp =  (wv- v)'*(wv - v);
      if(tmp < dis)
         ind = [i j];
         dis = tmp;
      end
   end
end

      

