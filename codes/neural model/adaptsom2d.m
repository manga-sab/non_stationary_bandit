function	wt = adaptsom2d(wt, v, nhrad, sig, eta, bdrycond)

[m,n,dim] = size(wt);

indnn = nstnbrind(wt, v);	%index of the nearest neighbor
%keyboard
[nbhood, nndist] = neighborhood2(m,n, indnn, nhrad, bdrycond);

[mx,nx] = size(nbhood);

for i = 1:mx,
   ind = nbhood(i,:);
   wv = reshape(wt(ind(1),ind(2), :), dim,1);
   diff = (v-wv);
%   diffi = (indnn - ind);
   dis = nndist(i)/(sig*sig);
   delta = eta.*(exp(-dis)*diff);
   %delta = eta.*(diff);
%	 delta = 2*eta*(exp(-dis)-exp(-dis/2))*diff;
	  
     wt(ind(1),ind(2), :) = wv + delta;
 %  keyboard
end
