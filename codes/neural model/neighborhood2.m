function [nbrind, nndist]= neighborhood2(m,n, ind, nhrad, bdrycond)

%Returns the indices of the nodes in the neighborhood
%bdrycond = 1: not periodic
%		= 2: periodic

if(bdrycond == 1)
   % Non-periodic
   nh2 = (2*nhrad+1)*(2*nhrad+1);
   nbrind = zeros(nh2,2);
    for i = 1:2*nhrad+1,
   	for j = 1:2*nhrad+1,
      	ix = (i-1)*(2*nhrad+1)+j;
      	nbrind(ix,:) = [i-nhrad-1, j-nhrad-1] + ind;
   	end
   end
   ii = find((nbrind(:,1)< 1)|(nbrind(:,2)< 1));
   if(length(ii)>0)
      nbrind(ii,:)=[];
   end
   ii = find((nbrind(:,1)> m)|(nbrind(:,2)> n));
   if(length(ii)>0)
      nbrind(ii,:)=[];
   end
   [mm,nn] = size(nbrind);
   nndist = zeros(mm,1);
   for i = 1:mm,
      diff = nbrind(i,:) - ind;
      nndist(i) = sum(diff.*diff);
   end
end

if(bdrycond == 2)
   % Periodic
   nh2 = (2*nhrad+1)*(2*nhrad+1);
   nbrind = zeros(nh2,2);
	for i = 1:2*nhrad+1,
   	for j = 1:2*nhrad+1,
      	ix = (i-1)*(2*nhrad+1)+j;
      	nbrind(ix,:) = [i-nhrad-1, j-nhrad-1] + ind;
   	end
   end
   [mi,ni] = size(nbrind);
   nndist = zeros(mi,1);
   for i = 1:mi,
      diff = nbrind(i,:) - ind;
      nndist(i) = sum(diff.*diff);
   end

   for i = 1:mi,
      if(nbrind(i,1) < 1)
         nbrind(i,1) = m + nbrind(i,1);
      end
      if(nbrind(i,1) > m)
         nbrind(i,1) = nbrind(i,1) - m;
      end
      if(nbrind(i,2) < 1)
         nbrind(i,2) = n + nbrind(i,2);
      end
      if(nbrind(i,2) > n)
         nbrind(i,2) = nbrind(i,2) - n;
      end
   end
end