function [wt_hlsom,wt_subsom] = trainsom2d(wt_hlsom,wt_subsom, X, niter, bdrycond)

[N,dim] = size(X);
sig_hlsom = 0.01;
eta_hlsom = 0.6;
nhrad_hlsom = 1;
sig_subsom = 0.1;
eta_subsom = 0.4;
nhrad_subsom = 1;
n_subiter=500;
for i = 1:niter
    i
    X=X(randperm(size(X,1)),:);
    for j=1:N
        v=X(j,:)';
        indnn = nstnbrind(wt_hlsom, v);
        X1=gen_st_act_pairs(v');
        wt1=reshape(wt_subsom(indnn(1),indnn(2),:,:,:),size(wt_subsom,3),size(wt_subsom,4),size(wt_subsom,5));
        for k=1:n_subiter
            X1=X1(randperm(size(X1,1)),:);
            for l=1:size(X1,1)
                v1 = X1(l,:)';
                wt1 = adaptsom2d(wt1, v1, nhrad_subsom, sig_subsom, eta_subsom, bdrycond);   
            end
        end
        wt_subsom(indnn(1),indnn(2),:,:,:)=wt1;
        wt_hlsom = adaptsom2d(wt_hlsom, v, nhrad_hlsom, sig_hlsom, eta_hlsom, bdrycond);
    end
   % disp(['Iteration: ',num2str(i)])
    
end

