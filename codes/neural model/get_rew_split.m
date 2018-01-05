function [r,cor_sel]=get_rew_split(rew_dist,at,st)
val=rew_dist(at);
cor_sel=0;
if at==find(rew_dist==max(rew_dist),1)
    cor_sel=1;
end
if rand()<=val
    r=20;
else
    r=0;
end
end