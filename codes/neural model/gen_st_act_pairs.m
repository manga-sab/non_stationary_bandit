function x=gen_st_act_pairs(st)
x=[];
no_act=2;
for i=1:size(st,1)
x=[x;[eye(no_act)]];
end
end