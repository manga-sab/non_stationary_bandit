
epi_learn_tot=[];
for l=1:1
    rng('default');
%Loads the congnitive and motor inputs as state and action inputs
load inps

%Setting the reward probabilities for the two actions
rp=0.05*(l-1);
rew_prob_4(:,1)=[1-rp;rp];
rew_prob_4(:,2)=[rp;1-rp];

%Setting the dimension of the Cognitive HLSOM and SUBSOM (Striosome and
%Matrisomes)
dim_cog_hlsom=[3,2];
dim_cog_subsom=[1,2];


%Number of trials in each session
no_trial=3000;

%Number of independent sessions
no_session=10;

%Number of episodes after which reward distribution is switched
switch_epi=500;

%Number of possible motor actions
no_act=2;

%Arbitrary initialisation of the probability vector for cognitive SOM for
%display
ps_cog=[1;0];

%Number of modules required in the task
no_modules=size(rew_prob_4,1);

%Track of rewards across trials and sessions
rews=[];

%Track of corrected selected option across trials and sessions
cor_sel=[];

%Track of Probability of each module across trials and sessions
pmods=[];


test_win_length=10;
epi_learn=[];


modus=[];
%Iterate across the number of sessions
for j=1:no_session
    
    epi_learn_1=[];
    sw_chk=0;
    j
    clc;
    
    %Train the cognitive striatum input representations
    disp('Training Cognitive SOMs');
    
    %Iterations to train the cognitive maps
    niter=10;
      wt_cog_hlsom = initsom2d_hlsom(dim_cog_hlsom(1),dim_cog_hlsom(2), st_inp, 2);
      wt_cog_subsom=initsom2d_subsom(dim_cog_hlsom(1),dim_cog_hlsom(2),dim_cog_subsom(1),dim_cog_subsom(2),gen_st_act_pairs(st_inp(1,:)), 2);
      [wt_cog_hlsom,wt_cog_subsom] = trainsom2d(wt_cog_hlsom,wt_cog_subsom, ac_inp,niter, 1);
%       
%       keyboard;
%    load ress
    clc;
    
    
    %Learning the task specific value functions
    disp('Learning');
    
    
    %Initialisation for the cognitive striatum
    %Initiasing the weights for the context identity signal
    wt_cog_p_calc_all=0.4*(1*rand(dim_cog_hlsom(1)*dim_cog_hlsom(2),no_modules)-0);
    
    %Initialising the weights for the state values
    wt_cog_v_calc_all=0.4*(1*rand(dim_cog_hlsom(1)*dim_cog_hlsom(2),no_modules)-0);
    
    %Initialising the weights for the action values
    wt_cog_q_calc_all=0.05*(1*rand(dim_cog_hlsom(1),dim_cog_hlsom(2),dim_cog_subsom(1)*dim_cog_subsom(2),no_modules));
    
    %Initialising the action values
    q_cog=zeros(no_act,1);
    
    %Setting the responsibility decay factor
    gam_mod_cog=0.95;
    
    %Initialising the responsibility signal
    lam_cog=zeros(no_modules,1);
    
    %Temperature factor for module selection
    beta_lam_cog=50;
    
    %Temperature Factor for selection of cognitive input
    beta_cog=10;
    
    %Learning rate of weights of the state value perceptron
    eta10_cog=0.05;
    
    %Learning rate of weights of the action value perceptron
    eta20_cog=0.0005;
    
    %Learning rate of weights of the context indicator perceptron
    eta_pt0_cog=0.05;
    eta_pt_cog=eta_pt0_cog;
    eta1_cog=eta10_cog;
    eta2_cog=eta20_cog;
    
    %Current context
    curr_mod_cog=1;
    
    %Storage and result variable initialisation
    mods_sel_cog=[];
    mods_sel_ratio_cog=[];
    epis_sel_cog=[];
    lam_disp_cog=zeros(no_modules,1);
    rew_dist=rew_prob_4(curr_mod_cog,:);
    rs=[];
    c_ss=[];
    
    pmod=[];
    
    sw_no=0;
    %Iterate over the trial
    for i=1:no_trial
        
        %Selecting two cognitive inputs
        cog_st=[1,0,0,1];
        
        
        
        disp(['Session no: ',num2str(j),' Trial no: ',num2str(i),' Responsibility cog: ',num2str(ps_cog(:)')])
        
        %Switch rewards after certain trials
        if mod(i,switch_epi)==0
            
            rew_dist=fliplr(rew_dist);
            sw_chk=0;
            sw_no=sw_no+1;
            if length(epi_learn_1)~=sw_no
                epi_learn_1=[epi_learn_1;switch_epi];
            end
        end
        mods_sel_cog=[];
        
        %Choosing a cognitive assumed context by using a softmax rule
        ps_cog=exp(beta_lam_cog*lam_cog/(max(abs(lam_cog))+1));
        ps_cog=ps_cog/sum(ps_cog);
        pmod=[pmod,ps_cog(:)'];
%         modu_cog=randsample_array(ps_cog',[1:no_modules]);
        modu_cog=find(ps_cog==max(ps_cog),1);
%          modu_cog=1;
%         modu_cog=modu_cog(randperm(length(modu_cog)));
%         modu_cog=modu_cog(1);
        
        mods_sel_ratio_cog=[mods_sel_ratio_cog;modu_cog];
        mods_sel_cog=[mods_sel_cog;modu_cog];
        
        %Extracting the required weights for the asummed context
        wt_cog_p_calc=wt_cog_p_calc_all(:,modu_cog);
        wt_cog_v_calc=wt_cog_v_calc_all(:,modu_cog);
        wt_cog_q_calc=reshape(wt_cog_q_calc_all(:,:,:,modu_cog),size(wt_cog_q_calc_all,1),size(wt_cog_q_calc_all,2),size(wt_cog_q_calc_all,3));
        
        %Generating all posibible input for cognitive action selection
        inps_cog=gen_st_act_pairs(cog_st);
        
        %Getting the state action values based on representation built
        indnn_cog = nstnbrind(wt_cog_hlsom, cog_st');
        y_hlsom_cog = respsom2d(cog_st,wt_cog_hlsom,0.01);
        vt_cog=y_hlsom_cog(:)'*wt_cog_v_calc;
        pt_cog=y_hlsom_cog(:)'*wt_cog_p_calc;
        wt_cog=reshape(wt_cog_subsom(indnn_cog(1),indnn_cog(2),:,:,:),size(wt_cog_subsom,3),size(wt_cog_subsom,4),size(wt_cog_subsom,5));
        for k=1:no_act
            y_subsom_cog = respsom2d(inps_cog(k,:),wt_cog,0.1);
            wt_cog_q_1=reshape(wt_cog_q_calc(indnn_cog(1),indnn_cog(2),:),dim_cog_subsom(1)*dim_cog_subsom(2),1);
            q_cog(k)=y_subsom_cog(:)'*wt_cog_q_1;
        end
        
        %Cognitive action selection based on softmax and Q values
        %calculated
        ptmp_cog=exp(q_cog*beta_cog);
        ptmp_cog=ptmp_cog/sum(ptmp_cog);
        
        at_cog=randsample_array(ptmp_cog(:)',[1:4]);
        
        %Get next state and reward
        [r,c_s]=get_rew_split(rew_dist,at_cog,cog_st);
        
        %Updating the weights using TD rule
        y_subsom_cog = respsom2d(inps_cog(at_cog,:),wt_cog,0.1);
        y_hlsom_cog_back=y_hlsom_cog;
        del_cog=r-vt_cog;
        del1_cog=r-pt_cog;
        
        
        lam_cog(modu_cog)=lam_cog(modu_cog)+0.05*(-1*lam_cog(modu_cog)-(0.1)*del1_cog.^2);
        
        
        %lam_disp_cog=lam_cog/sum(lam_cog);
        wvback_cog=wt_cog_v_calc;
        wt_cog_v_calc=update_wt_v(wt_cog_v_calc,y_hlsom_cog_back(:),del_cog,eta1_cog);
        
        
        wt_cog_v_calc(wt_cog_v_calc<0)=0;
        wt_cog_v_calc=20*wt_cog_v_calc/sum(wt_cog_v_calc);
        
        wt_cog_p_calc=update_wt_v(wt_cog_p_calc,y_hlsom_cog_back(:),del1_cog,eta_pt_cog);
        
        
        wt_cog_p_calc(wt_cog_p_calc<0)=0;
        wt_cog_p_calc=20*wt_cog_p_calc/sum(wt_cog_p_calc);
        
        
        wqback=wt_cog_q_1;
        wt_cog_q_1=update_wt_v(wt_cog_q_1,y_subsom_cog(:),del_cog,eta2_cog);
        
%         wt_cog_q_1(wt_cog_q_1<0)=0;
%         wt_cog_q_1=wt_cog_q_1/sum(wt_cog_q_1);
        
        wt_cog_q_calc(indnn_cog(1),indnn_cog(2),:)=wt_cog_q_1;
        wt_cog_v_calc_all(:,modu_cog)=wt_cog_v_calc;
        wt_cog_p_calc_all(:,modu_cog)=wt_cog_p_calc;
        wt_cog_q_calc_all(:,:,:,modu_cog)=wt_cog_q_calc;
        rs=[rs;r];
        c_ss=[c_ss;c_s];
        if i>(test_win_length+sw_no*switch_epi) && sw_chk==0
            if mean(c_ss(i-test_win_length+1:end))>=0.9
                epi_learn_1=[epi_learn_1;i-test_win_length-sw_no*switch_epi];
                sw_chk=1;
            end
        end
        %         eta_pt_cog=eta_pt0_cog-i*(eta_pt0_cog-1e-4)/no_trial;
        %         eta1_cog=eta10_cog-i*(eta10_cog-0.01)/no_trial;
        %         eta2_cog=eta20_cog-i*(eta20_cog-0.01)/no_trial;
        
    end
    modus=[modus;mods_sel_ratio_cog(:)'];
    rews=[rews;rs'];
    cor_sel=[cor_sel;c_ss'];
    pmods=[pmods;pmod];
    epi_learn=[epi_learn;epi_learn_1(:)'];
end


%
% %Plotting from the graphs

epi_learn_tot=[epi_learn_tot;mean(epi_learn)];
me=mean(cor_sel);
st=std(cor_sel);
% 
% 
figure(2*l+1);
shaded_plot(downsample([1:length(me)],1),downsample(me,1),downsample(st,1),'k');
% vline([switch_epi:switch_epi:no_trial],'r:','Switch');
saveas(2*l+1,['figure',num2str(2*l-1)],'png');
saveas(2*l+1,['figure',num2str(2*l-1)],'fig');
% % 
  save(['matte4_',num2str(l)]);
figure(2*l)
no_switch=switch_epi;
pc1=pmods(:,[1:2:2*no_trial]);
pc1(mean(pc1(:,no_switch:2*no_switch),2)<0.5,:)=1-pc1(mean(pc1(:,no_switch:2*no_switch),2)<0.5,:);
pc2=1-pc1;
subplot(2,1,1);
me=mean(pc1);
st=std(pc1);
shaded_plot(downsample([1:length(me)],1),downsample(me,1),downsample(st,1),'r');
hold on
ylim([-0.25,1.25]);
vline([no_switch:no_switch:no_trial],'k:','Switch');
hline(0.5,'k:');
xlabel('Trials');
ylabel('Probability');
title('Probability of context 1');
subplot(2,1,2);
me=mean(pc2);
st=std(pc2);
shaded_plot(downsample([1:length(me)],1),downsample(me,1),downsample(st,1),'r');
hold on
ylim([-0.25,1.25]);
vline([no_switch:no_switch:no_trial],'k:','Switch');
hline(0.5,'k:');
xlabel('Trials');
ylabel('Probability');
title('Probability of context 2');
saveas(2*l,['figure',num2str(2*l)],'png');
saveas(2*l,['figure',num2str(2*l)],'fig');
end
save('res5');