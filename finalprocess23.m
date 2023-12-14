% RESPONSE PERIOD

clc; clear all;
homedir = '/Users/rania/Downloads/MS_Project';
tab_total = [] ;
MS_total = [] ; 
sample=[];
sub={'S01','S02','S03','S04','S05','S06','S07','S08'};
subject = 'S01'; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
direction_name = {'VU_new','VL_new','HL_new','HR_new','LL_new','LR_new','UL_new','UR_new'};
ACC_W=[];
ACC_NW=[];
REC_W=[];
REC_NW=[];
Y_car=zeros(10,2);
Tab_w_cl=[];
Tab_nw_cl=[];
Tab_w_ncl=[];
Tab_nw_ncl=[];

acc_total_w=[];
acc_total_nw=[];
rec_total_w=[];
rec_total_nw=[];
d_total_w=[];
d_total_nw=[];
figpath='C:\Users\86186\Desktop\fig\new\new_pre';
color = {[127, 191, 123]	,[175, 141, 195]}	;

% color={[215,25,28],[253,174,97],[255,255,191],[171,217,233],[44,123,182]};	
load(fullfile(homedir,'reaction.mat'))

axe_y_car  = zeros([8,5]);
axe_y_obl  = zeros([8,5]);

for k = 1 : 8
    
    if k == 5         %|| kk == 7
        continue
    end
    
    Y_car=zeros(10,2);
    Y_obl=zeros(10,2);
    y_total=0;
    y_t = 0 ; 
    for i = 1 : 2
        if i ==1 & k > 6
            continue
        end
        for j = 1 : 8

            main_folder = fullfile(homedir,'Data_DI_wEYE', sub{k}, ...
                 'RawData', condition{i}, 'Block1');
            cd(fullfile(main_folder, 'eyedata'));
            edf_name = dir(sprintf('*%s*.edf', direction{j})).name;
            edf_path = fullfile(main_folder,'eyedata',edf_name);
            msg_filepath=replace(edf_path,'edf','msg');
            samplingRateData=findSamplingRate(msg_filepath);
            sample=[sample,samplingRateData];
            MATpath = fullfile(main_folder, 'eyedata','MATs');
    %         cd(fullfile(main_folder, 'eyedata','MATs'));
    %         ms_name = dir(sprintf('%s.mat', direction{j})).name;
            ms_path= fullfile(MATpath,sprintf('%s.mat', direction_name{j}) );
            load(ms_path);
    %         MS_TEMP(:,10)=MS_TEMP(:,10)+ ((i-1)*8 + (j-1))*800;   
    %         MS_total=[MS_total;MS_TEMP];
    %         MATpath = fullfile(main_folder, 'eyedata','MATs');
            tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
            load(tab_path)
    %         tab(:,1)=tab(:,1)+ ((i-1)*8 + (j-1))*800;  
    %         tab_total=[tab_total;tab];
    
            %mean_RT = mean(tab(:,14)); % if session specific RT
            
            load(fullfile(homedir,'reaction_median.mat'))
            
            mean_RT = react(k);
            
            [acc_w,acc_nw,rec_w,rec_nw,y,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,MS_TEMP,samplingRateData,[0,300]); %react(k)]);
            if j < 5 
                
                Y_car=Y_car+y;
            else
                Y_obl=Y_obl+y;
            end

        end


    end

    y_total=y_total+Y_car;
    y_t =y_t+Y_obl;
    
%     figure
    X = categorical({'0.5','1','2','4','8'});
    X = reordercats(X,{'0.5','1','2','4','8'});
    
    y_ax = zeros([5,1]);
    y_ax1 = zeros([5,1]);
    for or = 1 : 5 
        
        y_ax(or)  = y_total(or,1) + y_total(11-or,1);

        y_ax1(or) = y_t(or,1) + y_t(11-or,1);
        
    end
    axe_y_car(k,:) = y_ax;
    axe_y_obl(k,:) = y_ax1;
    
    
    
%     y_use = [y_ax,y_ax1]';
%     bs = bar(X,y_use','stacked');
% 
%     for sd = 1 : 2 
%         bs(1,sd).FaceColor = 'flat';
%         bs(1,sd).CData = color{sd}/255;
%     end
% %     
% %     
%     title_na = ['Example Subject ',sub{k}];
%     title(title_na)
%     xlabel('Tilt from standard (degrees)')
%     ylabel('% trials with microsaccade during stimulus')
%     ylim([0,ceil(max(y_ax))+2])
%     for s = 1 : 10
%         str1 = [num2str(y_ax(s)),'%'];
%         text(s,y_ax(s),str1,...   
%                'HorizontalAlignment','center',...    
%                'VerticalAlignment','bottom')
%     end  
%     name = [figpath,'/',sub{k},'-ase.png'];
%     title(name)
%    saveas(gca,name)
    
end

yy_car = zeros([1,40]);
yy_obl = zeros([1,40]);
nSubs = 0;

for ss = 1 : 8 
    if ss == 5
        continue
    end
    yy_car((ss-1)*5+1:(ss*5)) = axe_y_car(ss,:);
    yy_obl((ss-1)*5+1:(ss*5)) = axe_y_obl(ss,:);
    nSubs = nSubs+1;
end

yy_car = [yy_car(1:20), yy_car(26:end)];
yy_obl = [yy_obl(1:20), yy_obl(26:end)];

% re-order
card_order = reshape(yy_car, [5,nSubs]);
obl_order = reshape(yy_obl, [5,nSubs]);
total_order = [card_order; obl_order];
sumIdx = [7 1 2 6 4 3 5]; % based on stimPeriod order

newCard = reshape(card_order(:,sumIdx), [1,5*nSubs]);
newObl = reshape(obl_order(:,sumIdx), [1,5*nSubs]);

yy_car = newCard;
yy_obl = newObl;

figure
%xt = [0.7,0.85,1,1.15,1.3,1.7,1.85,2,2.15,2.3,2.7,2.85,3,3.15,3.3,3.7,3.85,4,4.15,4.3,4.7,4.85,5,5.15,5.3,5.7,5.85,6,6.15,6.3,6.7,6.85,7,7.15,7.3,7.7,7.85,8,8.15,8.3];
xt = [0.7,0.85,1,1.15,1.3,1.7,1.85,2,2.15,2.3,2.7,2.85,3,3.15,3.3,3.7,3.85,4,4.15,4.3,4.7,4.85,5,5.15,5.3,5.7,5.85,6,6.15,6.3,6.7,6.85,7,7.15,7.3]; %,7.7,7.85,8,8.15,8.3];
bs = bar(xt,[yy_car;yy_obl],'stacked');
for sd = 1 : 2 
    bs(1,sd).FaceColor = 'flat';
    bs(1,sd).CData = color{sd}/255;
end    
xtic = [1,2,3,4,5,6,7,8];
set(gca,'XTick',xtic); 
set(gca,'XTicklabels',sub); 
%title('Response period')
%xlabel('Subject')
ylabel('# of MS')

ax = gca;
ax.FontSize = 20;
ax.LineWidth = 2;
fig1 = gcf;
fig1.Position = [737 520 686 211];
ax.XTick = [];
box off
ylim([0 800])

% figure
% xt = [0.7,0.85,1,1.15,1.3,1.7,1.85,2,2.15,2.3,2.7,2.85,3,3.15,3.3,3.7,3.85,4,4.15,4.3,4.7,4.85,5,5.15,5.3,5.7,5.85,6,6.15,6.3,6.7,6.85,7,7.15,7.3,7.7,7.85,8,8.15,8.3];
% bs = bar(xt,[yy_car;yy_obl],'stacked');
% for sd = 1 : 2 
%     bs(1,sd).FaceColor = 'flat';
%     bs(1,sd).CData = color{sd}/255;
% end    
% xtic = [1,2,3,4,5,6,7,8];
% set(gca,'XTick',xtic); 
% set(gca,'XTicklabels',sub); 
% title('Reaction period')
% xlabel('Subject')
% ylabel('# of MS')
% 
% name='C:\Users\86186\Desktop\fig\new\new_pre\ms_count_re.png';
% saveas(gca,name)