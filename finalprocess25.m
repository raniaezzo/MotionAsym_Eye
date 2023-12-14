clc; clear all;
homedir = '/Users/rania/Downloads/MS_Project';
subject = {'S01','S02','S03','S04','S05','S06','S07','S08'}; condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
location = {'UR','LL','UL','LR','VU','VL','HR','HL'};
direction_name = {'VU_new','VL_new','HL_new','HR_new','LL_new','LR_new','UL_new','UR_new'};
oblique_cum = zeros(5,8);
cardinal_cum = zeros(5,8);
tilt_angle = [0.5,1,2,4,8];
time = [0,0];
color = {[127, 191, 123]	,[175, 141, 195]}	;
figure
for kk = 1 : 8 
    tab_nw_total = [];
    tab_w_total = [];

    c05_a = [];
    c1_a = [];
    c2_a = [];
    c4_a = [];
    c8_a = [];
    o05_a = [];
    o1_a = [];
    o2_a = [];
    o4_a = [];
    o8_a = [];
    
    c05_a_n = [];
    c1_a_n = [];
    c2_a_n = [];
    c4_a_n = [];
    c8_a_n = [];
    o05_a_n = [];
    o1_a_n = [];
    o2_a_n = [];
    o4_a_n = [];
    o8_a_n = [];
    
    
    
    for i = 1 : 2
        if kk > 6 & i == 1 
            continue
        end
        for j = 1 : 8
            if kk == 3 & i == 1 & j==1 
                continue
            end
            main_folder = fullfile(homedir,'Data_DI_wEYE', subject{kk}, ...
                 'RawData', condition{i}, 'Block1');
            cd(fullfile(main_folder, 'eyedata'));
            edf_name = dir(sprintf('*%s*.edf', direction{j})).name;
            edf_path = fullfile(main_folder,'eyedata',edf_name);
            msg_filepath=replace(edf_path,'edf','msg');
            samplingRateData=findSamplingRate(msg_filepath);
%            sample=[sample,samplingRateData];
            MATpath = fullfile(main_folder, 'eyedata','MATs');
            ms_path= fullfile(MATpath,sprintf('%s.mat', direction_name{j}) );
            load(ms_path);
            tab_path = fullfile(MATpath, replace(edf_name, '.edf', '_tab_new_outside_blink.mat'));
            load(tab_path)

            [acc_w,acc_nw,rec_w,rec_nw,y,tab_w_cl,tab_nw_cl,tab_w_ncl,tab_nw_ncl]=acc_rec_count(tab,MS_TEMP,samplingRateData,time);
            
            tab_nw = [tab_nw_cl;tab_nw_ncl];
            tab_w = [tab_w_cl;tab_w_ncl];
            tab_nw_total = [tab_nw_total;tab_nw];
            tab_w_total = [tab_w_total;tab_w];
            
            if j < 5 
                for ll = 1 : size(tab_nw,1)
                    switch tab_nw(ll,11)
                        case 0.5 
      
                           c05_a = [c05_a;tab_nw(ll,15)];
                        case 1 
                         
                           c1_a = [c1_a;tab_nw(ll,15)];
                        case 2 
                    
                           c2_a = [c2_a;tab_nw(ll,15)];
                        case 4 
                        
                           c4_a = [c4_a;tab_nw(ll,15)];
                        case 8 
                         
                           c8_a = [c8_a;tab_nw(ll,15)];
                    end
                end
                
                for lls = 1 : size(tab_w,1)
                    switch tab_w(lls,11)
                        case 0.5 
      
                           c05_a_n = [c05_a_n;tab_w(lls,15)];
                        case 1 
                         
                           c1_a_n = [c1_a_n;tab_w(lls,15)];
                        case 2 
                    
                           c2_a_n = [c2_a_n;tab_w(lls,15)];
                        case 4 
                        
                           c4_a_n = [c4_a_n;tab_w(lls,15)];
                        case 8 
                         
                           c8_a_n = [c8_a_n;tab_w(lls,15)];
                    end
                end
                
                
                
            else
                for ll = 1 : size(tab_nw,1)
                    switch tab_nw(ll,11)
                        case 0.5 
                           
                           o05_a = [o05_a;tab_nw(ll,15)];
                        case 1 
                          
                           o1_a = [o1_a;tab_nw(ll,15)];
                        case 2 
                          
                           o2_a = [o2_a;tab_nw(ll,15)];
                        case 4 
                          
                           o4_a = [o4_a;tab_nw(ll,15)];
                        case 8 
                          
                           o8_a = [o8_a;tab_nw(ll,15)];
                    end
                end
                
                for lls = 1 : size(tab_w,1)
                    switch tab_w(lls,11)
                        case 0.5 
      
                           o05_a_n = [o05_a_n;tab_w(lls,15)];
                        case 1 
                         
                           o1_a_n = [o1_a_n;tab_w(lls,15)];
                        case 2 
                    
                           o2_a_n = [o2_a_n;tab_w(lls,15)];
                        case 4 
                        
                           o4_a_n = [o4_a_n;tab_w(lls,15)];
                        case 8 
                         
                           o8_a_n = [o8_a_n;tab_w(lls,15)];
                    end
                end
                
                
            end
            
            
            
            
        end
    end
    

    card_ob_a = zeros([5,2]);
    card_ob_a_n = zeros([5,2]);

    card_ob_a(1,1) = mean(c05_a);
    card_ob_a_n(1,1) = mean(c05_a_n);
    
    card_ob_a(2,1) = mean(c1_a);
    card_ob_a_n(2,1) = mean(c1_a_n);
    
    card_ob_a(3,1) = mean(c2_a);
    card_ob_a_n(3,1) = mean(c2_a_n);
    
    card_ob_a(4,1) = mean(c4_a);
    card_ob_a_n(4,1) = mean(c4_a_n);
 
    card_ob_a(5,1) = mean(c8_a);
    card_ob_a_n(5,1) = mean(c8_a_n);
    

    card_ob_a(1,2) = mean(o05_a);
    card_ob_a_n(1,2) = mean(o05_a_n);

    card_ob_a(2,2) = mean(o1_a);
    card_ob_a_n(2,2) = mean(o1_a_n);

    card_ob_a(3,2) = mean(o2_a);
    card_ob_a_n(3,2) = mean(o2_a_n);

    card_ob_a(4,2) = mean(o4_a);
    card_ob_a_n(4,2) = mean(o4_a_n);

    card_ob_a(5,2) = mean(o8_a);
    card_ob_a_n(5,2) = mean(o8_a_n);
    
    y_c = 0.9:4.9;
    y_o = 1.1:5.1;
    scatter(card_ob_a(:,1)',y_c,[],color{1}/255)
    hold on     
    scatter(card_ob_a_n(:,1)',y_c,[],color{1}/255,'filled')
    hold on 
    scatter(card_ob_a(:,2)',y_o,[],color{2}/255)
    hold on 
    scatter(card_ob_a_n(:,2)',y_o,[],color{2}/255,'filled')
    hold on
%     if sum(card_ob_r(:,1) > 3) > 0 
%         print(kk)
%     end

%     if sum(card_ob_r<0.2) == 0 
%         kk
%     end
    
end

yticks(1:5)
yticklabels({'0.5','1','2','4','8'})
ylabel('tilt angle')
xlabel('Accuracy (%)')
name = 'C:\Users\86186\Desktop\fig\new\accur2.png';
legend({'Cardinal not with MS','Cardinal with MS','Oblique not with MS','Oblique with MS'},'Location','northeastoutside')
%saveas(gca,name) 