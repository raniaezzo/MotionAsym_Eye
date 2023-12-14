directiondegrees = [90,270,180,0,225,315,135,45];
subject = {'S01','S02','S03','S04','S05','S06'}; 
condition = {'Full_distance_non_radialtangential', 'Full_distance_radialtangential'}; 
direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
loction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
figpath='C:\Users\86186\Desktop\fig';
edges_plot=-0.125*pi:0.25*pi:1.875*pi;
for k = 1 :6
    for i = 1 : 8 
        name = [figpath,'/',subject{k},'_',direction{i},'_dir_count.png'];
        figure        
        title_name = [subject{k},'_ ',direction{i}];
        polarhistogram('BinEdges',edges_plot,'BinCounts',total_d_dir1(i,:,k),'FaceAlpha',.0,'LineWidth',1.5,'EdgeColor','r')
        hold on
        polarhistogram('BinEdges',edges_plot,'BinCounts',total_d_loc(i,:,k),'FaceAlpha',.0,'LineWidth',1.5,'EdgeColor','b')
%        legend('stimulus direction','stimulus location')
        title(title_name)
        saveas(gca,name)
    end
end



    for pp = 1 : 6
        title_name = [subject{pp}];
        figpath='C:\Users\86186\Desktop\fig';
        name = [figpath,'/','_',subject{pp},'_sumtotal.png'];
        figure
        polarhistogram('BinEdges',edges_plot,'BinCounts',loc(:,pp),'FaceAlpha',.0,'LineWidth',1.5)
        
        title(title_name)
        saveas(gca,name)
    end