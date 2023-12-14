%% draw figure (time from 0-2500)
function MS_new=drawfigure(MS,samplerate,tab,filepath)
locationids = 1:8; locationdegrees = {315,135,225,45,270,90,180,0};
locationlabels = strcat('loc',cellfun(@num2str,locationdegrees,'un',0));
[path,filename, ~] = fileparts(filepath);

time = samplerate/1000;

fil=MS(:,9)/time < 2500;
MS_new=MS(fil,:);
MS=MS(fil,:);
con1=tab(MS(1,10),9);
a=figure
hold on 
plot([0,0],[0,0],'b')
hold on
plot([0,0],[0,0],'r')
hold on 
plot([0,0],[0,0],'g')
hold on
plot([0,0],[0,0],'y')
% if con1 == 1 ||  con1 == 2   || con1 == 3 ||  con1 == 4   
%     legend(locationlabels{1},locationlabels{2},locationlabels{3},locationlabels{4})
% else
%     legend(locationlabels{5},locationlabels{6},locationlabels{7},locationlabels{8})
% end
for i = 1 : size(MS,1)
    con=tab(MS(i,10),9);
    hold on
    if con == 1 ||  con == 5       
        plot([MS(i,8)/time,MS(i,9)/time],[MS(i,10),MS(i,10)],'b')
%         if con == 1             
%             legend(locationlabels{1})
%         else
%             legend(locationlabels{5})
%         end
        
    elseif con == 2 ||  con == 6
        plot([MS(i,8)/time,MS(i,9)/time],[MS(i,10),MS(i,10)],'r')
%         if con == 2             
%             legend(locationlabels{2})
%         else
%             legend(locationlabels{6})
%         end
        
    elseif con == 3  ||  con == 7 
        plot([MS(i,8)/time,MS(i,9)/time],[MS(i,10),MS(i,10)],'g')
%         if con == 3             
%             legend(locationlabels{3})
%         else
%             legend(locationlabels{7})
%         end
        
    else
        plot([MS(i,8)/time,MS(i,9)/time],[MS(i,10),MS(i,10)],'y')
%         if con == 4             
%             legend(locationlabels{4})
%         else
%             legend(locationlabels{8})
%         end
        
    end
    hold on
    plot([1300,1300],[0,800],'k')
    hold on
    plot([1800,1800],[0,800],'k')
if con1 == 1 ||  con1 == 2   || con1 == 3 ||  con1 == 4   
    legend(locationlabels{1},locationlabels{2},locationlabels{3},locationlabels{4})
else
    legend(locationlabels{5},locationlabels{6},locationlabels{7},locationlabels{8})
end
%     hold on 
%     plot([0,0],[0,0],'b')
%     hold on
%     plot([0,0],[0,0],'r')
%     hold on 
%     plot([0,0],[0,0],'g')
%     hold on
%     plot([0,0],[0,0],'y')
%     if con == 1 ||  con == 2   || con == 3 ||  con == 4   
%         legend(locationlabels{1},locationlabels{2},locationlabels{3},locationlabels{4})
%     else
%         legend(locationlabels{5},locationlabels{6},locationlabels{7},locationlabels{8})
%     end
end

savefig(a,sprintf('%s/%s_figure1.fig',path,filename))
end