%% this function to add outside num clone into tab file
function new_tab_outside=add_outside(tab_path,outside_path)
        [path,filename, ~] = fileparts(tab_path);
        load(tab_path);
        load(outside_path);
        len=size(tab_new,1);
        bl=zeros(len,1);
        for i = 1 : len
            order=tab_new(i,2)+1300 <=outside(:,2) & tab_new(i,2)+1800 >=outside(:,2);
            c1=outside(:,1);
            num=sum(c1(order));
            bl(i)=num;
        end
        tab=[tab_new,bl];
        save(sprintf('%s/%s_outside.mat',path,filename),'tab');
        new_tab_outside=sprintf('%s/%s_outside.mat',path,filename);
end