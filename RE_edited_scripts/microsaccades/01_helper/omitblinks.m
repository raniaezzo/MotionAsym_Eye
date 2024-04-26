%% function omitblinks
%input: data_in: origin dat_all.m  blink: corresponding blink data
%output: data_out: data with nan
function omitblinks(data_in,blink_path,omit,samplerate)

         omit_num=omit*samplerate/1000; %how mach time point before and after blink need to change NAN
         [path,filename, ~] = fileparts(data_in);
         load(data_in)
         load(blink_path)

         % how many blinks during trials for entire session
         blink_len=size(blink,1);

         time=Dat_all(:,1);
         for i = 1 : blink_len
             blink_start=blink(i,1);
             blink_end=blink(i,2);
             order=blink_start-omit_num<Dat_all(:,1)& blink_end-omit_num>Dat_all(:,1);
     
             time(order)=nan;
             
         end

         % additionally make nans any time outside of stimulus period that
         % are blinks (-1, -1, 0) - RE added 4-13-24
         for di=1:length(Dat_all)
             if Dat_all(di,2:4) == [-1, -1, 0]
                 time(di) = nan;
             end
         end


         Dat_all=[Dat_all,time];
         %save(sprintf('%s/%s_blink.mat',path,filename),'Dat_all');
         %data_out=sprintf('%s/%s_blink.mat',path,filename);
         save(sprintf('%s/%s.mat',path,filename),'Dat_all');
end
