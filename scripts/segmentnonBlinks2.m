%% function segmentnonBlinks
%input trial data(5 column)
%output trial data(6 column)
function output=segmentnonBlinks2(trial)
         a=isnan(trial(:,5));
         b=size(trial,1);
         aaa=ones(b,1);
         t=1;
         for i = 2 : b
             if a(i)==a(i-1)
                 aaa(i)=t;
             else
                 t=t+1;
                 aaa(i)=t;
             end
         end
         output=[trial,aaa];    
                 
end