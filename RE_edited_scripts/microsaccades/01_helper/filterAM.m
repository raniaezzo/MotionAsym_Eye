%% this script is used to filter AM
%input original MS
%output MS after filter
function filter=filterAM(nonfilter,threshold)
         x=nonfilter(:,6);
         y=nonfilter(:,7);
         Am=sqrt((x.^2)+(y.^2));
         order= Am > threshold(1) & Am < threshold(2);
         filter= nonfilter(order,:);
         
end