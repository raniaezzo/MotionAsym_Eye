%% this function is used to compute velocity
function v = computevelocity(d,samplingRateData)
         v = zeros(size(d));      
         %%calculate mean velocity for x,y seperately, and then sampliing
         %%window = 3         (d(5:end,:) - d(2:end-3,:) + d(4:end-1,:) - d(1:end-4,:))); %
         %%calculate begin and end of segment, using window = 2
         v(3:length(d)-2,:) = samplingRateData/3*(0.5*(d(5:end,:) - d(2:end-3,:) + d(4:end-1,:) - d(1:end-4,:)));
         v(2,:) = samplingRateData/2*(d(3,:) - d(1,:));
         v(length(d)-1,:) = samplingRateData/2*(d(end,:) - d(end-2,:));
end