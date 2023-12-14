%to find point 
function [t1,t2,t3,y_fil] = infl(rate,time,pa_pass,temporalthresh,onset,stepsize, sub)


seglength = 49; % was 100 but few errors
dev_thresh = 10;
fac=50; % increase to make more liberal

if strcmp(sub, 'S01')
    seglength = 30;
elseif strcmp(sub, 'S02')
    seglength = 150;
    dev_thresh = 2;
    fac=50;
elseif strcmp(sub, 'S03')
    seglength = 50; %50? % 300?
    fac = 100; %80;
    %dev_thresh = 0;
elseif strcmp(sub, 'S04')
    fac = 75;
elseif strcmp(sub, 'S06')
    seglength = 150;
% elseif strcmp(sub, 'S03')
%     seglength = 100; %50?
%     fac = 80; %80;
%     dev_thresh = 4;
end

% time=[500,2000];
% pa_pass=0.001;
% temporalthresh=70;
% onset=1;
% stepsize=1;

% y_fil_old = lowpass(rate,pa_pass);
y_fil = imgaussfilt(rate, 50); %3);

time1 = time./(2500/length(rate));
temporalthresh=temporalthresh/(2500/length(rate));
m = min(rate);
d = nan(1, length(y_fil)/stepsize); % derivative
for i=1:(length(y_fil)/stepsize)
    try
        temp = mean(diff(y_fil(onset:onset+temporalthresh)));
        d(i) = temp;
        onset = onset + stepsize;
    catch
        d(i) = d(i);
    end
end

[~, idx] = min(d(time1(1):time1(1)+800));
tt=sqrt(var(rate(time1(1):time1(2))));
t1 = idx + time1(1);
% find where slope is 0 after this
temp = d(t1:time1(2))>0 & rate(t1:time1(2)) - m <tt/2;
t2 = find(temp, 1, 'first') + t1;
%temp2 = d(t2:time1(2))<0 & rate(t2:time1(2)) - m <tt/2; % works for S08

minderiv = ((tt/2)/fac);

if strcmp(sub, 'S01') || strcmp(sub, 'S06')
    minderiv=0;
end

temp2 = d(1300:time1(2))>minderiv & rate(1300:time1(2)) - m >tt/dev_thresh;
%temp2 = d(1300:time1(2))>0.01 & rate(1300:time1(2)) - m >tt/3; %optimal?
%t3 = find(temp2, 1, 'first') + 1300;

% rania added
x = temp2;
comps = bwconncomp(x);
list = comps.PixelIdxList;
%[~,finderVec] = max(cell2mat(cellfun(@size,list,'UniformOutput',false)));
detected_sizes = [];
for pp=1:numel(list)
    [ee, ~] = size(list{pp});
    detected_sizes = [detected_sizes ee];
end
isBig = detected_sizes>=seglength;
%isBig = detected_sizes>=100;

if ~all(isBig)
    %isBig = detected_sizes>=seglength; % use 50 thresh if no 100 
    [~,isBig] = max(detected_sizes);
end

finderVec = find(isBig, 1, 'first');
bing = list{finderVec}(1); % start of segment in timeseries
t3 = bing + 1300; % add the offset outside of the segement for abs measure

% temp3 = d(t3_a:time1(2))>0 & rate(t3_a:time1(2)) - m <tt/2;
% t3_b = find(temp3, 1, 'first') + t2;
% 
% t3 = t3_a+(t3_b-t3_a);

%seg = 200; % check if continually increasing for at least 200 ms


% figure
% plot(rate)
% hold on
% plot(y_fil_old, 'LineWidth',2)
% hold on
% plot(y_fil, 'LineWidth',2)
% hold on
% xline(t1)
% hold on
% xline(t2)
% hold on
% xline(t3, 'r')

% figure
% plot(d)

end