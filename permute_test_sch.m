
% these indices correspond to 1300-100 --> 1800+200
% vectorsubs_cardinal_easy = vectorsubs_cardinal_easy(:,200:1000);
% vectorsubs_oblique_hard = vectorsubs_oblique_hard(:,200:1000);
% vectorsubs_cardinal_hard = vectorsubs_cardinal_hard(:,200:1000);
% vectorsubs_oblique_easy = vectorsubs_oblique_easy(:,200:1000);

% after running preprocess3.m
figure
plot(vectorsubs_cardinal_easy', 'LineWidth',2)
hold on
plot(vectorsubs_oblique_hard', ':', 'LineWidth',2)
hold off

% cannot use imgaussfilt! It filters in both dimensions of the image
windowSize = 10; %50; 
b = (1/windowSize)*ones(1,windowSize); a = 1;
vectorsubs_cardinal_easy_sm = filter(b,a,vectorsubs_cardinal_easy'); %imgaussfilt(vectorsubs_cardinal_hard,3);
vectorsubs_cardinal_hard_sm = filter(b,a,vectorsubs_cardinal_hard'); %imgaussfilt(vectorsubs_cardinal_easy,3);
vectorsubs_oblique_easy_sm = filter(b,a,vectorsubs_oblique_easy'); %imgaussfilt(vectorsubs_oblique_hard,3);
vectorsubs_oblique_hard_sm = filter(b,a,vectorsubs_oblique_hard'); %imgaussfilt(vectorsubs_oblique_easy,3);

figure
plot(vectorsubs_cardinal_easy_sm, 'LineWidth',2)
hold on
plot(vectorsubs_oblique_hard_sm, ':', 'LineWidth',2)
hold off

vectorsubs_cardinal_easy_sm = vectorsubs_cardinal_easy_sm';
vectorsubs_oblique_hard_sm = vectorsubs_oblique_hard_sm';

%% ---------------------------------cluster permutation test----------------------------------------------%
% data1 and data2 are organized as n x m matrices
% in which n is "datapoints" and m is "subjects"

% cond1 = [vectorsubs_cardinal_easy + vectorsubs_oblique_easy]./2;
% cond2 = [vectorsubs_cardinal_hard + vectorsubs_oblique_hard]./2;
% 
% %cond1 = [vectorsubs_cardinal_easy + vectorsubs_cardinal_hard]./2;
% %cond2 = [vectorsubs_oblique_easy + vectorsubs_oblique_hard]./2;
%  
% cond1 = imgaussfilt(cond1,3);
% cond2 = imgaussfilt(cond2,3);

% % cardinal easy vs oblique hard
cond1 = vectorsubs_cardinal_easy_sm;
cond2 = vectorsubs_oblique_hard_sm;

figure
plot(cond1', 'LineWidth',2)
hold on
plot(cond2', ':', 'LineWidth',2)
hold on

data1 = cond1'; %all_MSrate_pre'
data2 = cond2'; %all_MSrate_post1'; 
% run cluster permutation test
dependent_samples = true;
p_threshold = 0.05;
num_permutations = 1000; % 1000
two_sided = false; %true;
num_clusters =[];
[clusters, p_values, t_sums, permutation_distribution ] = permutest( data1, data2, dependent_samples, ...
    p_threshold, num_permutations, two_sided, num_clusters );
% find significant cluster(s)
sig_cluster = clusters(p_values<p_threshold);% 2 clusters for 24 observers

xline(300) % this is due to the segment arrangement
hold on
xline(800)
hold on
% xline(100)
% hold on
% xline(600)

for pi=1:numel(sig_cluster)
    x = sig_cluster{pi};
    plot(sig_cluster{pi}, ones(1,length(x))*3, 'k', 'LineWidth', 2)
    hold on
end

ylim([0 3.2])
xlim([200 1000])
% xlim([])