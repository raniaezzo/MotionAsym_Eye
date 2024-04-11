%% create a new tab file, which contain behavior data
%filepath1= tab filepath way/ filepath2= csv filepath way/
%filepath3 = Dat_stim file path
function converge(filepath1,filepath2)
    [path,filename, ~] = fileparts(filepath1);
    tab_stim=load(filepath1);
    M_csv=csvread(filepath2);
    %Dat_stim=load(filepath3);
    tab=tab_stim.tab; 

    % alternative ways to filter (in case experiment ended unexpectedly)
    try
        % alternative: index the tab file based in the behavioral file
        % (this used to be the catch but it's less prone to issues)
%         M_csv = M_csv(M_csv(:,2)>0,:);
%         a = M_csv(:,2);
%         tab=tab(a,:);
        filter = M_csv(:,2)>0;
        tab=tab(filter,:);
        M_csv=M_csv(filter,:);

    catch
        % index the behavioral file based on the tab file
        a=tab(:,7)==0; % get indices for fixation / trial breaks
        %a = (tab(:, 6) ~= 0) & (tab(:, 7) == 0); % extra logic 
        tab=tab(a,:); % create tab file without fixation break segments
        M_csv=M_csv(a,:); % filter out rows from csv that correspond to fixation break in tab file
    end
    
    c=M_csv(:,3);
    d=M_csv(:,4);
    f=M_csv(:,6);
    j=M_csv(:,10);
    m=M_csv(:,13);
    n=M_csv(:,14);
    %tab_new=[tab c d f j m n];
    tab=[tab c d f j m n];
    %savename = fullfile(outputpath, sprintf('%s_new.mat',filename));
    savename = fullfile(path, filename);
    save(savename,'tab');
end

