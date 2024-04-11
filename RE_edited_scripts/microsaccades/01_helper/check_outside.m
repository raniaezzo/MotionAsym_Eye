%% check the outside points 
%this function can save the clear data and output the number of points for
%each trial, and the input is the filepath for new Dat_stim data
function check_outside(filepath, screenCenter, dvaPerPx)

    thresh = 1.5/dvaPerPx; % this labels all data outside of 1.5 degrees from center

    [path, filename,~] = fileparts(filepath);
    load(filepath);
    aa=size(Dat_stim,1);
    outside=zeros(aa,2);
    for jj = 1 : aa 
        outside(jj,1)=(Dat_stim(jj,2)-screenCenter(1))^2+(Dat_stim(jj,3)-screenCenter(2))^2 > thresh^2;
    end
    outside(:,2)=Dat_stim(:,1);

    filelabel = strrep(filename, 'Dat_stim', 'outside');
    savename = fullfile(path, filelabel);
    save(savename,'outside');
end