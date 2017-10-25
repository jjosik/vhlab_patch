function [ data_segment,tdata_segment ] = anaylze_iterations(highest_level,ignore_trial,top_folder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

comp_ignore_trials = strtrim(ignore_trial);

if highest_level == 1,
    %addpath(fullfile(pwd));     %begin at cell list under approp. day/slice #
    main_folder = uigetdir();   %choose sequence directory, e.g. VC_10sec_1_002
    d = dir(main_folder);
    isub = [d(:).isdir];
    sel_Fnames = {d(isub).name}';
    sel_Fnames(ismember(sel_Fnames,{'.','..','Patch'})) = [];
    [s,v] = listdlg('PromptString','Select folders:',...
                'SelectionMode','multiple',...
                'ListString',sel_Fnames);

    new_directories = char(sel_Fnames(s));
    cd(main_folder);
else
    sel_Fnames = [];
    for i=1:size(top_folder,1),
        main_folder = top_folder;   %choose sequence directory, i.e. VC_10sec_1_002
        d = dir(main_folder);
        isub = [d(:).isdir];
        unform_Fnames = {d(isub).name}';
        unform_Fnames(ismember(unform_Fnames,{'.','..','Patch','hide'})) = [];
        upper_fols = {repmat({main_folder},numel(unform_Fnames),1)};
        fsep_symb = char('/');
        upper_fols_comp = {strcat(upper_fols{:,1},repmat(fsep_symb,length(upper_fols{:,1}),1))};
        sel_Fnames = [sel_Fnames;strcat(upper_fols_comp{i},unform_Fnames)];
        cd(main_folder);
    end
end
for i=1:length(sel_Fnames),
    if i == 1,
        new_directories{i,1} = char(sel_Fnames{i,1});
    else
        new_directories{i,1} = char(unform_Fnames{i,1});
    end
end
%new_directories = char(sel_Fnames{:});

data_segment = [];
tdata_segment = [];
rep_count = 0;
for sub_d=1:size(new_directories,1),
    d_flag = 0;
    subd_string = num2str(sub_d);
    if ~isempty(regexp(comp_ignore_trials,subd_string,'once')),
        continue
    else
    end
    try
        cd(new_directories{sub_d,1});
        d_flag = 1;
    catch
        cd(char(unform_Fnames{sub_d,1}));
        d_flag = 1;
    end
    search_d = dir(pwd);
    isub_D = [search_d(:).isdir];
    sub_s = {search_d(isub_D).name}';
    sub_s(ismember(sub_s,{'.','..','Patch','hide'})) = [];
    rep_count = rep_count + 1;
    if ~isempty(find(isub_D)),
        for j = 1:length(sub_s),
            subs_string = char(sub_s(j,1));
            cd(subs_string);
            tdata = hdf5read('Clamp1.ma','/info/1/values');
            data = hdf5read('Clamp1.ma','/data');
            sampling_rate = hdf5read('Clamp1.ma','info/2/DAQ/primary/rate');
            data_segment = [data_segment;data(:,2)];
            spacer = repmat((tdata(end,1)*(j-1)),length(tdata),1);
            if rep_count == 1,
                trial_spacer = zeros(length(tdata),1);
            else
            end
            tdata_segment = [tdata_segment;(trial_spacer+spacer+tdata(:))];
            if j == length(sub_s),
                trial_mark = tdata_segment(end,1);
                trial_spacer = repmat(trial_mark,length(tdata),1);
            else
            end
            clear spacer;
            cd ../
        end
        if d_flag == 1,
            cd ../
        end
    else
        tdata = hdf5read('Clamp1.ma','/info/1/values');
        data = hdf5read('Clamp1.ma','/data');
        sampling_rate = hdf5read('Clamp1.ma','info/2/DAQ/primary/rate');
        data_segment = [data_segment;data(:,2)];
        tdata_segment = [tdata_segment;tdata(:)];
        if d_flag == 1,
            cd ../
        end
    end

end


    
    
end


