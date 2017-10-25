function [ seq_total_E_charge,seq_total_I_charge,seq_durationE,seq_durationI,count_E,count_I ] = sequence_syn_input_charge_A...
    ( highest_level,plot_it,ignore_trial,varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

comp_ignore_trials = strtrim(ignore_trial);

if highest_level == 1,
    %addpath(fullfile(pwd));     %begin at cell list under approp. day/slice #
    main_folder = uigetdir();   %choose sequence directory, i.e. VC_10sec_1_002
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
    
    for i=1:length(varargin),
        main_folder = varargin{i}; %provides a list of vc sequence dirs. to be analyzed
        d = dir(cell2mat(main_folder));
        isub = [d(:).isdir];
        unform_Fnames = {d(isub).name}';
        unform_Fnames(ismember(unform_Fnames,{'.','..','Patch','hide'})) = [];
        upper_fols = {repmat({main_folder},numel(unform_Fnames),1)};
        fsep_symb = char('/');
        upper_fols_comp = {strcat(upper_fols{:,1},repmat(fsep_symb,length(upper_fols{:,1}),1))};
        sel_Fnames = [sel_Fnames;strcat(upper_fols_comp{i},unform_Fnames)];
    end
    for i=1:length(sel_Fnames),
        new_directories{i,1} = char(cell2mat(sel_Fnames{i,1}));
    end
    %new_directories = char(sel_Fnames{:});
    
end

loop_count = 0;
count_E = 0;
count_I = 0;
seq_total_E_charge = [];
seq_total_I_charge = [];
seq_durationE = [];
seq_durationI = [];
samp_rate_array = [];
for fd =1:size(new_directories,1),
    fd_string = num2str(fd);
    if ~isempty(regexp(comp_ignore_trials,fd_string,'once')),
        continue
    else
    end
    loop_count = loop_count+1;
    next_level = strtrim(new_directories(fd,:));
    alt_next_level = strtrim(char(unform_Fnames(fd,:)));
    if loop_count == 1,
        dn = dir(char(next_level));
    else
        dn = dir(alt_next_level);
    end
    minsub = [dn(:).isdir];
    min_Fnames = {dn(minsub).name}';
    min_Fnames(ismember(min_Fnames,{'.','..','Patch','hide'})) = [];
    min_dirs = char(min_Fnames{:});
    if loop_count == 1,
        cd(char(next_level));
    else
        cd(alt_next_level);
    end
    
    for j = 1:numel(min_Fnames),
        cd(strtrim(min_dirs(j,:)));
        if exist([pwd filesep 'Clamp1_uncomp.ma'],'file'),
            use_medwin = 1; filter_win = 5.0;
            [epoch_charge,sampling_rate,is_E_test,acq_duration] = compute_syn_input_TRAP_medwin(plot_it,use_medwin,filter_win);
            samp_rate_array(end+1) = sampling_rate;
            if is_E_test == 1,
                seq_total_E_charge(end+1,1) = epoch_charge;
                count_E = count_E + 1;
                seq_durationE(end+1,1) = acq_duration;
            else
                seq_total_I_charge(end+1,1) = epoch_charge;
                count_I = count_I + 1;
                seq_durationI(end+1,1) = acq_duration;
            end
        else
            fol_string = strcat(strcat(new_directories(fd,:),fsep_symb),min_dirs(j,:));
            qstring = (['Error finding analysis file in folder: ', fol_string,'.  Exclude this directory from analysis?']);
            choice = questdlg(qstring,'Resolution','Yes','No','Yes');
            if choice == 'Yes',
                samp_rate_array(end+1) = NaN;
                cd ..
                continue;
            else
                break;
            end
        end
        cd ..
    end
    cd ..
end

%if ~isnan(samp_rate_array),
%    if ~all(samp_rate_array == samp_rate_array(1)),
%        fprintf('%s \n','WARNING: sampling rate mismatch in one or more analyzed files.');
%        samp_warning = 1;
%    else
%        sequence_sampling_rate = samp_rate_array(end);
%    end
%else
%end



end


