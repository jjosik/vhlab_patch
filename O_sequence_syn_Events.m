function [E_sequence,I_sequence,count_E,count_I,seq_durationE,seq_durationI,num_devs_E,num_devs_I] = ...
    O_sequence_syn_Events( highest_level,plot_it,ignore_trial,filter_window,varargin )
%ignore_trial - string input specifies which trial numbers should be
%excluded from analysis; note, use whitespace to separate each entry and if 
%trial number >= 10 (i.e. contains more than a single digit) 
%the entry within the string variable must be in square brackets,
%e.g. ignore_trial = '1 9 [10]' will analyze all trial runs except 1, 9 and
%10.
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
    
    
end

loop_count = 0;
count_E = 0;
count_I = 0;
E_sequence = struct('collated_counts',[],'collated_means',[],'collated_IEI',[],...
    'collated_amps',[],'collated_alt_amps',[],'collated_th_durations',[],'collated_fwhm',[],...
    'alt_collated_counts',[],'alt_collated_means',[],'alt_collated_IEI',[],'collated_risetimes',[],...
    'collated_decay_taus',[]);
I_sequence = struct('collated_counts',[],'collated_means',[],'collated_IEI',[],...
    'collated_amps',[],'collated_alt_amps',[],'collated_th_durations',[],'collated_fwhm',[],...
    'alt_collated_counts',[],'alt_collated_means',[],'alt_collated_IEI',[],'collated_risetimes',[],...
    'collated_decay_taus',[]);
seq_durationE = [];
seq_durationI = [];
samp_rate_array = [];
for fd = 1:size(new_directories,1),
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
            save_it = 1;  %consider passing these from analyze_cell... file (enter as input there)
            bin_sizes = [0.01,0.1,1];  %consider passing from analyze_cell...
            num_devs_E = 1.5;
            num_devs_I = 1.5;
            s_margin = 0.002;
            windowed_baseline = 1;
            th_follows_base = 1;        %***IMPORTANT - following function modified to _EDIT version 5/4/17
            [frequency_data,amplitude_data,duration_data,sampling_rate,is_E_test,acq_duration] = ...
                O_PSC_deconvolution_EDIT(plot_it,save_it,bin_sizes,num_devs_E,num_devs_I,s_margin,...
                windowed_baseline,th_follows_base,filter_window);
            samp_rate_array(end+1) = sampling_rate;
            if is_E_test == 1,
                E_sequence.collated_counts{1,end+1} = frequency_data.binned_counts;
                E_sequence.alt_collated_counts{1,end+1} = frequency_data.alt_binned_counts;
                E_sequence.collated_means{1,end+1} = frequency_data.binned_means;
                E_sequence.alt_collated_means{1,end+1} = frequency_data.alt_binned_means;
                E_sequence.collated_IEI{end+1} = frequency_data.epoch_IEI;
                E_sequence.alt_collated_IEI{end+1} = frequency_data.alt_epoch_IEI;
                E_sequence.collated_amps{end+1} = amplitude_data.deconv_trial_amplitudes;
                E_sequence.collated_alt_amps{end+1} = amplitude_data.alt_trial_amplitudes;
                E_sequence.collated_th_durations{end+1} = duration_data.threshold_crossing_durations;
                E_sequence.collated_fwhm{end+1} = duration_data.fwhm_durations;
                E_sequence.collated_risetimes{end+1} = duration_data.risetimes;
                E_sequence.collated_decay_taus{end+1} = duration_data.decay_taus;
                count_E = count_E + 1;
                seq_durationE(end+1,1) = acq_duration;
            else
                I_sequence.collated_counts{1,end+1} = frequency_data.binned_counts;
                I_sequence.alt_collated_counts{1,end+1} = frequency_data.alt_binned_counts;
                I_sequence.collated_means{1,end+1} = frequency_data.binned_means;
                I_sequence.alt_collated_means{1,end+1} = frequency_data.alt_binned_means;
                I_sequence.collated_IEI{end+1} = frequency_data.epoch_IEI;
                I_sequence.alt_collated_IEI{end+1} = frequency_data.alt_epoch_IEI;
                I_sequence.collated_amps{end+1} = amplitude_data.deconv_trial_amplitudes;
                I_sequence.collated_alt_amps{end+1} = amplitude_data.alt_trial_amplitudes;
                I_sequence.collated_th_durations{end+1} = duration_data.threshold_crossing_durations;
                I_sequence.collated_fwhm{end+1} = duration_data.fwhm_durations;
                I_sequence.collated_risetimes{end+1} = duration_data.risetimes;
                I_sequence.collated_decay_taus{end+1} = duration_data.decay_taus;
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