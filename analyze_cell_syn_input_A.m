function [cell_chargedata] = analyze_cell_syn_input_A( plot_it,ignore_trial,save_matfile )
%ANALYZE_CELL_SYN_INPUT_A  Executes the main interface initiating analysis of total
%   collated charge from cellular E/I datasets.
%
%   Inputs:     PLOT_IT  - when set equal to '1' will generate analysis plot outputs
%               IGNORE_TRIAL - documents the subfolders under the cell folder that should
%                   be excluded from analysis (typically I = 0 checks).  E.g., format is '1 [10]',
%                   in this case subfolders 1 and 10 are marked for exclusion.
%               SAVE_MATFILE - when set equal to '1' will save to local directory the outputs, plots,
%                   and workspace for the current analysis.
%   Outputs:    CELL_CHARGEDATA -
%
%
%   NOTE - This function chains to other functions: sequence_syn_input_charge_A.m and compute_syn_input_TRAP_medwin.m
%   REQUISITES:  must begin already cd'd to appropriate date folder (displays all slice
%   options)
%Alectryon:
top = ('/Users/vhlab/Documents/Rig_3/data');
%osik_mac:
%top = ('Users/osik_mac/Documents/MATLAB/data');
p = addpath(genpath(top));

slice_folder = uigetdir();
d = dir(slice_folder);
isub = [d(:).isdir];
sel_Snames = {d(isub).name}';
sel_Snames(ismember(sel_Snames,{'.','..','Patch'})) = [];
[s,v] = listdlg('PromptString','Select folders:',...
    'SelectionMode','multiple',...
    'ListString',sel_Snames);

directories = char(sel_Snames(s));
cd(slice_folder);

highest_level = 0;

for i = 1:size(directories,1),
    cell_chargedata{i,1} = struct('cell_name',[],...
        'cell_ages',[],...
        'cell_experience',[],...
        'cell_total_E_Q',[],...
        'cell_total_I_Q',[],...
        'num_Eseqs',[],...
        'num_Iseqs',[],...
        'E_total_int_time',[],...
        'I_total_int_time',[],...
        'cell_total_E_Q_dir',[],...
        'cell_total_I_Q_dir',[],...
        'num_Eseqs_dir',[],...
        'num_Iseqs_dir',[],...
        'E_total_int_time_dir',[],...
        'I_total_int_time_dir',[]),...
        'ignore_trial';
end

for i = 1:size(directories,1),
    prompt = ['Enter age of animal for cell ',fullfile(top,slice_folder,directories(i,:)),':'];
    a = inputdlg(prompt);
    prompt_ = ['Enter duration of animal experience associated with cell ',...
        fullfile(top,slice_folder,directories(i,:)),':'];
    e = inputdlg(prompt_);
    cell_chargedata{i,1}.cell_name = [slice_folder,'_',directories(i,:)];
    cell_chargedata{i,1}.cell_ages = str2num(cell2mat(a));
    cell_chargedata{i,1}.cell_experience = str2num(cell2mat(e));
    [seq_total_E_charge,seq_total_I_charge,seq_durationE,seq_durationI,count_E,count_I]=sequence_syn_input_charge_A...
        (highest_level,plot_it,ignore_trial,{strtrim(directories(i,:))});
    cell_chargedata{i,1}.num_Eseqs_dir{i,1} = count_E;
    cell_chargedata{i,1}.num_Iseqs_dir{i,1} = count_I;
    cell_chargedata{i,1}.E_total_int_time_dir{i,1} = seq_durationE;
    cell_chargedata{i,1}.I_total_int_time_dir{i,1} = seq_durationI;
    cell_chargedata{i,1}.cell_total_E_Q_dir{i,1} = seq_total_E_charge;
    cell_chargedata{i,1}.cell_total_I_Q_dir{i,1} = seq_total_I_charge;
    cell_chargedata{i,1}.num_Eseqs = sum(cell2mat(cell_chargedata{i,1}.num_Eseqs_dir));
    cell_chargedata{i,1}.num_Iseqs = sum(cell2mat(cell_chargedata{i,1}.num_Iseqs_dir));
    cell_chargedata{i,1}.E_total_int_time = sum(cell2mat(cell_chargedata{i,1}.E_total_int_time_dir));
    cell_chargedata{i,1}.I_total_int_time = sum(cell2mat(cell_chargedata{i,1}.I_total_int_time_dir));
    cell_chargedata{i,1}.cell_total_E_Q = sum(cell2mat(cell_chargedata{i,1}.cell_total_E_Q_dir));
    cell_chargedata{i,1}.cell_total_I_Q = sum(cell2mat(cell_chargedata{i,1}.cell_total_I_Q_dir));
    cell_chargedata{i,1}.ignore_trial = ignore_trial;
end
%if ~(samp_warning == 1),
%    cell_sampling_rate = sequence_sampling_rate;
%else
%    cell_sampling_rate = NaN;
%end
 
 
%cell_chargedata{i,1}.num_Eseqs = sum(cell2mat(cell_chargedata{i,1}.num_Eseqs_dir));
%cell_chargedata{i,1}.num_Iseqs = sum(cell2mat(cell_chargedata{i,1}.num_Iseqs_dir));
%cell_chargedata{i,1}.E_total_int_time = sum(cell2mat(cell_chargedata{i,1}.E_total_int_time_dir));
%cell_chargedata{i,1}.I_total_int_time = sum(cell2mat(cell_chargedata{i,1}.I_total_int_time_dir));
%cell_chargedata{i,1}.cell_total_E_Q = sum(cell2mat(cell_chargedata{i,1}.cell_total_E_Q_dir));
%cell_chargedata{i,1}.cell_total_I_Q = sum(cell2mat(cell_chargedata{i,1}.cell_total_I_Q_dir));

if save_matfile == 1,
    curr_filename = [cell_chargedata{i,1}.cell_name,'_charge_vars.mat'];
    save(curr_filename,'cell_chargedata');
else
end

end