function [cell_frequencydata,cell_amplitudedata,cell_durationdata] = O_analyze_cell_syn_Events( plot_it,ignore_trial,save_matfile,filter_window )

%ignore_trial - see input info under sequence_syn_frequency.m
%Alectryon:
top = ('/Users/vhlab/Documents/Rig_3/data');
%osik_mac:
%top = ('Users/osik_mac/Documents/MATLAB/data');
%p = addpath(genpath(top));
%save_location = ('/Users/vhlab/Documents/Rig_3/data/completed_analysis');

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
    cell_frequencydata{i,1} = struct('cell_name',[],...
        'cell_ages',[],...
        'cell_experience',[],...
        'cell_total_E_fmean',[],...         %derived from deconvolution method
        'alt_cell_total_E_fmean',[],...     %derived from first deriv method
        'cell_total_I_fmean',[],...         %derived from deconvolution method
        'alt_cell_total_I_fmean',[],...     %derived from first deriv method
        'cell_total_E_IEI',[],...
        'alt_cell_total_E_IEI',[],...
        'cell_total_I_IEI',[],...
        'alt_cell_total_I_IEI',[],...
        'num_Eseqs_dir',[],...
        'num_Iseqs_dir',[],...
        'E_total_int_time_dir',[],...
        'I_total_int_time_dir',[]);
    cell_amplitudedata{i,1} = struct('cell_name',[],...
        'cell_ages',[],...
        'cell_experience',[],...
        'E_trial_deconv_amps',[],...      %column 1 means, column 2 medians
        'E_trial_alt_amps',[],...           %see line 44 comment
        'I_trial_deconv_amps',[],...
        'I_trial_alt_amps',[]);
    cell_durationdata{i,1} = struct('cell_name',[],...
        'cell_ages',[],...
        'cell_experience',[],...
        'E_trial_th_durations',[],...     %see line 44 comment
        'E_trial_fwhm',[],...               %see line 44 comment
        'I_trial_th_durations',[],...
        'I_trial_fwhm',[],...
        'E_trial_risetime',[],...
        'I_trial_risetime',[],...
        'E_trial_decaytau',[],...
        'I_trial_decaytau',[]);
    end

for i = 1:size(directories,1),
    prompt = ['Enter age of animal for cell ',fullfile(top,slice_folder,directories(i,:)),':'];
    a = inputdlg(prompt);
    prompt_ = ['Enter duration of animal experience associated with cell ',...
        fullfile(top,slice_folder,directories(i,:)),':'];
    e = inputdlg(prompt_);
    cell_frequencydata{i,1}.cell_name = [slice_folder,'_',directories(i,:)];
    cell_frequencydata{i,1}.cell_ages = str2num(cell2mat(a));
    cell_frequencydata{i,1}.cell_experience = str2num(cell2mat(e));
    cell_amplitudedata{i,1}.cell_name = [slice_folder,'_',directories(i,:)];
    cell_amplitudedata{i,1}.cell_ages = str2num(cell2mat(a));
    cell_amplitudedata{i,1}.cell_experience = str2num(cell2mat(e));
    cell_durationdata{i,1}.cell_name = [slice_folder,'_',directories(i,:)];
    cell_durationdata{i,1}.cell_ages = str2num(cell2mat(a));
    cell_durationdata{i,1}.cell_experience = str2num(cell2mat(e));
    [E_sequence,I_sequence,count_E,count_I,seq_durationE,seq_durationI,num_devs_E,num_devs_I]=...
        O_sequence_syn_Events(highest_level,plot_it,ignore_trial,filter_window,{strtrim(directories(i,:))});
    cell_frequencydata{i,1}.num_Eseqs_dir = count_E;
    cell_frequencydata{i,1}.num_Iseqs_dir{i,1} = count_I;
    cell_frequencydata{i,1}.E_total_int_time_dir = sum(seq_durationE);
    cell_frequencydata{i,1}.I_total_int_time_dir = sum(seq_durationI);
    cell_total_E_fmean_A = cell2mat(E_sequence.collated_means);
    cell_total_I_fmean_A = cell2mat(I_sequence.collated_means);
    alt_cell_total_E_fmean_A = cell2mat(E_sequence.alt_collated_means);
    alt_cell_total_I_fmean_A = cell2mat(I_sequence.alt_collated_means);
    cell_frequencydata{i,1}.cell_total_E_fmean = cell_total_E_fmean_A(:,1);
    cell_frequencydata{i,1}.cell_total_E_fmean = cell_total_E_fmean_A(:,(2:2:end));
    cell_frequencydata{i,1}.cell_total_I_fmean = cell_total_I_fmean_A(:,1);
    cell_frequencydata{i,1}.cell_total_I_fmean = cell_total_I_fmean_A(:,(2:2:end));
    cell_frequencydata{i,1}.alt_cell_total_E_fmean = alt_cell_total_E_fmean_A(:,1);
    cell_frequencydata{i,1}.alt_cell_total_E_fmean = alt_cell_total_E_fmean_A(:,(2:2:end));
    cell_frequencydata{i,1}.alt_cell_total_I_fmean = alt_cell_total_I_fmean_A(:,1);
    cell_frequencydata{i,1}.alt_cell_total_I_fmean = alt_cell_total_I_fmean_A(:,(2:2:end));
    cell_frequencydata{i,1}.cell_total_E_IEI = cell2mat(reshape(E_sequence.collated_IEI,length(E_sequence.collated_IEI),1));
    cell_frequencydata{i,1}.cell_total_I_IEI = cell2mat(reshape(I_sequence.collated_IEI,length(I_sequence.collated_IEI),1));
    cell_frequencydata{i,1}.alt_cell_total_E_IEI = cell2mat(reshape(E_sequence.alt_collated_IEI,length(E_sequence.alt_collated_IEI),1));
    cell_frequencydata{i,1}.alt_cell_total_I_IEI = cell2mat(reshape(I_sequence.alt_collated_IEI,length(I_sequence.alt_collated_IEI),1));
    cell_amplitudedata{i,1}.E_trial_deconv_amps = cell2mat(reshape(E_sequence.collated_amps,length(E_sequence.collated_amps),1));
    cell_amplitudedata{i,1}.I_trial_deconv_amps = cell2mat(reshape(I_sequence.collated_amps,length(I_sequence.collated_amps),1));
    cell_amplitudedata{i,1}.E_trial_alt_amps = cell2mat(reshape(E_sequence.collated_alt_amps,length(E_sequence.collated_alt_amps),1));
    cell_amplitudedata{i,1}.I_trial_alt_amps = cell2mat(reshape(I_sequence.collated_alt_amps,length(I_sequence.collated_alt_amps),1));
    cell_durationdata{i,1}.E_trial_th_durations = cell2mat(reshape(E_sequence.collated_th_durations,length(E_sequence.collated_th_durations),1));
    cell_durationdata{i,1}.I_trial_th_durations = cell2mat(reshape(I_sequence.collated_th_durations,length(I_sequence.collated_th_durations),1));
    cell_durationdata{i,1}.E_trial_fwhm = cell2mat(reshape(E_sequence.collated_fwhm,length(E_sequence.collated_fwhm),1));
    cell_durationdata{i,1}.I_trial_fwhm = cell2mat(reshape(I_sequence.collated_fwhm,length(I_sequence.collated_fwhm),1));
    cell_durationdata{i,1}.E_trial_risetime = cell2mat(reshape(E_sequence.collated_risetimes,length(E_sequence.collated_risetimes),1));
    cell_durationdata{i,1}.I_trial_risetime = cell2mat(reshape(I_sequence.collated_risetimes,length(I_sequence.collated_risetimes),1));
    cell_durationdata{i,1}.E_trial_decaytau = cell2mat(reshape(E_sequence.collated_decay_taus,length(E_sequence.collated_decay_taus),1));
    cell_durationdata{i,1}.I_trial_decaytau = cell2mat(reshape(I_sequence.collated_decay_taus,length(I_sequence.collated_decay_taus),1));
    %cell_total_E_fcount = 
    %cell_total_I_fcount = 
end

%if plot_fanofactor == 1,   %%no justification for synaptic current inputs
%    FF_timesteps = 0.001;
%    FF_maxwindow = 1.0;

for j = 1:size(cell_frequencydata),
    bin_edges = [0:0.001:0.1];
    D_E = histc(cell_frequencydata{j,1}.cell_total_E_IEI,bin_edges);
    D_I = histc(cell_frequencydata{j,1}.cell_total_I_IEI,bin_edges);
    bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2;
    D_E = D_E(1:end-1);
    D_I = D_I(1:end-1);
    M_y = max([D_E;D_I]);
    f = figure;
    subplot(1,2,1);
    bar(bin_centers,D_E,'b');
    hold on;
    xlabel('Excitatory Inter-event interval (in sec.)');
    ylabel('Count');
    ylim([0 M_y+10]);
    subplot(1,2,2);
    bar(bin_centers,D_I,'r');
    xlabel('Inhibitory Inter-event interval (in sec.)');
    ylabel('Count');
    ylim([0 M_y+10]);
    title(['IEI histogram: Deconvolution Method, E detect:',num2str(num_devs_E),', I detect:',num2str(num_devs_I)]);
    hold off;
    %need to add save file commands here
    saveas(gcf,[cell_frequencydata{j,1}.cell_name,'_cell_deconv_IEIhist','filter_w_',num2str(filter_window)],'fig');
    close(f);
end
for k = 1:size(cell_frequencydata),
    bin_edges = [0:0.001:0.1];
    D_E = histc(cell_frequencydata{k,1}.alt_cell_total_E_IEI,bin_edges);
    D_I = histc(cell_frequencydata{k,1}.alt_cell_total_I_IEI,bin_edges);
    bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2;
    D_E = D_E(1:end-1);
    D_I = D_I(1:end-1);
    M_y = max([D_E;D_I]);
    f = figure;
    subplot(1,2,1);
    bar(bin_centers,D_E,'b');
    hold on;
    xlabel('Excitatory Inter-event interval (in sec.)');
    ylabel('Count');
    ylim([0 M_y+10]);
    subplot(1,2,2);
    bar(bin_centers,D_I,'r');
    xlabel('Inhibitory Inter-event interval (in sec.)');
    ylabel('Count');
    ylim([0 M_y+10]);
    title(['IEI histogram: Derivative Method, E detect:',num2str(num_devs_E),', I detect:',num2str(num_devs_I)]);
    hold off;
    %need to add save file commands here
    saveas(gcf,[cell_frequencydata{j,1}.cell_name,'_cell_cderiv_IEIhist','filter_w_',num2str(filter_window)],'fig');
    close(f);
end

for m = 1:size(cell_durationdata),
    bin_edges = [0:0.001:0.05];
    DR_E = histc(cell_durationdata{m,1}.E_trial_risetime,bin_edges);
    DR_I = histc(cell_durationdata{m,1}.I_trial_risetime,bin_edges);
    bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2;
    DR_E = DR_E(1:end-1);
    DR_I = DR_I(1:end-1);
    MR_y = max([DR_E;DR_I]);
    f= figure;
    subplot(1,2,1);
    bar(bin_centers,DR_E,'b');
    hold on;
    xlabel('Excitatory 10-90% risetime (sec.)');
    ylabel('Count');
    ylim([0 MR_y+10]);
    subplot(1,2,2);
    bar(bin_centers,DR_I,'r');
    xlabel('Inhibitory 10-90% risetime (sec.)');
    ylabel('Count');
    ylim([0 MR_y+10]);
    title('sPSC event rise times');
    hold off;
    saveas(gcf,[cell_durationdata{m,1}.cell_name,'_cell_1090_risetime.fig']);
    close(f);
    bin_edges_ = [0:0.001:0.1];
    DD_E = histc(cell_durationdata{m,1}.E_trial_decaytau,bin_edges_);
    DD_I = histc(cell_durationdata{m,1}.I_trial_decaytau,bin_edges_);
    bin_centers_ = (bin_edges_(1:end-1)+bin_edges_(2:end))/2;
    DD_E = DD_E(1:end-1);
    DD_I = DD_I(1:end-1);
    MD_y = max([DD_E;DD_I]);
    f2 = figure;
    subplot(1,2,1);
    bar(bin_centers_,DD_E,'b');
    hold on;
    xlabel('Excitatory decay time constant (sec.)');
    ylabel('Count');
    ylim([0 MD_y+10]);
    subplot(1,2,2);
    bar(bin_centers_,DD_I,'r');
    xlabel('Inhibitory decay time constant (sec.)');
    ylabel('Count');
    ylim([0 MD_y+10]);
    title('sPSC event decay time constants');
    hold off;
    saveas(gcf,[cell_durationdata{m,1}.cell_name,'_cell_decaytau.fig']);
    close(f2);
end

if save_matfile == 1,
    curr_filename = [cell_frequencydata{j,1}.cell_name,'_adddur_analysis_vars','_filter_w_',num2str(filter_window),'.mat'];
    save(curr_filename,'cell_frequencydata');
    save(curr_filename,'cell_amplitudedata','-append');
    save(curr_filename,'cell_durationdata','-append');
else
end 
    
end
    

