function [sorted_ISIs,sorted_spiketimes,sorted_peakVm,sorted_apThreshold,...
    sorted_apMaxslope,sorted_apMaxslopeLoc,sorted_vmbaseline,sorted_apAhp,...
    alt_command,current_rep] = get_Intprops(data_folder)

step_folders = dir(data_folder);
step_folders(ismember({step_folders.name}',{'.','..','Patch','hide'})) = [];
step_folders(find([step_folders.isdir]==0)) = [];
current_rep = [];
root = pwd;

step_intStruct = {};
for i=1:length(step_folders),
    if (step_folders(i).isdir)&&~strcmp(step_folders(i).name,'.')&&~strcmp(step_folders(i).name,'..')&&~strcmp(step_folders(i).name,'Patch'),
        current_step = step_folders(i).name;
        cd([root filesep data_folder filesep current_step]);
        
        data = hdf5read('Clamp1.ma','/data');
        tdata = hdf5read('Clamp1.ma','/info/1/values');
        sampling_rate = hdf5read('Clamp1.ma','/info/2/DAQ/primary/rate');
        
        current_rep(i,1) = mean(data((round(0.01*sampling_rate)+1):(round((0.01+0.3)*sampling_rate)),1));
        step_intStruct{end+1,1} = open('step_AP_properties.mat');
    else
    end
end

step_ISIs = {};
step_spiketimes = {};
step_peakVm = {};
step_apThreshold = {};
step_apMaxslope = {};
step_apMaxslopeLoc = {};
step_vmbaseline = [];
step_apAhp = [];
for j=1:length(step_intStruct),
    step_ISIs{end+1,1} = step_intStruct{j,1}.step_AP_properties.train_ISIs;
    step_spiketimes{end+1,1} = step_intStruct{j,1}.step_AP_properties.spike_times;
    step_peakVm{end+1,1} = step_intStruct{j,1}.step_AP_properties.peak_vm;
    step_apThreshold{end+1,1} = step_intStruct{j,1}.step_AP_properties.ap_threshold;
    step_apMaxslope{end+1,1} = step_intStruct{j,1}.step_AP_properties.ap_maxslope;
    step_apMaxslopeLoc{end+1,1} = step_intStruct{j,1}.step_AP_properties.ap_maxslope_loc;
    step_vmbaseline(end+1,1) = step_intStruct{j,1}.step_AP_properties.vm_baseline;
    step_apAhp(end+1,1) = step_intStruct{j,1}.step_AP_properties.ap_ahp;
end
clear step_intStruct;
%output_Struct = struct('sorted_ISIs',{},'sorted_spiketimes',{},'sorted_peakVm',...
%{},'sorted_apThreshold',{},'sorted_apMaxslope',{},'sorted_apMaxslopeLoc',{},...
%    'sorted_vmbaseline',[],'sorted_apAhp',[],'alt_command',[]);
%sort stim. order
[sorted_command,sort_indices] = sort(current_rep);
sorted_ISIs = {};
sorted_spiketimes = {};
sorted_peakVm = {};
sorted_apThreshold = {};
sorted_apMaxslope = {};
sorted_apMaxslopeLoc = {};
sorted_vmbaseline = {};
sorted_apAhp = {};
for k=1:length(sort_indices),
    sorted_ISIs{end+1,1} = step_ISIs{sort_indices(k,1),1};
    sorted_spiketimes{end+1,1} = step_spiketimes{sort_indices(k,1),1};
    sorted_peakVm{end+1,1} = step_peakVm{sort_indices(k,1),1};
    sorted_apThreshold{end+1,1} = step_apThreshold{sort_indices(k,1),1};
    sorted_apMaxslope{end+1,1} = step_apMaxslope{sort_indices(k,1),1};
    sorted_apMaxslopeLoc{end+1,1} = step_apMaxslopeLoc{sort_indices(k,1),1};
    sorted_vmbaseline{end+1,1} = step_vmbaseline(sort_indices(k,1),1);
    sorted_apAhp{end+1,1} = step_apAhp(sort_indices(k,1),1);
end
alt_command = 0:(400/(length(sorted_command)-1)):400;

%output_Struct.sorted_ISIs = int_sorted_ISIs;
%output_Struct.sorted_spiketimes = int_sorted_spiketimes;
%output_Struct.sorted_peakVm = int_sorted_peakVm;
%output_Struct.sorted_apThreshold = int_sorted_apThreshold;
%output_Struct.sorted_apMaxslope = int_sorted_apMaxslope;
%output_Struct.sorted_apMaxslopeLoc = int_sorted_apMaxslopeLoc;
%output_Struct.sorted_vmbaseline = int_sorted_vmbaseline;
%output_Struct.sorted_apAhp = int_sorted_apAhp;
%output_Struct.alt_command = int_alt_command;


%cd ../..
% cell plots
%interval_bins = 0:0.001:3;
%peak_bins = -0.01:0.001:0.05;
%th_bins = -0.06:0.001:-0.02;
%base_bins = -0.08:0.001:0.05;
%ahp_bins = 0.001:0.001:0.1;
%cell_ISI_counts = histc(cell2mat(step_ISIs{

end