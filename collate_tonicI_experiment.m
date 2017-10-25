function [ cell_data_vector,cell_tdata_vector,p,on_time,off_time ] = collate_tonicI_experiment( plot_it,ignore_trial )

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

duration = 10.01;
iteration_num = 14;
cell_data_vector = [];
cell_tdata_vector = [];
for i = 1:size(directories,1),
    top_folder = cell2mat({strtrim(directories(i,:))});
    [ data_segment,tdata_segment ] = anaylze_iterations(highest_level,ignore_trial,top_folder);
    if i > 1,
        cell_data_vector = [cell_data_vector;data_segment];
        translate_tdata = tdata_segment(:,1)+(duration*iteration_num*(i-1));
        cell_tdata_vector = [cell_tdata_vector;translate_tdata];
        clear data_segment tdata_segment translate_tdata;
    else
        cell_data_vector = data_segment;
        cell_tdata_vector = tdata_segment;
    end
end

if plot_it == 1,
    prompt_a = ['Enter age of animal for this cell:'];
    a = str2num(char(inputdlg(prompt_a)));
    prompt_b = ['Enter number directory when GABA was turned ON:'];
    b = str2num(char(inputdlg(prompt_b)));
    prompt_c = ['Enter number of iteration under this directory when GABA was turned ON:'];
    c = str2num(char(inputdlg(prompt_c)));
    prompt_d = ['Enter number of directory when GABA was turned OFF:'];
    d = str2num(char(inputdlg(prompt_d)));
    prompt_e = ['Enter number of iteration under this directory when GABA was turned OFF:'];
    e = str2num(char(inputdlg(prompt_e)));
    on_time = (b-1)*(duration*iteration_num)+(duration*c);
    off_time = (d-1)*(duration*iteration_num)+(duration*e);
    f = figure;
    p = plot(cell_tdata_vector(:,1),cell_data_vector(:,1));
    xlabel('Time (seconds)');
    ylabel('Current (Amps)');
    %title('
    hold on;
    line([on_time off_time],[max(cell_data_vector)+(0.1*max(cell_data_vector)) max(cell_data_vector)+(0.1*max(cell_data_vector))]);
    hold off;
else
end



    
    

end

