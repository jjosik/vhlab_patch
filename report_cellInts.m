function[] = report_cellInts()

%addpath(pwd);
pathFolder = uigetdir();

d = dir(pathFolder);
isub = [d(:).isdir];
test_subFold = {d(isub).name}';
test_subFold(ismember(test_subFold,{'.','..'})) = [];

[s,v] = listdlg('PromptString','Select folders:',...
    'SelectionMode','multiple',...
    'ListString',test_subFold);

cell_directories = test_subFold(s);
matlab_root = pwd;

col_apThreshold = [];
col_apMaxslope = [];

for j = 1:length(cell_directories),
    cd(pathFolder);
    current_dir = cell_directories{j};
    set_dir = dir(current_dir);
    cd([pathFolder filesep current_dir]);
    for jj = 1:length(set_dir),
        if (set_dir(jj).isdir)&&strcmp(set_dir(jj).name,'Patch'),
            [sorted_ISIs,sorted_spiketimes,sorted_peakVm,sorted_apThreshold,...
                sorted_apMaxslope,sorted_apMaxslopeLoc,sorted_vmbaseline,sorted_apAhp,...
                alt_command,current_rep] = get_Intprops(set_dir);
            col_apThreshold = [col_apThreshold;cell2mat(sorted_apThreshold)];
            col_apMaxslope = [col_apMaxslope;cell2mat(sorted_apMaxslope)];
            
        else
        end
    end
end

end
