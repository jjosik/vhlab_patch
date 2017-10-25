pathFolder = uigetdir();

d = dir(pathFolder);
isub = [d(:).isdir];
test_subFold = {d(isub).name}';
test_subFold(ismember(test_subFold,{'.','..'})) = [];

[s,v] = listdlg('PromptString','Select folders:',...
    'SelectionMode','multiple',...
    'ListString',test_subFold);

cell_directories = test_subFold(s);

col_apThreshold = [];
col_apMaxslope = [];
cd(pathFolder);

for j = 1:length(cell_directories),
    cd(cell_directories{j});
    cd('cell_int_props');
    load('cell_int_props.mat');
    col_apThreshold = [col_apThreshold;mean(cell2mat(cell_apThreshold))];
    col_apMaxslope = [col_apMaxslope;mean(cell2mat(cell_apMaxslope))];
    cd ../../
end
    