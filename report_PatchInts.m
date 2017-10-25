function[s_cap,s_ir,s_ar] = report_PatchInts()

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

cap_array = {};
IR_array = {};
AR_array = {};
for j = 1:length(cell_directories),
    cd(pathFolder);
    current_dir = cell_directories{j};
    set_dir = dir(current_dir);
    cd([pathFolder filesep current_dir]);
    for jj = 1:length(set_dir),
        if (set_dir(jj).isdir)&&strcmp(set_dir(jj).name,'Patch'),
            [capacitance,input_R,access_R] = get_Patchprops(set_dir(jj).name);
            cap_array{end+1,1} = capacitance;
            IR_array{end+1,1} = input_R;
            AR_array{end+1,1} = access_R;
        else
        end
    end
    
end

s_cap = [];
s_ir = [];
s_ar = [];
for i=1:length(cap_array),
    s_cap(end+1,1) = cap_array{i,1}.value;
    s_ir(end+1,1) = IR_array{i,1}.value;
    s_ar(end+1,1) = AR_array{i,1}.value;
end
