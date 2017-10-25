function [capacitance,input_R,access_R] = get_Patchprops(data_folder)

cd([pwd filesep data_folder]);
input_R = struct('cell',[],'value',[],'flag',[]);
access_R = struct('cell',[],'value',[],'flag',[]);
capacitance = struct('cell',[],'value',[],'flag',[]);

max_file = [];
try
    data = hdf5read('Clamp1_000.ma','/data');
    tdata = hdf5read('Clamp1_000.ma','/info/0/values');
catch
    data = hdf5read('Clamp1_001.ma','/data');
    tdata = hdf5read('Clamp1_001.ma','/info/0/values');
end
if exist('Clamp1_001.ma','file'),
    try
        data2 = hdf5read('Clamp1_001.ma','/data');
        tdata2 = hdf5read('Clamp1_001.ma','/info/0/values');
    catch
        data2 = [];
        tdata2 = [];
    end
else
    max_file = 1;
end
%capacitance
firstpass_cap = mean(data(1,:));
if ~isempty(data2),
    secondpass_cap = mean(data2(1,:));
    if abs(secondpass_cap-firstpass_cap) > (0.15*firstpass_cap),
        cell_capac = secondpass_cap;
        flag_capac = 1;
    else
        cell_capac = mean([firstpass_cap;secondpass_cap]);
        flag_capac = 0;
    end
else
    cell_capac = firstpass_cap;
    flag_capac = 0;
end
%input resistance
firstpass_ir = mean(data(4,:));
if ~isempty(data2),
    secondpass_ir = mean(data2(4,:));
    if abs(secondpass_ir-firstpass_ir) > (0.15*firstpass_ir),
        cell_inputR = secondpass_ir;
        flag_inputR = 1;
    else
        cell_inputR = mean([firstpass_ir;secondpass_ir]);
        flag_inputR = 0;
    end
else
    cell_inputR = firstpass_ir;
    flag_inputR = 0;
end
%access resistance
firstpass_ar = mean(data(3,:));
if ~isempty(data2),
    secondpass_ar = mean(data2(3,:));
    if abs(secondpass_ar-firstpass_ar) > (0.2*firstpass_ar),
        cell_accessR = secondpass_ar;
        flag_accessR = 1;
    else
        cell_accessR = mean([firstpass_ar;secondpass_ar]);
        flag_accessR = 0;
    end
else
    cell_accessR = firstpass_ar;
    flag_accessR = 0;
end

capacitance.cell = pwd;
capacitance.value = cell_capac;
capacitance.flag = flag_capac;
input_R.cell = pwd;
input_R.value = cell_inputR;
input_R.flag = flag_inputR;
access_R.cell = pwd;
access_R.value = cell_accessR;
access_R.flag = flag_accessR;

cd ..
