% in directory: data/collected_IP/closed 

cellFiles = dir('*.mat');
numCells = length(cellFiles);

for n = 1:numCells,
    s{n} = load(cellFiles(n).name);
end

collated_ISIs = cell(size(s{1,1}.cell_ISIs,1),numCells);
collated_apThreshold = cell(size(s{1,1}.cell_apThreshold,1),numCells);
collated_apMaxslope = cell(size(s{1,1}.cell_apMaxslope,1),numCells);
collated_spiketimes = cell(size(s{1,1}.cell_spiketimes,1),numCells);
collated_peakVm = cell(size(s{1,1}.cell_peakVm,1),numCells);
collated_vmbaseline = cell(size(s{1,1}.cell_vmbaseline,1),numCells);
collated_ap_Ahp = cell(size(s{1,1}.cell_ap_Ahp,1),numCells);

for n = 1:numCells,
    for m = 1:size(s{1,1}.cell_ISIs,1),
        collated_ISIs{m,n} = s{1,n}.cell_ISIs{m,1};
        collated_spiketimes{m,n} = s{1,n}.cell_spiketimes{m,1};
        collated_apThreshold{m,n} = s{1,n}.cell_apThreshold{m,1};
        collated_apMaxslope{m,n} = s{1,n}.cell_apMaxslope{m,1};
        collated_peakVm{m,n} = s{1,n}.cell_peakVm{m,1};
        collated_vmbaseline{m,n} = s{1,n}.cell_vmbaseline{m,1};
        collated_ap_Ahp{m,n} = s{1,n}.cell_ap_Ahp{m,1};
    end
end
emptyCells = cellfun(@isempty,collated_ISIs);
for i = 1:size(emptyCells,1),
    for j = 1:size(emptyCells,2),
        if emptyCells(i,j) == 1,
            collated_ISIs{i,j} = NaN;
        else
        end
    end
end
for a = 1:size(collated_ISIs,1),
    for b = 1:size(collated_ISIs,2),
        c_CV_ISIs(a,b) = std(collated_ISIs{a,b})./nanmean(collated_ISIs{a,b});
    end
end
mc_CV_ISIs = nanmean(c_CV_ISIs,2);
ec_CV_ISIs = nanstd(c_CV_ISIs,1,2)./sqrt(size(c_CV_ISIs,2));
        

for ii = 1:size(collated_ISIs,1),
    for jj = 1:size(collated_ISIs,2),
        m_ISIs{ii,jj} = nanmean(collated_ISIs{ii,jj});
        med_ISIs{ii,jj} = nanmedian(collated_ISIs{ii,jj});
        m_apTH{ii,jj} = nanmean(collated_apThreshold{ii,jj});
        m_peakVm{ii,jj} = nanmean(collated_peakVm{ii,jj});
        m_vmbaseline{ii,jj} = nanmean(collated_vmbaseline{ii,jj});
        m_Ahp{ii,jj} = nanmean(collated_ap_Ahp{ii,jj});
        
    end
end
c_mxs_ISIs = nanmean(cell2mat(m_ISIs),2);
c_Exs_ISIs = std(cell2mat(m_ISIs),1,2)./sqrt(size(m_ISIs,2));
c_m_med_ISIs = nanmean(cell2mat(med_ISIs),2);
c_E_med_ISIs = std(cell2mat(m_ISIs),1,2)./sqrt(size(m_ISIs,2));
c_mxs_apTH = nanmean(cell2mat(m_apTH),2);
c_Exs_apTH = std(cell2mat(m_apTH),1,2)./sqrt(size(m_apTH,2));
c_mxs_peakVm = nanmean(cell2mat(m_peakVm),2);
c_Exs_peakVm = std(cell2mat(m_peakVm),1,2)./sqrt(size(m_peakVm,2));
c_mxs_vmbaseline = nanmean(cell2mat(m_vmbaseline),2);
c_Exs_vmbaseline = std(cell2mat(m_vmbaseline),1,2)./sqrt(size(m_vmbaseline,2));
c_mxs_Ahp = nanmean(cell2mat(m_Ahp),2);
c_Exs_Ahp = std(cell2mat(m_Ahp),1,2)./sqrt(size(m_Ahp,2));

save('closed_processedIP.mat');
clear s
clear collated_ISIs
clear collated_spiketimes
clear collated_apThreshold
clear collated_apMaxslope
clear collated_peakVm
clear collated_vmbaseline
clear colated_ap_Ahp
clear m_ISIs
clear med_ISIs
clear m_apTH
clear m_peakVm
clear m_vmbaseline
clear m_Ahp

cd ../open

cellFiles = dir('*.mat');
numCells = length(cellFiles);

for n = 1:numCells,
    s{n} = load(cellFiles(n).name);
end

collated_ISIs = cell(size(s{1,1}.cell_ISIs,1),numCells);
collated_apThreshold = cell(size(s{1,1}.cell_apThreshold,1),numCells);
collated_apMaxslope = cell(size(s{1,1}.cell_apMaxslope,1),numCells);
collated_spiketimes = cell(size(s{1,1}.cell_spiketimes,1),numCells);
collated_peakVm = cell(size(s{1,1}.cell_peakVm,1),numCells);
collated_vmbaseline = cell(size(s{1,1}.cell_vmbaseline,1),numCells);
collated_ap_Ahp = cell(size(s{1,1}.cell_ap_Ahp,1),numCells);

for n = 1:numCells,
    for m = 1:size(s{1,1}.cell_ISIs,1),
        collated_ISIs{m,n} = s{1,n}.cell_ISIs{m,1};
        collated_spiketimes{m,n} = s{1,n}.cell_spiketimes{m,1};
        collated_apThreshold{m,n} = s{1,n}.cell_apThreshold{m,1};
        collated_apMaxslope{m,n} = s{1,n}.cell_apMaxslope{m,1};
        collated_peakVm{m,n} = s{1,n}.cell_peakVm{m,1};
        collated_vmbaseline{m,n} = s{1,n}.cell_vmbaseline{m,1};
        collated_ap_Ahp{m,n} = s{1,n}.cell_ap_Ahp{m,1};
    end
end
emptyCells = cellfun(@isempty,collated_ISIs);
for i = 1:size(emptyCells,1),
    for j = 1:size(emptyCells,2),
        if emptyCells(i,j) == 1,
            collated_ISIs{i,j} = NaN;
        else
        end
    end
end
for a = 1:size(collated_ISIs,1),
    for b = 1:size(collated_ISIs,2),
        o_CV_ISIs(a,b) = std(collated_ISIs{a,b})./nanmean(collated_ISIs{a,b});
    end
end
mo_CV_ISIs = nanmean(o_CV_ISIs,2);
eo_CV_ISIs = nanstd(o_CV_ISIs,1,2)./sqrt(size(o_CV_ISIs,2));

for ii = 1:size(collated_ISIs,1),
    for jj = 1:size(collated_ISIs,2),
        m_ISIs{ii,jj} = nanmean(collated_ISIs{ii,jj});
        med_ISIs{ii,jj} = nanmedian(collated_ISIs{ii,jj});
        m_apTH{ii,jj} = nanmean(collated_apThreshold{ii,jj});
        m_peakVm{ii,jj} = nanmean(collated_peakVm{ii,jj});
        m_vmbaseline{ii,jj} = nanmean(collated_vmbaseline{ii,jj});
        m_Ahp{ii,jj} = nanmean(collated_ap_Ahp{ii,jj});
        
    end
end
o_mxs_ISIs = nanmean(cell2mat(m_ISIs),2);
o_Exs_ISIs = std(cell2mat(m_ISIs),1,2)./sqrt(size(m_ISIs,2));
o_m_med_ISIs = nanmean(cell2mat(med_ISIs),2);
o_E_med_ISIs = std(cell2mat(m_ISIs),1,2)./sqrt(size(m_ISIs,2));
o_mxs_apTH = nanmean(cell2mat(m_apTH),2);
o_Exs_apTH = std(cell2mat(m_apTH),1,2)./sqrt(size(m_apTH,2));
o_mxs_peakVm = nanmean(cell2mat(m_peakVm),2);
o_Exs_peakVm = std(cell2mat(m_peakVm),1,2)./sqrt(size(m_peakVm,2));
o_mxs_vmbaseline = nanmean(cell2mat(m_vmbaseline),2);
o_Exs_vmbaseline = std(cell2mat(m_vmbaseline),1,2)./sqrt(size(m_vmbaseline,2));
o_mxs_Ahp = nanmean(cell2mat(m_Ahp),2);
o_Exs_Ahp = std(cell2mat(m_Ahp),1,2)./sqrt(size(m_Ahp,2));

save('opened_processedIP.mat');
clear s
clear collated_ISIs
clear collated_spiketimes
clear collated_apThreshold
clear collated_apMaxslope
clear collated_peakVm
clear collated_vmbaseline
clear colated_ap_Ahp
clear m_ISIs
clear med_ISIs
clear m_apTH
clear m_peakVm
clear m_vmbaseline
clear m_Ahp

cd ../DR

cellFiles = dir('*.mat');
numCells = length(cellFiles);

for n = 1:numCells,
    s{n} = load(cellFiles(n).name);
end

collated_ISIs = cell(size(s{1,1}.cell_ISIs,1),numCells);
collated_apThreshold = cell(size(s{1,1}.cell_apThreshold,1),numCells);
collated_apMaxslope = cell(size(s{1,1}.cell_apMaxslope,1),numCells);
collated_spiketimes = cell(size(s{1,1}.cell_spiketimes,1),numCells);
collated_peakVm = cell(size(s{1,1}.cell_peakVm,1),numCells);
collated_vmbaseline = cell(size(s{1,1}.cell_vmbaseline,1),numCells);
collated_ap_Ahp = cell(size(s{1,1}.cell_ap_Ahp,1),numCells);

for n = 1:numCells,
    for m = 1:size(s{1,1}.cell_ISIs,1),
        collated_ISIs{m,n} = s{1,n}.cell_ISIs{m,1};
        collated_spiketimes{m,n} = s{1,n}.cell_spiketimes{m,1};
        collated_apThreshold{m,n} = s{1,n}.cell_apThreshold{m,1};
        collated_apMaxslope{m,n} = s{1,n}.cell_apMaxslope{m,1};
        collated_peakVm{m,n} = s{1,n}.cell_peakVm{m,1};
        collated_vmbaseline{m,n} = s{1,n}.cell_vmbaseline{m,1};
        collated_ap_Ahp{m,n} = s{1,n}.cell_ap_Ahp{m,1};
    end
end
emptyCells = cellfun(@isempty,collated_ISIs);
for i = 1:size(emptyCells,1),
    for j = 1:size(emptyCells,2),
        if emptyCells(i,j) == 1,
            collated_ISIs{i,j} = NaN;
        else
        end
    end
end
for a = 1:size(collated_ISIs,1),
    for b = 1:size(collated_ISIs,2),
        dr_CV_ISIs(a,b) = std(collated_ISIs{a,b})./nanmean(collated_ISIs{a,b});
    end
end
mdr_CV_ISIs = nanmean(dr_CV_ISIs,2);
edr_CV_ISIs = nanstd(dr_CV_ISIs,1,2)./sqrt(size(dr_CV_ISIs,2));

for ii = 1:size(collated_ISIs,1),
    for jj = 1:size(collated_ISIs,2),
        m_ISIs{ii,jj} = nanmean(collated_ISIs{ii,jj});
        med_ISIs{ii,jj} = nanmedian(collated_ISIs{ii,jj});
        m_apTH{ii,jj} = nanmean(collated_apThreshold{ii,jj});
        m_peakVm{ii,jj} = nanmean(collated_peakVm{ii,jj});
        m_vmbaseline{ii,jj} = nanmean(collated_vmbaseline{ii,jj});
        m_Ahp{ii,jj} = nanmean(collated_ap_Ahp{ii,jj});
        
    end
end

dr_mxs_ISIs = nanmean(cell2mat(m_ISIs),2);
dr_Exs_ISIs = std(cell2mat(m_ISIs),1,2)./sqrt(size(m_ISIs,2));
dr_m_med_ISIs = nanmean(cell2mat(med_ISIs),2);
dr_E_med_ISIs = std(cell2mat(m_ISIs),1,2)./sqrt(size(m_ISIs,2));
dr_mxs_apTH = nanmean(cell2mat(m_apTH),2);
dr_Exs_apTH = std(cell2mat(m_apTH),1,2)./sqrt(size(m_apTH,2));
dr_mxs_peakVm = nanmean(cell2mat(m_peakVm),2);
dr_Exs_peakVm = std(cell2mat(m_peakVm),1,2)./sqrt(size(m_peakVm,2));
dr_mxs_vmbaseline = nanmean(cell2mat(m_vmbaseline),2);
dr_Exs_vmbaseline = std(cell2mat(m_vmbaseline),1,2)./sqrt(size(m_vmbaseline,2));
dr_mxs_Ahp = nanmean(cell2mat(m_Ahp),2);
dr_Exs_Ahp = std(cell2mat(m_Ahp),1,2)./sqrt(size(m_Ahp,2));

save('DR_processedIP.mat');
clear s
clear collated_ISIs
clear collated_spiketimes
clear collated_apThreshold
clear collated_apMaxslope
clear collated_peakVm
clear collated_vmbaseline
clear colated_ap_Ahp
clear m_ISIs
clear med_ISIs
clear m_apTH
clear m_peakVm
clear m_vmbaseline
clear m_Ahp

%plots
f = figure;
for i = 1:length(dr_mxs_ISIs),
    subplot(5,3,i);
    errorbar(1.0,c_mxs_ISIs(i,1),c_Exs_ISIs(i,1));
    hold on;
    errorbar(2.0,dr_mxs_ISIs(i,1),dr_Exs_ISIs(i,1));
    errorbar(3.0,o_mxs_ISIs(i,1),o_Exs_ISIs(i,1));
    bar(1.0,c_mxs_ISIs(i,1));
    bar(2.0,dr_mxs_ISIs(i,1));
    bar(3.0,o_mxs_ISIs(i,1));
    local_max = max([c_mxs_ISIs(i,1),dr_mxs_ISIs(i,1),o_mxs_ISIs(i,1)]);
    if local_max>0.12,
        ylim([0 0.1*local_max+local_max]);
    else
        ylim([0 0.12]);
    end
    ylabel('mean ISIs (seconds)');
    title('1 = closed, 2 = dark-reared, 3 = open');
end

f1 = figure;
for i = 1:length(dr_m_med_ISIs),
    subplot(5,3,i);
    errorbar(1.0,c_m_med_ISIs(i,1),c_E_med_ISIs(i,1));
    hold on;
    errorbar(2.0,dr_m_med_ISIs(i,1),dr_E_med_ISIs(i,1));
    errorbar(3.0,o_m_med_ISIs(i,1),o_E_med_ISIs(i,1));
    bar(1.0,c_m_med_ISIs(i,1));
    bar(2.0,dr_m_med_ISIs(i,1));
    bar(3.0,o_m_med_ISIs(i,1));
    local_max = max([c_m_med_ISIs(i,1),dr_m_med_ISIs(i,1),o_m_med_ISIs(i,1)]);
    if local_max>0.12,
        ylim([0 0.1*local_max+local_max]);
    else
        ylim([0 0.12]);
    end
    ylabel('mean median-ISIs (sec)');
    title('1 = closed, 2 = dark-reared, 3 = open');
end

f2 = figure;
for i = 1:length(dr_mxs_apTH),
    subplot(5,3,i);
    errorbar(1.0,c_mxs_apTH(i,1),c_Exs_apTH(i,1));
    hold on;
    errorbar(2.0,dr_mxs_apTH(i,1),dr_Exs_apTH(i,1));
    errorbar(3.0,o_mxs_apTH(i,1),o_Exs_apTH(i,1));
    bar(1.0,c_mxs_apTH(i,1));
    bar(2.0,dr_mxs_apTH(i,1));
    bar(3.0,o_mxs_apTH(i,1));
    local_max = max([c_mxs_apTH(i,1),dr_mxs_apTH(i,1),o_mxs_apTH(i,1)]);
    if local_max>0.12,
        ylim([0 0.1*local_max+local_max]);
    else
        ylim([-0.06 0.01]);
    end
    ylabel('spike threshold (V)');
    title('1 = closed, 2 = dark-reared, 3 = open');
end

f3 = figure;
for i = 1:length(dr_mxs_peakVm),
    subplot(5,3,i);
    errorbar(1.0,c_mxs_peakVm(i,1),c_Exs_peakVm(i,1));
    hold on;
    errorbar(2.0,dr_mxs_peakVm(i,1),dr_Exs_peakVm(i,1));
    errorbar(3.0,o_mxs_peakVm(i,1),o_Exs_peakVm(i,1));
    bar(1.0,c_mxs_peakVm(i,1));
    bar(2.0,dr_mxs_peakVm(i,1));
    bar(3.0,o_mxs_peakVm(i,1));
    local_max = max([c_mxs_peakVm(i,1),dr_mxs_peakVm(i,1),o_mxs_peakVm(i,1)]);
    if local_max>0.08,
        ylim([0 0.1*local_max+local_max]);
    else
        ylim([0 0.08]);
    end
    ylabel('mean peak Vm (V)');
    title('1 = closed, 2 = dark-reared, 3 = open');
end

f4 = figure;
for i = 1:length(dr_mxs_vmbaseline),
    subplot(5,3,i);
    errorbar(1.0,c_mxs_vmbaseline(i,1),c_Exs_vmbaseline(i,1));
    hold on;
    errorbar(2.0,dr_mxs_vmbaseline(i,1),dr_Exs_vmbaseline(i,1));
    errorbar(3.0,o_mxs_vmbaseline(i,1),o_Exs_vmbaseline(i,1));
    bar(1.0,c_mxs_vmbaseline(i,1));
    bar(2.0,dr_mxs_vmbaseline(i,1));
    bar(3.0,o_mxs_vmbaseline(i,1));
    local_min = min([c_mxs_vmbaseline(i,1),dr_mxs_vmbaseline(i,1),o_mxs_vmbaseline(i,1)]);
    if local_min<-0.09,
        ylim([0 0.1*local_max+local_max]);
    else
        ylim([-0.09 0]);
    end
    ylabel('mean Vm baseline (V)');
    title('1 = closed, 2 = dark-reared, 3 = open');
end

f5 = figure;
for i = 1:length(dr_mxs_Ahp),
    subplot(5,3,i);
    errorbar(1.0,c_mxs_Ahp(i,1),c_Exs_Ahp(i,1));
    hold on;
    errorbar(2.0,dr_mxs_Ahp(i,1),dr_Exs_Ahp(i,1));
    errorbar(3.0,o_mxs_Ahp(i,1),o_Exs_Ahp(i,1));
    bar(1.0,c_mxs_Ahp(i,1));
    bar(2.0,dr_mxs_Ahp(i,1));
    bar(3.0,o_mxs_Ahp(i,1));
    local_min = min([c_mxs_Ahp(i,1),dr_mxs_Ahp(i,1),o_mxs_Ahp(i,1)]);
    if local_min<-0.07,
        ylim([0 0.1*local_max+local_max]);
    else
        ylim([-0.07 -0.06]);
    end
    ylabel('mean AHP trough (V)');
    title('1 = closed, 2 = dark-reared, 3 = open');
end

f6 = figure;
for i = 1:length(mdr_CV_ISIs),
    subplot(5,3,i);
    errorbar(1.0,mc_CV_ISIs(i,1),ec_CV_ISIs(i,1));
    hold on;
    errorbar(2.0,mdr_CV_ISIs(i,1),edr_CV_ISIs(i,1));
    errorbar(3.0,mo_CV_ISIs(i,1),eo_CV_ISIs(i,1));
    bar(1.0,mc_CV_ISIs(i,1));
    bar(2.0,mdr_CV_ISIs(i,1));
    bar(3.0,mo_CV_ISIs(i,1));
    local_max = max([mc_CV_ISIs(i,1),mdr_CV_ISIs(i,1),mo_CV_ISIs(i,1)]);
    if local_max>1.0,
        ylim([0 0.1*local_max+local_max]);
    else
        ylim([0 1.0]);
    end
    ylabel('mean CV of ISIs');
    title('1 = closed, 2 = dark-reared, 3 = open');
end
