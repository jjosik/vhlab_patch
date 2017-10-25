
%NOTE - if using function to analyze rise/decay times, need to "uncomment"
%the relevant variable names in the file loop - these elements were
%introduced later to use the same code to analyze the "risedecay"
%directory.  Leave unchanged to run on "IEI" directory.

cell_results = uigetdir();
files = dir(fullfile(cell_results,'*.mat'));
use_age = 0;  %recommend setting as function input

collated_E_IEI = [];
collated_I_IEI = [];
collated_alt_E_IEI = [];
collated_alt_I_IEI = [];
collated_E_fmean = [];
collated_I_fmean = [];
collated_alt_E_fmean = [];
collated_alt_I_fmean = [];
collated_E_deconv_amps = [];
collated_E_alt_amps = [];
collated_I_deconv_amps = [];
collated_I_alt_amps = [];
collated_E_th_durs = [];
collated_E_fwhm = [];
collated_I_th_durs = [];
collated_I_fwhm = [];
collated_E_risetime = [];
collated_I_risetime = [];
collated_E_decaytau = [];
collated_I_decaytau = [];
cell_E_alt_amps = cell(length(files),1);
cell_I_alt_amps = cell(length(files),1);
cell_E_th_durs = cell(length(files),1);
cell_I_th_durs = cell(length(files),1);
cell_E_fwhm = cell(length(files),1);
cell_I_fwhm = cell(length(files),1);
cell_E_risetime = cell(length(files),1);
cell_I_risetime = cell(length(files),1);
cell_E_decaytau = cell(length(files),1);
cell_I_decaytau = cell(length(files),1);
matched_order_cell_experience = [];
matched_order_age = [];
matched_order_f_experience_E = [];
matched_order_f_experience_I = [];
matched_order_d_experience_E = [];
matched_order_d_experience_I = [];
matched_order_a_experience_E = [];
matched_order_a_experience_I = [];
f_setnum_E = cell(length(files),1);
f_setnum_I = cell(length(files),1);
%d_setnum = cell(length(files),1);
%a_setnum = cell(length(files),1);

for i = 1:length(files),
    if files(i).isdir,
        continue
    else
    end
    cell_frequencydata = load(files(i).name,'cell_frequencydata');
    cell_amplitudedata = load(files(i).name,'cell_amplitudedata');
    cell_durationdata = load(files(i).name,'cell_durationdata');
    for j = 1:size(cell_frequencydata.cell_frequencydata,1),
        %frequency data
        collated_E_IEI = [collated_E_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI];
        collated_I_IEI = [collated_I_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI];
        collated_alt_E_IEI = [collated_alt_E_IEI;cell_frequencydata.cell_frequencydata{j,1}.alt_cell_total_E_IEI];
        collated_alt_I_IEI = [collated_alt_I_IEI;cell_frequencydata.cell_frequencydata{j,1}.alt_cell_total_I_IEI];
        collated_E_fmean = cat(2,collated_E_fmean,cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean);
        collated_alt_E_fmean = cat(2,collated_alt_E_fmean,cell_frequencydata.cell_frequencydata{j,1}.alt_cell_total_E_fmean);
        collated_I_fmean = cat(2,collated_I_fmean,cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean);
        collated_alt_I_fmean = cat(2,collated_alt_I_fmean,cell_frequencydata.cell_frequencydata{j,1}.alt_cell_total_I_fmean);
        matched_order_age = [matched_order_age,repmat(cell_frequencydata.cell_frequencydata{j,1}.cell_ages,1,...
            length(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean))];
        matched_order_f_experience_E = [matched_order_f_experience_E,repmat(cell_frequencydata.cell_frequencydata{j,1}.cell_experience,...
            1,length(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean))];
        matched_order_f_experience_I = [matched_order_f_experience_I,repmat(cell_frequencydata.cell_frequencydata{j,1}.cell_experience,...
            1,length(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean))];
        f_setnum_E{i,1} = length(collated_alt_E_IEI);
        f_setnum_I{i,1} = length(collated_alt_I_IEI);
        f_Esets = length(cell_frequencydata.cell_frequencydata{j,1}.alt_cell_total_E_fmean);
        f_Isets = length(cell_frequencydata.cell_frequencydata{j,1}.alt_cell_total_I_fmean);
        %amplitude data
        cell_E_alt_amps{i,1} = cell_amplitudedata.cell_amplitudedata{j,1}.E_trial_alt_amps;
        cell_I_alt_amps{i,1} = cell_amplitudedata.cell_amplitudedata{j,1}.I_trial_alt_amps;
        collated_E_deconv_amps = [collated_E_deconv_amps;cell_amplitudedata.cell_amplitudedata{j,1}.E_trial_deconv_amps];
        collated_E_alt_amps = [collated_E_alt_amps;cell_amplitudedata.cell_amplitudedata{j,1}.E_trial_alt_amps];
        collated_I_deconv_amps = [collated_E_deconv_amps;cell_amplitudedata.cell_amplitudedata{j,1}.I_trial_deconv_amps];
        collated_I_alt_amps = [collated_I_alt_amps;cell_amplitudedata.cell_amplitudedata{j,1}.I_trial_alt_amps];
        a_setnum_E{i,1} = length(collated_E_alt_amps);
        a_setnum_I{i,1} = length(collated_I_alt_amps);
        matched_order_a_experience_E = [matched_order_a_experience_E,repmat(cell_amplitudedata.cell_amplitudedata{j,1}.cell_experience,...
            1,length(cell_amplitudedata.cell_amplitudedata{j,1}.E_trial_alt_amps))];
        matched_order_a_experience_I = [matched_order_a_experience_I,repmat(cell_amplitudedata.cell_amplitudedata{j,1}.cell_experience,...
            1,length(cell_amplitudedata.cell_amplitudedata{j,1}.I_trial_alt_amps))];
        cell_names{i,1} = cell_amplitudedata.cell_amplitudedata{1,1}(:);
        %duration data
        collated_E_th_durs = [collated_E_th_durs;cell_durationdata.cell_durationdata{j,1}.E_trial_th_durations];
        collated_E_fwhm = [collated_E_fwhm;cell_durationdata.cell_durationdata{j,1}.E_trial_fwhm];
        collated_I_th_durs = [collated_I_th_durs;cell_durationdata.cell_durationdata{j,1}.I_trial_th_durations];
        collated_I_fwhm = [collated_I_fwhm;cell_durationdata.cell_durationdata{j,1}.I_trial_fwhm];
        cell_E_th_durs{i,1} = cell_durationdata.cell_durationdata{j,1}.E_trial_th_durations;
        cell_I_th_durs{i,1} = cell_durationdata.cell_durationdata{j,1}.I_trial_th_durations;
        cell_E_fwhm{i,1} = cell_durationdata.cell_durationdata{j,1}.E_trial_fwhm;
        cell_I_fwhm{i,1} = cell_durationdata.cell_durationdata{j,1}.I_trial_fwhm;
        %collated_E_risetime = [collated_E_risetime;cell_durationdata.cell_durationdata{j,1}.E_trial_risetime];
        %collated_I_risetime = [collated_I_risetime;cell_durationdata.cell_durationdata{j,1}.I_trial_risetime];
        %collated_E_decaytau = [collated_E_decaytau;cell_durationdata.cell_durationdata{j,1}.E_trial_decaytau];
        %collated_I_decaytau = [collated_I_decaytau;cell_durationdata.cell_durationdata{j,1}.I_trial_decaytau];
        %cell_E_risetime{i,1} = cell_durationdata.cell_durationdata{j,1}.E_trial_risetime;
        %cell_I_risetime{i,1} = cell_durationdata.cell_durationdata{j,1}.I_trial_risetime;
        %cell_E_decaytau{i,1} = cell_durationdata.cell_durationdata{j,1}.E_trial_decaytau;
        %cell_I_decaytau{i,1} = cell_durationdata.cell_durationdata{j,1}.I_trial_decaytau;
        matched_order_cell_experience = [matched_order_cell_experience;cell_durationdata.cell_durationdata{j,1}.cell_experience];
        matched_order_d_experience_E = [matched_order_d_experience_E,repmat(cell_durationdata.cell_durationdata{j,1}.cell_experience,...
            1,length(cell_durationdata.cell_durationdata{j,1}.E_trial_fwhm))];
        matched_order_d_experience_I = [matched_order_d_experience_I,repmat(cell_durationdata.cell_durationdata{j,1}.cell_experience,...
            1,length(cell_durationdata.cell_durationdata{j,1}.I_trial_fwhm))];
        %d_setnum = cumsum(d_setnum)+;
    end
    
end  

%matched_order_experience = matched_order_experience(1:end-1);
%*************
%****PLOTS****
%*************
%scatter plots of binned mean event rates
bin_sizes = [0.01,0.1,1.0,10];
for z = 1:size(cell_frequencydata.cell_frequencydata{1,1}.alt_cell_total_E_fmean,1),
    f = figure;
    scatter((matched_order_f_experience_E)+0.1,collated_alt_E_fmean(z,:),'b');
    hold on;
    scatter((matched_order_f_experience_I)-0.1,collated_alt_I_fmean(z,:),'r');
    xlabel('Days of experience');
    ylabel(['Mean frequency ', num2str(bin_sizes(1,z)), ' sec. bins']);
    title(['Trial mean PSC events per second vs. experience,',num2str(bin_sizes(1,z)),' sec. bins']);
    hold off;
end
%scatterplots of trial mean event durations
%NOTE - collated_E_fwhm is collected individual event values, need to sort
%into trial means
f = figure;
scatter((matched_order_d_experience_E)+0.1,collated_E_fwhm,'b');
hold on;
scatter((matched_order_d_experience_I)-0.1,collated_I_fwhm,'r');
xlabel('Days of experience');
ylabel('Mean FWHM');
title('Trial mean event duration (FWHM)');
hold off;

%Rise/decay times
pre_E_rise = [];
pre_I_rise = [];
pre_E_decay = [];
pre_I_decay = [];
dur_E_rise = [];
dur_I_rise = [];
dur_E_decay = [];
dur_I_decay = [];
post_E_rise = [];
post_I_rise = [];
post_E_decay = [];
post_I_decay = [];
for i = 1:length(matched_order_cell_experience)
    if matched_order_cell_experience(i,1) <= 0,
        pre_E_rise(end+1,1) = nanmedian(cell_E_risetime{i,1});
        pre_I_rise(end+1,1) = nanmedian(cell_I_risetime{i,1});
        pre_E_decay(end+1,1) = nanmedian(cell_E_decaytau{i,1});
        pre_I_decay(end+1,1) = nanmedian(cell_I_decaytau{i,1});
    elseif (matched_order_cell_experience(i,1) > 0)&&(matched_order_cell_experience(i,1) < 5),
        dur_E_rise(end+1,1) = nanmedian(cell_E_risetime{i,1});
        dur_I_rise(end+1,1) = nanmedian(cell_I_risetime{i,1});
        dur_E_decay(end+1,1) = nanmedian(cell_E_decaytau{i,1});
        dur_I_decay(end+1,1) = nanmedian(cell_I_decaytau{i,1});
    elseif matched_order_cell_experience(i,1) >= 5,
        post_E_rise(end+1,1) = nanmedian(cell_E_risetime{i,1});
        post_I_rise(end+1,1) = nanmedian(cell_I_risetime{i,1});
        post_E_decay(end+1,1) = nanmedian(cell_E_decaytau{i,1});
        post_I_decay(end+1,1) = nanmedian(cell_I_decaytau{i,1});
    end
end
f = figure;
boxwhisker(pre_E_rise.*1000,0.75,[0 0 1]);
hold on;
boxwhisker(pre_I_rise.*1000,1.25,[1 0 0]);
boxwhisker(dur_E_rise.*1000,2.75,[0 0 1]);
boxwhisker(dur_I_rise.*1000,3.25,[1 0 0]);
boxwhisker(post_E_rise.*1000,4.75,[0 0 1]);
boxwhisker(post_I_rise.*1000,5.25,[1 0 0]);
set(gca,'XTick',[1,3,5],'XTickLabels',{'Pre-opening','0-5 days','>5 days'});
xlabel('Experience');
ylabel('Rise time (ms)');
xlim([0 6]);
scatter(0.75.*ones(length(pre_E_rise),1),pre_E_rise.*1000,'k');
scatter(1.25.*ones(length(pre_I_rise),1),pre_I_rise.*1000,'k');
scatter(2.75.*ones(length(dur_E_rise),1),dur_E_rise.*1000,'k');
scatter(3.25.*ones(length(dur_I_rise),1),dur_I_rise.*1000,'k');
scatter(4.75.*ones(length(post_E_rise),1),post_E_rise.*1000,'k');
scatter(5.25.*ones(length(post_I_rise),1),post_I_rise.*1000,'k');
hold off;
%kruskal-wallis omnibus test
dataE = [pre_E_rise;dur_E_rise;post_E_rise];
groupsE = [ones(length(pre_E_rise),1);2.*ones(length(dur_E_rise),1);3.*ones(length(post_E_rise),1)];
dataI = [pre_I_rise;dur_I_rise;post_I_rise];
groupsI = [ones(length(pre_I_rise),1);2.*ones(length(dur_I_rise),1);3.*ones(length(post_I_rise),1)];
[p_E,anovatabE,stats_E] = kruskalwallis(dataE,groupsE);
[p_I,anovatabI,stats_I] = kruskalwallis(dataI,groupsI);
%B-c Mann-Whitney
alpha = 0.05;
group_N = 3;
corrected_alpha = alpha/((group_N*(group_N-1))/2);
p_E_early =

f2 = figure;
boxwhisker(pre_E_decay.*1000,0.75,[0 0 1]);
hold on;
boxwhisker(pre_I_decay.*1000,1.25,[1 0 0]);
boxwhisker(dur_E_decay.*1000,2.75,[0 0 1]);
boxwhisker(dur_I_decay.*1000,3.25,[1 0 0]);
boxwhisker(post_E_decay.*1000,4.75,[0 0 1]);
boxwhisker(post_I_decay.*1000,5.25,[1 0 0]);
set(gca,'XTick',[1,3,5],'XTickLabels',{'Pre-opening','0-5 days','>5 days'});
xlabel('Experience');
ylabel('Decay tau (ms)');
xlim([0 6]);
scatter(0.75.*ones(length(pre_E_decay),1),pre_E_decay.*1000,'k');
scatter(1.25.*ones(length(pre_I_decay),1),pre_I_decay.*1000,'k');
scatter(2.75.*ones(length(dur_E_decay),1),dur_E_decay.*1000,'k');
scatter(3.25.*ones(length(dur_I_decay),1),dur_I_decay.*1000,'k');
scatter(4.75.*ones(length(post_E_decay),1),post_E_decay.*1000,'k');
scatter(5.25.*ones(length(post_I_decay),1),post_I_decay.*1000,'k');
%K_W omnibus
ddataE = [pre_E_decay;dur_E_decay;post_E_decay];
dgroupsE = [ones(length(pre_E_decay),1);2.*ones(length(dur_E_decay),1);3.*ones(length(post_E_decay),1)];
ddataI = [pre_I_decay;dur_I_decay;post_I_decay];
dgroupsI = [ones(length(pre_I_decay),1);2.*ones(length(dur_I_decay),1);3.*ones(length(post_I_decay),1)];
[dp_E,anovatabdE,stats_dE] = kruskalwallis(ddataE,dgroupsE);
[dp_I,anovatabdI,stats_dI] = kruskalwallis(ddataI,dgroupsI);
%B-c Mann-Whitney
alpha = 0.05;
groups = 3;
corrected_alpha = alpha/((groups*(groups-1))/2);
p_E_early = ranksum(pre_E_decay,dur_E_decay);
p_E_late = ranksum(dur_E_decay,post_E_decay);
p_E_distal = ranksum(pre_E_decay,post_E_decay);
p_I_early = ranksum(pre_I_decay,dur_I_decay);
p_I_late = ranksum(dur_I_decay,post_I_decay);
p_I_distal = ranksum(pre_I_decay,post_I_decay);

%Duration (fwhm) cdfs
%POOLED
f = figure;
plot_step = [0:0.0001:0.2];
bin_centers_d = (plot_step(1:end-1)+plot_step(2:end))/2;
pre_E_fwhm = [];
dur_E_fwhm = [];
post_E_fwhm = [];
pre_I_fwhm = [];
dur_I_fwhm = [];
post_I_fwhm = [];
for i = 1:length(matched_order_d_experience_E),
    if matched_order_d_experience_E(1,i) <= 0,
        pre_E_fwhm(end+1,1) = collated_E_fwhm(i,1);
    elseif (matched_order_d_experience_E(1,i) > 0)&&(matched_order_d_experience_E(1,i) < 5),
        dur_E_fwhm(end+1,1) = collated_E_fwhm(i,1);
    elseif matched_order_d_experience_E(1,i) > 5,
        post_E_fwhm(end+1,1) = collated_E_fwhm(i,1);
    else
    end
end
for j = 1:length(matched_order_d_experience_I),
    if matched_order_d_experience_I(1,j) <= 0,
        pre_I_fwhm(end+1,1) = collated_I_fwhm(j,1);
    elseif (matched_order_d_experience_I(1,j) > 0)&&(matched_order_d_experience_I(1,j) < 5),
        dur_I_fwhm(end+1,1) = collated_I_fwhm(j,1);
    elseif matched_order_d_experience_I(1,j) > 5,
        post_I_fwhm(end+1,1) = collated_I_fwhm(j,1);
    else
    end
end
pre_E_dur = cumsum(histc(pre_E_fwhm,bin_centers_d));
dur_E_dur = cumsum(histc(dur_E_fwhm,bin_centers_d));
post_E_dur = cumsum(histc(post_E_fwhm,bin_centers_d));
pre_I_dur = cumsum(histc(pre_I_fwhm,bin_centers_d));
dur_I_dur = cumsum(histc(dur_I_fwhm,bin_centers_d));
post_I_dur = cumsum(histc(post_I_fwhm,bin_centers_d));
subplot(1,2,1);
plot(bin_centers_d,(pre_E_dur/max(pre_E_dur)),'b--');
hold on;
plot(bin_centers_d,(dur_E_dur/max(dur_E_dur)),'b-');
plot(bin_centers_d,(post_E_dur/max(post_E_dur)),'b:');
xlabel('Spontaneous event FWHM durations (sec.)');
ylabel('Proportion');
if use_age == 0,
    legend('pre-opening','0-5 days open','>5 days open','Location','East');
else
    legend('<p30','p30-35','>p35','Location','East');
end
hold off
subplot(1,2,2);
plot(bin_centers_d,(pre_I_dur/max(pre_I_dur)),'r--');
hold on;
plot(bin_centers_d,(dur_I_dur/max(dur_I_dur)),'r-');
plot(bin_centers_d,(post_I_dur/max(post_I_dur)),'r:');
xlabel('Spontaneous event FWHM durations (sec.)');
ylabel('Proportion');
if use_age == 0,
    legend('pre-opening','0-5 days open','>5 days open','Location','East');
else
    legend('<p30','p30-35','>p35','Location','East');
end
hold off;
%BY CELL
cpre_E_fwhm = [];
cpre_I_fwhm = [];
cpre_E_th_durs = [];
cpre_I_th_durs = [];
cdur_E_fwhm = [];
cdur_I_fwhm = [];
cdur_E_th_durs = [];
cdur_I_th_durs = [];
cpost_E_fwhm = [];
cpost_I_fwhm = [];
cpost_E_th_durs = [];
cpost_I_th_durs = [];
for i = 1:length(matched_order_cell_experience)
    if matched_order_cell_experience(i,1) <= 0,
        cpre_E_fwhm(end+1,1) = nanmedian(cell_E_fwhm{i,1});
        cpre_I_fwhm(end+1,1) = nanmedian(cell_I_fwhm{i,1});
        cpre_E_th_durs(end+1,1) = nanmedian(cell_E_th_durs{i,1});
        cpre_I_th_durs(end+1,1) = nanmedian(cell_I_th_durs{i,1});
    elseif (matched_order_cell_experience(i,1) > 0)&&(matched_order_cell_experience(i,1) < 5),
        cdur_E_fwhm(end+1,1) = nanmedian(cell_E_fwhm{i,1});
        cdur_I_fwhm(end+1,1) = nanmedian(cell_I_fwhm{i,1});
        cdur_E_th_durs(end+1,1) = nanmedian(cell_E_th_durs{i,1});
        cdur_I_th_durs(end+1,1) = nanmedian(cell_I_th_durs{i,1});
    elseif matched_order_cell_experience(i,1) >= 5,
        cpost_E_fwhm(end+1,1) = nanmedian(cell_E_fwhm{i,1});
        cpost_I_fwhm(end+1,1) = nanmedian(cell_I_fwhm{i,1});
        cpost_E_th_durs(end+1,1) = nanmedian(cell_E_th_durs{i,1});
        cpost_I_th_durs(end+1,1) = nanmedian(cell_I_th_durs{i,1});
    end
end
f = figure;
boxwhisker(cpre_E_fwhm.*1000,0.75,[0 0 1]);
hold on;
boxwhisker(cpre_I_fwhm.*1000,1.25,[1 0 0]);
boxwhisker(cdur_E_fwhm.*1000,2.75,[0 0 1]);
boxwhisker(cdur_I_fwhm.*1000,3.25,[1 0 0]);
boxwhisker(cpost_E_fwhm.*1000,4.75,[0 0 1]);
boxwhisker(cpost_I_fwhm.*1000,5.25,[1 0 0]);
xlabel('Experience');
set(gca,'XTick',[1,3,5],'XTickLabels',{'Pre-opening','0-5 days','>5 days'});
xlim([0 6]);
ylabel('FWHM duration (ms)');
scatter(0.75.*ones(length(cpre_E_fwhm),1),cpre_E_fwhm.*1000,'k');
scatter(1.25.*ones(length(cpre_E_fwhm),1),cpre_I_fwhm.*1000,'k');
scatter(2.75.*ones(length(cdur_E_fwhm),1),cdur_E_fwhm.*1000,'k');
scatter(3.25.*ones(length(cdur_I_fwhm),1),cdur_I_fwhm.*1000,'k');
scatter(4.75.*ones(length(cpost_E_fwhm),1),cpost_E_fwhm.*1000,'k');
scatter(5.25.*ones(length(cpost_I_fwhm),1),cpost_I_fwhm.*1000,'k');
%K-W omnibus
fdataE = [cpre_E_fwhm;cdur_E_fwhm;cpost_E_fwhm];
fgroupsE = [ones(length(cpre_E_fwhm),1);2.*ones(length(cdur_E_fwhm),1);3.*ones(length(cpost_E_fwhm),1)];
fdataI = [cpre_I_fwhm;cdur_I_fwhm;cpost_I_fwhm];
fgroupsI = [ones(length(cpre_I_fwhm),1);2.*ones(length(cdur_I_fwhm),1);3.*ones(length(cpost_I_fwhm),1)];
[fp_E,anovatabfE,stats_fE] = kruskalwallis(fdataE,fgroupsE);
[fp_I,anovatabfI,stats_fI] = kruskalwallis(fdataI,fgroupsI);
%B-c Mann-Whitney post hoc
alpha = 0.05;
groups = 3;
corrected_alpha = alpha/((groups*(groups-1))/2);
fp_E_early = ranksum(cpre_E_fwhm,cdur_E_fwhm);
fp_E_late = ranksum(cdur_E_fwhm,cpost_E_fwhm);
fp_E_distal = ranksum(cpre_E_fwhm,cpost_E_fwhm);
fp_I_early = ranksum(cpre_I_fwhm,cdur_I_fwhm);
fp_I_late = ranksum(cdur_I_fwhm,cpost_I_fwhm);
fp_I_distal = ranksum(cpre_I_fwhm,post_I_fwhm);
%TH-to-TH
f2 = figure;
boxwhisker(cpre_E_th_durs.*1000,0.75,[0 0 1]);
hold on;
boxwhisker(cpre_I_th_durs.*1000,1.25,[1 0 0]);
boxwhisker(cdur_E_th_durs.*1000,2.75,[0 0 1]);
boxwhisker(cdur_I_th_durs.*1000,3.25,[1 0 0]);
boxwhisker(cpost_E_th_durs.*1000,4.75,[0 0 1]);
boxwhisker(cpost_I_th_durs.*1000,5.25,[1 0 0]);
xlabel('Experience');
set(gca,'XTick',[1,3,5],'XTickLabels',{'Pre-opening','0-5 days','>5 days'});
xlim([0 6]);
ylabel('Theshold-to-threshold duration (ms)');
scatter(0.75.*ones(length(cpre_E_th_durs),1),cpre_E_th_durs.*1000,'k');
scatter(1.25.*ones(length(cpre_E_th_durs),1),cpre_I_th_durs.*1000,'k');
scatter(2.75.*ones(length(cdur_E_th_durs),1),cdur_E_th_durs.*1000,'k');
scatter(3.25.*ones(length(cdur_I_th_durs),1),cdur_I_th_durs.*1000,'k');
scatter(4.75.*ones(length(cpost_E_th_durs),1),cpost_E_th_durs.*1000,'k');
scatter(5.25.*ones(length(cpost_I_th_durs),1),cpost_I_th_durs.*1000,'k');
%K-W omnibus
tdataE = [cpre_E_th_durs;cdur_E_th_durs;cpost_E_th_durs];
tgroupsE = [ones(length(cpre_E_th_durs),1);2.*ones(length(cdur_E_th_durs),1);3.*ones(length(cpost_E_th_durs),1)];
tdataI = [cpre_I_th_durs;cdur_I_th_durs;cpost_I_th_durs];
tgroupsI = [ones(length(cpre_I_th_durs),1);2.*ones(length(cdur_I_th_durs),1);3.*ones(length(cpost_I_th_durs),1)];
[tp_E,anovatabtE,stats_tE] = kruskalwallis(tdataE,tgroupsE);
[tp_I,anovatabtI,stats_tI] = kruskalwallis(tdataI,tgroupsI);
%B-c Mann-Whitney post hoc
alpha = 0.05;
groups = 3;
corrected_alpha = alpha/((groups*(groups-1))/2);
tp_E_early = ranksum(cpre_E_th_durs,cdur_E_th_durs);
tp_E_late = ranksum(cdur_E_th_durs,cpost_E_th_durs);
tp_E_distal = ranksum(cpre_E_th_durs,cpost_E_th_durs);
tp_I_early = ranksum(cpre_I_th_durs,cdur_I_th_durs);
tp_I_late = ranksum(cdur_I_th_durs,cpost_I_th_durs);
tp_I_distal = ranksum(cpre_I_th_durs,cpost_I_th_durs);


%Frequency (cdfs)
f_setnum_E = cell2mat(f_setnum_E);
f_setnum_I = cell2mat(f_setnum_I);
f = figure;
plot_step = [0:0.001:0.5];
E_matched_cell_exp_IEI = matched_order_f_experience_E(1,1:f_Esets:end);
I_matched_cell_exp_IEI = matched_order_f_experience_I(1,1:f_Isets:end);
bin_centers_f = (plot_step(1:end-1)+plot_step(2:end))/2;
pre_E_freq = cell(length(E_matched_cell_exp_IEI),1);
dur_E_freq = cell(length(E_matched_cell_exp_IEI),1);
post_E_freq = cell(length(E_matched_cell_exp_IEI),1);
pre_I_freq = cell(length(I_matched_cell_exp_IEI),1);
dur_I_freq = cell(length(I_matched_cell_exp_IEI),1);
post_I_freq = cell(length(I_matched_cell_exp_IEI),1);
initial_setE = 0;
initial_setI = 0;
for k = 1:length(E_matched_cell_exp_IEI),
    if E_matched_cell_exp_IEI(1,k) <= 0,
        pre_E_freq{k,1} = collated_alt_E_IEI((initial_setE+1):f_setnum_E(k,1),1);
    elseif (E_matched_cell_exp_IEI(1,k) > 0)&&(E_matched_cell_exp_IEI(1,k) < 5),
        dur_E_freq{k,1} = collated_alt_E_IEI((initial_setE+1):f_setnum_E(k,1),1);
    elseif E_matched_cell_exp_IEI(1,k) > 5,
        post_E_freq{k,1} = collated_alt_E_IEI((initial_setE+1):f_setnum_E(k,1),1);
    else
    end
    initial_setE = f_setnum_E(k,1);
end
for kk = 1:length(I_matched_cell_exp_IEI),
    if I_matched_cell_exp_IEI(1,kk) <= 0,
        pre_I_freq{kk,1} = collated_alt_I_IEI((initial_setI+1):f_setnum_I(kk,1),1);
    elseif (I_matched_cell_exp_IEI(1,kk) > 0)&&(I_matched_cell_exp_IEI(1,kk) < 5),
        dur_I_freq{kk,1} = collated_alt_I_IEI((initial_setI+1):f_setnum_I(kk,1),1);
    elseif I_matched_cell_exp_IEI(1,kk) > 5,
        post_I_freq{kk,1} = collated_alt_I_IEI((initial_setI+1):f_setnum_I(kk,1),1);
    else
    end
    initial_setI = f_setnum_I(kk,1);
end
pre_E_freq = cell2mat(pre_E_freq(~(cellfun('isempty',pre_E_freq))));
dur_E_freq = cell2mat(dur_E_freq(~(cellfun('isempty',dur_E_freq))));
post_E_freq = cell2mat(post_E_freq(~(cellfun('isempty',post_E_freq))));
pre_I_freq = cell2mat(pre_I_freq(~(cellfun('isempty',pre_I_freq))));
dur_I_freq = cell2mat(dur_I_freq(~(cellfun('isempty',dur_I_freq))));
post_I_freq = cell2mat(post_I_freq(~(cellfun('isempty',post_I_freq))));
pre_E = cumsum(histc(pre_E_freq,bin_centers_f));
dur_E = cumsum(histc(dur_E_freq,bin_centers_f));
post_E = cumsum(histc(post_E_freq,bin_centers_f));
pre_I = cumsum(histc(pre_I_freq,bin_centers_f));
dur_I = cumsum(histc(dur_I_freq,bin_centers_f));
post_I = cumsum(histc(post_I_freq,bin_centers_f));
subplot(1,2,1);
plot(bin_centers_f,(pre_E/max(pre_E)),'b--');
hold on;
plot(bin_centers_f,(dur_E/max(dur_E)),'b-');
plot(bin_centers_f,(post_E/max(post_E)),'b:');
xlabel('Excitatory IEIs');
ylabel('Proportion');
if use_age == 0,
    legend('pre-opening','0-5 days open','>5 days open','Location','East');
else
    legend('<p30','p30-35','>p35','Location','East');
end
hold off
subplot(1,2,2);
plot(bin_centers_f,(pre_I/max(pre_I)),'r--');
hold on;
plot(bin_centers_f,(dur_I/max(dur_I)),'r-');
plot(bin_centers_f,(post_I/max(post_I)),'r:');
xlabel('Inhibitory IEIs');
ylabel('Proportion');
if use_age == 0,
    legend('pre-opening','0-5 days open','>5 days open','Location','East');
else
    legend('<p30','p30-35','>p35','Location','East');
end
hold off;
%Alt. barplots:
%CASE 1 -- how NOT TO DO IT: means of medians
figure;
bar(0.75,mean(mean(PRE_med_E_IEI,2)),0.3);
hold on;
errorbar(0.75,mean(mean(PRE_med_E_IEI,2)),std(mean(PRE_med_E_IEI,2))/size(PRE_med_E_IEI,1));
boxwhisker(PRE_med_E_IEI,0.75);
bar(1.25,mean(mean(PRE_med_I_IEI,2)),0.3);
errorbar(1.25,mean(mean(PRE_med_I_IEI,2)),std(mean(PRE_med_I_IEI,2))/size(PRE_med_I_IEI,1));
boxwhisker(PRE_med_I_IEI,0.75);
bar(2.75,mean(mean(DUR_med_E_IEI,2)),0.3);
errorbar(2.75,mean(mean(DUR_med_E_IEI,2)),std(mean(DUR_med_E_IEI,2))/size(DUR_med_E_IEI,1));
boxwhisker(DUR_med_E_IEI,2.75);
bar(3.25,mean(mean(DUR_med_I_IEI,2)),0.3);
errorbar(3.25,mean(mean(DUR_med_I_IEI,2)),std(mean(DUR_med_I_IEI,2))/size(DUR_med_I_IEI,1));
boxwhisker(DUR_med_I_IEI,3.25);
bar(4.75,mean(mean(POST_med_E_IEI,2)),0.3);
errorbar(4.75,mean(mean(POST_med_E_IEI,2)),std(mean(POST_med_E_IEI,2))/size(POST_med_E_IEI,1));
boxwhisker(POST_med_E_IEI,4.75);
bar(5.25,mean(mean(POST_med_I_IEI,2)),0.3);
errorbar(5.25,mean(mean(POST_med_I_IEI,2)),std(mean(POST_med_I_IEI,2))/size(POST_med_I_IEI,1));
boxwhisker(POST_med_I_IEI,5.25);
%case 2 -- pooled across all cells
figure;
%scatter(0.75+0.2.*(rand(length(PRE_collated_E_IEI),1)-0.5),PRE_collated_E_IEI,'b');
boxwhisker(PRE_collated_E_IEI,0.75,[0 0 1]);
hold on;
%scatter(1.25+0.2.*(rand(length(PRE_collated_I_IEI),1)-0.5),PRE_collated_I_IEI,'r');
boxwhisker(PRE_collated_I_IEI,1.25,[1 0 0]);
%scatter(2.75+0.2.*(rand(length(DUR_collated_E_IEI),1)-0.5),DUR_collated_E_IEI,'b');
boxwhisker(DUR_collated_E_IEI,2.75,[0 0 1]);
%scatter(3.25+0.2.*(rand(length(DUR_collated_I_IEI),1)-0.5),DUR_collated_I_IEI,'r');
boxwhisker(DUR_collated_I_IEI,3.25,[1 0 0]);
%scatter(4.75+0.2.*(rand(length(POST_collated_E_IEI),1)-0.5),POST_collated_E_IEI,'b');
boxwhisker(POST_collated_E_IEI,4.75,[0 0 1]);
%scatter(5.25+0.2.*(rand(length(POST_collated_I_IEI),1)-0.5),POST_collated_I_IEI,'r');
boxwhisker(POST_collated_I_IEI,5.25,[1 0 0]);
ax = gca;
ax.XTick = [1,3,5];
ax.XTickLabels = {'0 days','0-5 days','>5 days'};
ylim([0 0.1]);
xlabel('Experience');
ylabel('IEI (seconds)');

%Amplitudes (cdf)
f = figure;
plot_step = [0:0.1:1500];
bin_centers_a = (plot_step(1:end-1)+plot_step(2:end))/2;
pre_E_amps = [];
dur_E_amps = [];
post_E_amps = [];
pre_I_amps = [];
dur_I_amps = [];
post_I_amps = [];
for m = 1:length(matched_order_a_experience_E),
    if matched_order_a_experience_E(1,m) <= 0,
        pre_E_amps(end+1,1) = collated_E_alt_amps(m,1);
    elseif (matched_order_a_experience_E(1,m) > 0)&&(matched_order_a_experience_E(1,m) < 5),
        dur_E_amps(end+1,1) = collated_E_alt_amps(m,1);
    elseif matched_order_a_experience_E(1,m) > 5,
        post_E_amps(end+1,1) = collated_E_alt_amps(m,1);
    else
    end
end
for n = 1:length(matched_order_a_experience_I),
    if matched_order_a_experience_I(1,n) <= 0,
        pre_I_amps(end+1,1) = collated_I_alt_amps(n,1);
    elseif (matched_order_a_experience_I(1,n) > 0)&&(matched_order_a_experience_I(1,n) < 5),
        dur_I_amps(end+1,1) = collated_I_alt_amps(n,1);
    elseif matched_order_a_experience_I(1,n) > 5,
        post_I_amps(end+1,1) = collated_I_alt_amps(n,1);
    else
    end
end
c_factor = 1e12;
pre_E_amps = cumsum(histc((c_factor.*pre_E_amps),bin_centers_a));
dur_E_amps = cumsum(histc((c_factor.*dur_E_amps),bin_centers_a));
post_E_amps = cumsum(histc((c_factor.*post_E_amps),bin_centers_a));
pre_I_amps = cumsum(histc((c_factor.*pre_I_amps),bin_centers_a));
dur_I_amps = cumsum(histc((c_factor.*dur_I_amps),bin_centers_a));
post_I_amps = cumsum(histc((c_factor.*post_I_amps),bin_centers_a));
subplot(1,2,1);
plot(bin_centers_a,(pre_E_amps/max(pre_E_amps)),'b--');
hold on;
plot(bin_centers_a,(dur_E_amps/max(dur_E_amps)),'b-');
plot(bin_centers_a,(post_E_amps/max(post_E_amps)),'b:');
xlabel('Spontaneous event amplitudes (picoAmps)');
ylabel('Proportion');
if use_age == 0,
    legend('pre-opening','0-5 days open','>5 days open','Location','East');
else
    legend('<p30','p30-35','>p35','Location','East');
end
hold off
subplot(1,2,2);
plot(bin_centers_a,(pre_I_amps/max(pre_I_amps)),'r--');
hold on;
plot(bin_centers_a,(dur_I_amps/max(dur_I_amps)),'r-');
plot(bin_centers_a,(post_I_amps/max(post_I_amps)),'r:');
xlabel('Spontaneous event amplitudes (picoAmps)');
ylabel('Proportion');
if use_age == 0,
    legend('pre-opening','0-5 days open','>5 days open','Location','East');
else
    legend('<p30','p30-35','>p35','Location','East');
end
hold off;
%Scatter amplitudes
f = figure;
scatter((matched_order_a_experience_E)+0.1,(collated_E_alt_amps.*c_factor),'b');
hold on;
scatter((matched_order_a_experience_I)-0.1,(collated_I_alt_amps.*c_factor),'r');
xlabel('Days of experience');
ylabel('Mean Amplitude (pA)');
title('Trial mean event amplitudes');
legend('Exc.','Inh.','Location','NorthEast');
hold off;

%BY CELL
pre_GE_est = [];
pre_GI_est = [];
dur_GE_est = [];
dur_GI_est = [];
post_GE_est = [];
post_GI_est = [];
for i = 1:length(matched_order_cell_experience)
    if matched_order_cell_experience(i,1) <= 0,
        pre_GE_est(end+1,1) = nanmedian(cell_E_alt_amps{i,1}./E_drive(i,1));
        pre_GI_est(end+1,1) = nanmedian(cell_I_alt_amps{i,1}./I_drive(i,1));
    elseif (matched_order_cell_experience(i,1) > 0)&&(matched_order_cell_experience(i,1) < 5),
        dur_GE_est(end+1,1) = nanmedian(cell_E_alt_amps{i,1}./E_drive(i,1));
        dur_GI_est(end+1,1) = nanmedian(cell_I_alt_amps{i,1}./I_drive(i,1));
    elseif matched_order_cell_experience(i,1) >= 5,
        post_GE_est(end+1,1) = nanmedian(cell_E_alt_amps{i,1}./E_drive(i,1));
        post_GI_est(end+1,1) = nanmedian(cell_I_alt_amps{i,1}./I_drive(i,1));
    end
end
f = figure;
boxwhisker(pre_GE_est,0.75,[0 0 1]);
hold on;
boxwhisker(pre_GI_est,1.25,[1 0 0]);
boxwhisker(dur_GE_est,2.75,[0 0 1]);
boxwhisker(dur_GI_est,3.25,[1 0 0]);
boxwhisker(post_GE_est,4.75,[0 0 1]);
boxwhisker(post_GI_est,5.25,[1 0 0]);
scatter(0.75.*ones(length(pre_GE_est),1),pre_GE_est,'k');
scatter(1.25.*ones(length(pre_GI_est),1),pre_GI_est,'k');
scatter(2.75.*ones(length(dur_GE_est),1),dur_GE_est,'k');
scatter(3.25.*ones(length(dur_GI_est),1),dur_GI_est,'k');
scatter(4.75.*ones(length(post_GE_est),1),post_GE_est,'k');
scatter(5.25.*ones(length(post_GI_est),1),post_GI_est,'k');
set(gca,'XTick',[1,3,5],'XTickLabels',{'Pre-opening','0-5 days','>5 days'});
xlabel('Experience');
ylabel('Current amplitude / Driving force (nS)');
xlim([0 6]);
%K-W omnibus
gEdata = [pre_GE_est;dur_GE_est;post_GE_est];
gIdata = [pre_GI_est;dur_GI_est;post_GI_est];
ggroups = [ones(length(pre_GE_est),1);2.*ones(length(dur_GE_est),1);3.*ones(length(post_GE_est),1)];
[gEp,anovatabGE,stats_GE] = kruskalwallis(gEdata,ggroups);
[gIp,anovatabGI,stats_GI] = kruskalwallis(gIdata,ggroups);
%B-c Mann-Whitney
alpha = 0.05;
groups = 3;
corrected_alpha = alpha/((groups*(groups-1))/2);
G_E_early = ranksum(pre_GE_est,dur_GE_est);
G_E_late = ranksum(dur_GE_est,post_GE_est);
G_E_distal = ranksum(pre_GE_est,post_GE_est);
G_I_early = ranksum(pre_GI_est,dur_GI_est);
G_I_late = ranksum(dur_GI_est,post_GI_est);
G_I_distal = ranksum(pre_GI_est,post_GI_est);
%in time
med_E_G = [];
med_I_G = [];
for i = 1:length(matched_order_cell_experience),
    med_E_G = [med_E_G;nanmedian(cell_E_alt_amps{i,1}./E_drive(i,1))];
    med_I_G = [med_I_G;nanmedian(cell_I_alt_amps{i,1}./I_drive(i,1))];
end
f2 = figure;
scatter(matched_order_cell_experience,med_E_G,'b');
hold on;
scatter(matched_order_cell_experience,med_I_G,'r');
xlabel('Experience (in days)');
ylabel('Amplitude / Driving force (nS)');
xlim([0 20]);



%SKETCH AREA**************

pre_E_amps = cumsum(histc(pre_E_amps,bin_centers_a));
dur_E_amps = cumsum(histc(dur_E_amps,bin_centers_a));
post_E_amps = cumsum(histc(post_E_amps,bin_centers_a));
pre_I_amps = cumsum(histc(pre_I_amps,bin_centers_a));
dur_I_amps = cumsum(histc(dur_I_amps,bin_centers_a));
post_I_amps = cumsum(histc(post_I_amps,bin_centers_a));

f = figure;
scatter((matched_order_a_experience_E)+0.1,(collated_E_alt_amps),'b');
hold on;
scatter((matched_order_a_experience_I)-0.1,(collated_I_alt_amps),'r');
xlabel('Days of experience');
ylabel('Mean Amplitude / Driving force');
title('Estimated trial mean peak conductances');
legend('Exc.','Inh.','Location','NorthEast');
hold off;
 
pre_E_ = [];
dur_E_ = [];
post_E_ = [];
for i = 1:length(matched_order_a_experience_E),
    if matched_order_a_experience_E(1,i) <= 0,
        pre_E_(end+1,1) = collated_E_alt_amps(i,1);
    elseif (matched_order_a_experience_E(1,i) > 0)&&(matched_order_a_experience_E(1,i) <= 5),
        dur_E_(end+1,1) = collated_E_alt_amps(i,1);
    elseif matched_order_a_experience_E(1,i) > 5,
        post_E_(end+1,1) = collated_E_alt_amps(i,1);
    else
    end
end

            
            
            