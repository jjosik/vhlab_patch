
%
%
use_age = 0;  %recommend setting as function input
opening_onset_age = 30;
opening_estab_age = 35;
cell_results = uigetdir();
files = dir(fullfile(cell_results,'*.mat'));

PRE_collated_E_IEI = [];
PRE_collated_I_IEI = [];
DUR_collated_E_IEI = [];
DUR_collated_I_IEI = [];
POST_collated_E_IEI = [];
POST_collated_I_IEI = [];
PRE_med_E_IEI = [];  %first column is median, second is 25th p-tile, third is 75th
PRE_med_I_IEI = [];
DUR_med_E_IEI = [];
DUR_med_I_IEI = [];
POST_med_E_IEI = [];
POST_med_I_IEI = [];
PRE_collated_E_fmean = [];
PRE_collated_I_fmean = [];
DUR_collated_E_fmean = [];
DUR_collated_I_fmean = [];
POST_collated_E_fmean = [];
POST_collated_I_fmean = [];
PRE_cell_E_IEI = [];
PRE_cell_I_IEI = [];
DUR_cell_E_IEI = [];
DUR_cell_I_IEI = [];
POST_cell_E_IEI = [];
POST_cell_I_IEI = [];
for i = 1:length(files),
    cell_frequencydata = load(files(i).name,'cell_frequencydata');
    for j = 1:size(cell_frequencydata.cell_frequencydata,1),
        if use_age == 1,
            if cell_frequencydata.cell_frequencydata{j,1}.cell_ages < opening_onset_age,
                PRE_collated_E_IEI = [PRE_collated_E_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI];
                PRE_collated_I_IEI = [PRE_collated_I_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI];
                PRE_med_E_IEI(end+1,1) = median(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI);
                PRE_med_I_IEI(end+1,1) = median(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI);
                PRE_med_E_IEI(end,2) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI,25);
                PRE_med_I_IEI(end,2) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI,25);
                PRE_med_E_IEI(end,3) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI,75);
                PRE_med_I_IEI(end,3) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI,75);
                for z = 1:size(cell_frequencydata.cell_frequencydata{1,1}.cell_total_E_fmean,1),
                    if z == 1,
                        PRE_collated_E_fmean(end+1,1) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean(z,:));
                        PRE_collated_I_fmean(end+1,1) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean(z,:));
                    else
                        PRE_collated_E_fmean(end,z) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean(z,:));
                        PRE_collated_I_fmean(end,z) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean(z,:));
                    end
                end
            elseif (cell_frequencydata.cell_frequencydata{j,1}.cell_ages >= opening_onset_age)&&...
                    (cell_frequencydata.cell_frequencydata{j,1}.cell_ages <= opening_estab_age),
                DUR_collated_E_IEI = [DUR_collated_E_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI];
                DUR_collated_I_IEI = [DUR_collated_I_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI];
                DUR_med_E_IEI(end+1,1) = median(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI);
                DUR_med_I_IEI(end+1,1) = median(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI);
                DUR_med_E_IEI(end,2) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI,25);
                DUR_med_I_IEI(end,2) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI,25);
                DUR_med_E_IEI(end,3) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI,75);
                DUR_med_I_IEI(end,3) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI,75);
                for z = 1:size(cell_frequencydata.cell_frequencydata{1,1}.cell_total_E_fmean,1),
                    if z == 1,
                        DUR_collated_E_fmean(end+1,1) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean(z,:));
                        DUR_collated_I_fmean(end+1,1) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean(z,:));
                    else
                        DUR_collated_E_fmean(end,z) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean(z,:));
                        DUR_collated_I_fmean(end,z) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean(z,:));
                    end
                end
            elseif cell_frequencydata.cell_frequencydata{j,1}.cell_ages > opening_estab_age,
                POST_collated_E_IEI = [POST_collated_E_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI];
                POST_collated_I_IEI = [POST_collated_I_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI];
                POST_med_E_IEI(end+1,1) = median(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI);
                POST_med_I_IEI(end+1,1) = median(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI);
                POST_med_E_IEI(end,2) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI,25);
                POST_med_I_IEI(end,2) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI,25);
                POST_med_E_IEI(end,3) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI,75);
                POST_med_I_IEI(end,3) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI,75);
                for z = 1:size(cell_frequencydata.cell_frequencydata{1,1}.cell_total_E_fmean,1),
                    if z == 1,
                        POST_collated_E_fmean(end+1,1) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean(z,:));
                        POST_collated_I_fmean(end+1,1) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean(z,:));
                    else
                        POST_collated_E_fmean(end,z) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean(z,:));
                        POST_collated_I_fmean(end,z) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean(z,:));
                    end
                end
            end
        else
            if cell_frequencydata.cell_frequencydata{j,1}.cell_experience == 0,
                PRE_collated_E_IEI = [PRE_collated_E_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI];
                PRE_collated_I_IEI = [PRE_collated_I_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI];
                PRE_med_E_IEI(end+1,1) = median(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI);
                PRE_med_I_IEI(end+1,1) = median(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI);
                PRE_med_E_IEI(end,2) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI,25);
                PRE_med_I_IEI(end,2) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI,25);
                PRE_med_E_IEI(end,3) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI,75);
                PRE_med_I_IEI(end,3) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI,75);
                for z = 1:size(cell_frequencydata.cell_frequencydata{1,1}.cell_total_E_fmean,1),
                    if z == 1,
                        PRE_collated_E_fmean(end+1,1) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean(z,:));
                        PRE_collated_I_fmean(end+1,1) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean(z,:));
                    else
                        PRE_collated_E_fmean(end,z) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean(z,:));
                        PRE_collated_I_fmean(end,z) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean(z,:));
                    end
                end
            elseif (cell_frequencydata.cell_frequencydata{j,1}.cell_experience > 0)&&...
                    (cell_frequencydata.cell_frequencydata{j,1}.cell_experience <= 5),
                DUR_collated_E_IEI = [DUR_collated_E_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI];
                DUR_collated_I_IEI = [DUR_collated_I_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI];
                DUR_med_E_IEI(end+1,1) = median(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI);
                DUR_med_I_IEI(end+1,1) = median(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI);
                DUR_med_E_IEI(end,2) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI,25);
                DUR_med_I_IEI(end,2) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI,25);
                DUR_med_E_IEI(end,3) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI,75);
                DUR_med_I_IEI(end,3) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI,75);
                for z = 1:size(cell_frequencydata.cell_frequencydata{1,1}.cell_total_E_fmean,1),
                    if z == 1,
                        DUR_collated_E_fmean(end+1,1) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean(z,:));
                        DUR_collated_I_fmean(end+1,1) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean(z,:));
                    else
                        DUR_collated_E_fmean(end,z) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean(z,:));
                        DUR_collated_I_fmean(end,z) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean(z,:));
                    end
                end
            else
                POST_collated_E_IEI = [POST_collated_E_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI];
                POST_collated_I_IEI = [POST_collated_I_IEI;cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI];
                POST_med_E_IEI(end+1,1) = median(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI);
                POST_med_I_IEI(end+1,1) = median(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI);
                POST_med_E_IEI(end,2) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI,25);
                POST_med_I_IEI(end,2) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI,25);
                POST_med_E_IEI(end,3) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_IEI,75);
                POST_med_I_IEI(end,3) = prctile(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_IEI,75);
                for z = 1:size(cell_frequencydata.cell_frequencydata{1,1}.cell_total_E_fmean,1),
                    if z == 1,
                        POST_collated_E_fmean(end+1,1) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean(z,:));
                        POST_collated_I_fmean(end+1,1) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean(z,:));
                    else
                        POST_collated_E_fmean(end,z) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_E_fmean(z,:));
                        POST_collated_I_fmean(end,z) = mean(cell_frequencydata.cell_frequencydata{j,1}.cell_total_I_fmean(z,:));
                    end
                end
            end
        end
    end
end

%kruskal-wallis
E_dists = [PRE_collated_E_IEI;DUR_collated_E_IEI;POST_collated_E_IEI];
E_dists = [E_dists,[ones(length(PRE_collated_E_IEI),1);(2*ones(length(DUR_collated_E_IEI),1));(3*ones(length(POST_collated_E_IEI),1))]];
I_dists = [PRE_collated_I_IEI;DUR_collated_I_IEI;POST_collated_I_IEI];
I_dists = [I_dists,[ones(length(PRE_collated_I_IEI),1);(2*ones(length(DUR_collated_I_IEI),1));(3*ones(length(POST_collated_I_IEI),1))]];
[p_E,anovatab_E,stats_e] = kruskalwallis(E_dists(:,1),E_dists(:,2));
[p_I,anovatab_I,stats_I] = kruskalwallis(I_dists(:,1),I_dists(:,2));

%bonferroni-corrected Mann-Whitney
alpha = 0.05;
groups = 3;
corrected_alpha = alpha/((groups*(groups-1))/2);
p_E_early = ranksum(PRE_collated_E_IEI,DUR_collated_E_IEI);
p_E_late = ranksum(DUR_collated_E_IEI,POST_collated_E_IEI);
p_E_distal = ranksum(PRE_collated_E_IEI,POST_collated_E_IEI);
p_I_early = ranksum(PRE_collated_I_IEI,DUR_collated_I_IEI);
p_I_late = ranksum(DUR_collated_I_IEI,POST_collated_I_IEI);
p_I_distal = ranksum(PRE_collated_I_IEI,POST_collated_I_IEI);

%anova on medians
E_dists = [PRE_med_E_IEI(:,1);DUR_med_E_IEI(:,1);POST_med_E_IEI(:,1)];
E_dists = [E_dists,[ones(length(PRE_med_E_IEI),1);(2*ones(length(DUR_med_E_IEI),1));(3*ones(length(POST_med_E_IEI),1))]];
I_dists = [PRE_med_I_IEI(:,1);DUR_med_I_IEI(:,1);POST_med_I_IEI(:,1)];
I_dists = [I_dists,[ones(length(PRE_med_I_IEI),1);(2*ones(length(DUR_med_I_IEI),1));(3*ones(length(POST_med_I_IEI),1))]];
[p_m_E,anovatab_m_E,stats_E_] = anova1(E_dists(:,1),E_dists(:,2));
[p_m_I,anovatab_m_I,stats_I_] = anova1(I_dists(:,1),I_dists(:,2));

%bonferroni-corrected t-test on medians (unclear if normality assumption is
%appropriate)
alpha = 0.05;
groups = 3;
corrected_alpha = alpha/((groups*(groups-1))/2);
p_med_E_early = ttest(PRE_med_E_IEI(:,1),


%test of power for anova on medians

%*********************************
%***********PLOTS*****************
%*********************************
%cdfs
f = figure;
plot_step = [0:0.001:1.0];
bin_centers = (plot_step(1:end-1)+plot_step(2:end))/2;
pre_E = cumsum(histc(PRE_collated_E_IEI,bin_centers));
dur_E = cumsum(histc(DUR_collated_E_IEI,bin_centers));
post_E = cumsum(histc(POST_collated_E_IEI,bin_centers));
pre_I = cumsum(histc(PRE_collated_I_IEI,bin_centers));
dur_I = cumsum(histc(DUR_collated_I_IEI,bin_centers));
post_I = cumsum(histc(POST_collated_I_IEI,bin_centers));
subplot(1,2,1);
plot(bin_centers,(pre_E/max(pre_E)),'b--');
hold on;
plot(bin_centers,(dur_E/max(dur_E)),'b-');
plot(bin_centers,(post_E/max(post_E)),'b:');
xlim([-0.025 0.2]);
xlabel('Excitatory IEIs');
ylabel('Proportion');
if use_age == 0,
    legend('pre-opening','0-5 days open','>5 days open','Location','East');
else
    legend('<p30','p30-35','>p35','Location','East');
end
hold off
subplot(1,2,2);
plot(bin_centers,(pre_I/max(pre_I)),'r--');
hold on;
plot(bin_centers,(dur_I/max(dur_I)),'r-');
plot(bin_centers,(post_I/max(post_I)),'r:');
xlim([-0.025 0.2]);
xlabel('Inhibitory IEIs');
ylabel('Proportion');
if use_age == 0,
    legend('pre-opening','0-5 days open','>5 days open','Location','East');
else
    legend('<p30','p30-35','>p35','Location','East');
end
hold off;

%pooled medians
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

%average of cell median IEIs barplot (to match anova on medians)
f1 = figure;
root1 = 1.5;
root2 = 2.5;
root3 = 3.5;
shift = 0.1;
e_1 = bar((root1-shift),mean(PRE_med_E_IEI(:,1)),0.2,'b');
hold on;
i_1 = bar((root1+shift),mean(PRE_med_I_IEI(:,1)),0.2,'r');
e_2 = bar((root2-shift),mean(DUR_med_E_IEI(:,1)),0.2,'b');
i_2 = bar((root2+shift),mean(DUR_med_I_IEI(:,1)),0.2,'r');
e_3 = bar((root3-shift),mean(POST_med_E_IEI(:,1)),0.2,'b');
i_3 = bar((root3+shift),mean(POST_med_I_IEI(:,1)),0.2,'r');
set(gca,'XTick',[]);
ylabel('Median PSC Inter-event intervals (in sec.)');
legend('Excitatory','Inhibitory','Location','EastOutside');
scatter(repmat((root1-shift),length(PRE_med_E_IEI),1),PRE_med_E_IEI(:,1),'k');
scatter(repmat((root1+shift),length(PRE_med_I_IEI),1),PRE_med_I_IEI(:,1),'k');
scatter(repmat((root2-shift),length(DUR_med_E_IEI),1),DUR_med_E_IEI(:,1),'k');
scatter(repmat((root2+shift),length(DUR_med_I_IEI),1),DUR_med_I_IEI(:,1),'k');
scatter(repmat((root3-shift),length(POST_med_E_IEI),1),POST_med_E_IEI(:,1),'k');
scatter(repmat((root3+shift),length(POST_med_I_IEI),1),POST_med_I_IEI(:,1),'k');
lower_e_1 = scatter((root1-shift),mean(PRE_med_E_IEI(:,2)),'+','k');
upper_e_1 = scatter((root1-shift),mean(PRE_med_E_IEI(:,3)),'+','k');
lower_i_1 = scatter((root1+shift),mean(PRE_med_I_IEI(:,2)),'+','k');
upper_i_1 = scatter((root1+shift),mean(PRE_med_I_IEI(:,3)),'+','k');
lower_e_2 = scatter((root2-shift),mean(DUR_med_E_IEI(:,2)),'+','k');
upper_e_2 = scatter((root2-shift),mean(DUR_med_E_IEI(:,3)),'+','k');
lower_i_2 = scatter((root2+shift),mean(DUR_med_I_IEI(:,2)),'+','k');
upper_i_2 = scatter((root2+shift),mean(DUR_med_I_IEI(:,3)),'+','k');
lower_e_3 = scatter((root3-shift),mean(POST_med_E_IEI(:,2)),'+','k');
upper_e_3 = scatter((root3-shift),mean(POST_med_E_IEI(:,3)),'+','k');
lower_i_3 = scatter((root3+shift),mean(POST_med_I_IEI(:,2)),'+','k');
upper_i_3 = scatter((root3+shift),mean(POST_med_I_IEI(:,3)),'+','k');
hold off;
%median of cell median IEIs barplot (to match kruskal-wallis on cell medians)
f3 = figure;
root1 = 1.5;
root2 = 2.5;
root3 = 3.5;
shift = 0.1;
e_1 = bar((root1-shift),median(PRE_med_E_IEI(:,1)),0.2,'b');
hold on;
i_1 = bar((root1+shift),median(PRE_med_I_IEI(:,1)),0.2,'r');
e_2 = bar((root2-shift),median(DUR_med_E_IEI(:,1)),0.2,'b');
i_2 = bar((root2+shift),median(DUR_med_I_IEI(:,1)),0.2,'r');
e_3 = bar((root3-shift),median(POST_med_E_IEI(:,1)),0.2,'b');
i_3 = bar((root3+shift),median(POST_med_I_IEI(:,1)),0.2,'r');
set(gca,'XTick',[]);
ylabel('Median PSC Inter-event intervals (in sec.)');
legend('Excitatory','Inhibitory','Location','EastOutside');
scatter(repmat((root1-shift),length(PRE_med_E_IEI),1),PRE_med_E_IEI(:,1),'k');
scatter(repmat((root1+shift),length(PRE_med_I_IEI),1),PRE_med_I_IEI(:,1),'k');
scatter(repmat((root2-shift),length(DUR_med_E_IEI),1),DUR_med_E_IEI(:,1),'k');
scatter(repmat((root2+shift),length(DUR_med_I_IEI),1),DUR_med_I_IEI(:,1),'k');
scatter(repmat((root3-shift),length(POST_med_E_IEI),1),POST_med_E_IEI(:,1),'k');
scatter(repmat((root3+shift),length(POST_med_I_IEI),1),POST_med_I_IEI(:,1),'k');
lower_e_1 = scatter((root1-shift),median(PRE_med_E_IEI(:,2)),'+','k');
upper_e_1 = scatter((root1-shift),median(PRE_med_E_IEI(:,3)),'+','k');
lower_i_1 = scatter((root1+shift),median(PRE_med_I_IEI(:,2)),'+','k');
upper_i_1 = scatter((root1+shift),median(PRE_med_I_IEI(:,3)),'+','k');
lower_e_2 = scatter((root2-shift),median(DUR_med_E_IEI(:,2)),'+','k');
upper_e_2 = scatter((root2-shift),median(DUR_med_E_IEI(:,3)),'+','k');
lower_i_2 = scatter((root2+shift),median(DUR_med_I_IEI(:,2)),'+','k');
upper_i_2 = scatter((root2+shift),median(DUR_med_I_IEI(:,3)),'+','k');
lower_e_3 = scatter((root3-shift),median(POST_med_E_IEI(:,2)),'+','k');
upper_e_3 = scatter((root3-shift),median(POST_med_E_IEI(:,3)),'+','k');
lower_i_3 = scatter((root3+shift),median(POST_med_I_IEI(:,2)),'+','k');
upper_i_3 = scatter((root3+shift),median(POST_med_I_IEI(:,3)),'+','k');
hold off;
%boxplots by cell
f20 = figure;
boxwhisker(PRE_med_E_IEI(:,1).*1000,0.75,[0 0 1]);
hold on;
boxwhisker(PRE_med_I_IEI(:,1).*1000,1.25,[1 0 0]);
boxwhisker(DUR_med_E_IEI(:,1).*1000,2.75,[0 0 1]);
boxwhisker(DUR_med_I_IEI(:,1).*1000,3.25,[1 0 0]);
boxwhisker(POST_med_E_IEI(:,1).*1000,4.75,[0 0 1]);
boxwhisker(POST_med_I_IEI(:,1).*1000,5.25,[1 0 0]);
xlabel('Experience');
set(gca,'XTick',[1,3,5],'XTickLabels',{'Pre-opening','0-5 days','>5 days'});
ylabel('Inter-event intervals (ms)');
scatter(0.75.*ones(length(PRE_med_E_IEI),1),PRE_med_E_IEI(:,1).*1000,'k');
scatter(1.25.*ones(length(PRE_med_I_IEI),1),PRE_med_I_IEI(:,1).*1000,'k');
scatter(2.75.*ones(length(DUR_med_E_IEI),1),DUR_med_E_IEI(:,1).*1000,'k');
scatter(3.25.*ones(length(DUR_med_I_IEI),1),DUR_med_I_IEI(:,1).*1000,'k');
scatter(4.75.*ones(length(POST_med_E_IEI),1),POST_med_E_IEI(:,1).*1000,'k');
scatter(5.25.*ones(length(POST_med_I_IEI),1),POST_med_I_IEI(:,1).*1000,'k');
xlim([0 6]);
%K-W omnibus
idataE = [PRE_med_E_IEI(:,1);DUR_med_E_IEI(:,1);POST_med_E_IEI(:,1)];
idataI = [PRE_med_I_IEI(:,1);DUR_med_I_IEI(:,1);POST_med_I_IEI(:,1)];
igroups = [ones(length(PRE_med_E_IEI),1);2.*ones(length(DUR_med_E_IEI),1);3.*ones(length(POST_med_E_IEI),1)];
[ipE,anovatabiE,stats_iE] = kruskalwallis(idataE,igroups);
[ipI,anovatabiI,stats_iI] = kruskalwallis(idataI,igroups);
%B-c Mann-Whitney
alpha = 0.05;
groups = 3;
corrected_alpha = alpha/((groups*(groups-1))/2);
pi_E_early = ranksum(PRE_med_E_IEI(:,1),DUR_med_E_IEI(:,1));
pi_E_late = ranksum(DUR_med_E_IEI(:,1),POST_med_E_IEI(:,1));
pi_E_distal = ranksum(PRE_med_E_IEI(:,1),POST_med_E_IEI(:,1));
pi_I_early = ranksum(PRE_med_I_IEI(:,1),DUR_med_I_IEI(:,1));
pi_I_late = ranksum(DUR_med_I_IEI(:,1),POST_med_I_IEI(:,1));
pi_I_distal = ranksum(PRE_med_I_IEI(:,1),POST_med_I_IEI(:,1));
%anova on medians


%average PSC event rate, differential binning
%normalization factor (1/bin_size) should be removed if cell-by-cell
%analysis is run again in future (change was implemented in lowest level code
%7/13/15, but was quicker to make adjustment here)
f2 = figure;
root1e = 1.0;
root1i = 1.5;
root2e = 2.5;
root2i = 3.0;
root3e = 4.0;
root3i = 4.5;
shift = 0.1;
e1_b1 = bar((root1e-shift),mean(PRE_collated_E_fmean(:,1))*(1/0.01),0.1,'FaceColor',[0 0 0.4]);
hold on;
e1_b2 = bar(root1e,mean(PRE_collated_E_fmean(:,2))*(1/0.1),0.1,'FaceColor',[0 0 0.7]);
e1_b3 = bar((root1e+shift),mean(PRE_collated_E_fmean(:,3)),0.1,'FaceColor',[0 0 1]);
i1_b1 = bar((root1i-shift),mean(PRE_collated_I_fmean(:,1))*(1/0.01),0.1,'FaceColor',[0.4 0 0]);
i1_b2 = bar(root1i,mean(PRE_collated_I_fmean(:,2))*(1/0.1),0.1,'FaceColor',[0.7 0 0]);
i1_b3 = bar((root1i+shift),mean(PRE_collated_I_fmean(:,3)),0.1,'FaceColor',[1 0 0]);
e2_b1 = bar((root2e-shift),mean(DUR_collated_E_fmean(:,1))*(1/0.01),0.1,'FaceColor',[0 0 0.4]);
e2_b2 = bar(root2e,mean(DUR_collated_E_fmean(:,2))*(1/0.1),0.1,'FaceColor',[0 0 0.7]);
e2_b3 = bar((root2e+shift),mean(DUR_collated_E_fmean(:,3)),0.1,'FaceColor',[0 0 1]);
i2_b1 = bar((root2i-shift),mean(DUR_collated_I_fmean(:,1))*(1/0.01),0.1,'FaceColor',[0.4 0 0]);
i2_b2 = bar(root2i,mean(DUR_collated_I_fmean(:,2))*(1/0.1),0.1,'FaceColor',[0.7 0 0]);
i2_b3 = bar((root2i+shift),mean(DUR_collated_I_fmean(:,3)),0.1,'FaceColor',[1 0 0]);
e3_b1 = bar((root3e-shift),mean(POST_collated_E_fmean(:,1))*(1/0.01),0.1,'FaceColor',[0 0 0.4]);
e3_b2 = bar(root3e,mean(POST_collated_E_fmean(:,2))*(1/0.1),0.1,'FaceColor',[0 0 0.7]);
e3_b3 = bar((root3e+shift),mean(POST_collated_E_fmean(:,3)),0.1,'FaceColor',[0 0 1]);
i3_b1 = bar((root3i-shift),mean(POST_collated_I_fmean(:,1))*(1/0.01),0.1,'FaceColor',[0.4 0 0]);
i3_b2 = bar(root3i,mean(POST_collated_I_fmean(:,2))*(1/0.1),0.1,'FaceColor',[0.7 0 0]);
i3_b3 = bar((root3i+shift),mean(POST_collated_I_fmean(:,3)),0.1,'FaceColor',[1 0 0]);
set(gca,'XTick',[1.25,2.75,4.25],'XTickLabels',{'Pre_opening','0-5 days','>5 days'});
ylabel('Mean PSC events per second');
legend('Excitatory, 10 ms bins','Excitatory, 100 ms bins','Excitatory, 1 sec',...
    'Inhibitory, 10 ms bins','Inhibitory, 100 ms bins','Inhibitory, 1 sec bins');
hold off;
f4 = figure;












