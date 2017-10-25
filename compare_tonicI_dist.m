function [f,H,p] = compare_tonicI_dist(cell_data_vector,cell_tdata_vector,on_time,off_time)

[y1,I_s] = min(abs(on_time-cell_tdata_vector));
[y2,I_e] = min(abs(off_time-cell_tdata_vector));
gaba_indices = I_s:I_e;
off_gaba_indices = [1:I_s-1,I_e:length(cell_tdata_vector)];
gaba_points = cell_data_vector(gaba_indices);
gaba_tpoints = cell_tdata_vector(gaba_indices);
off_gaba_points = cell_data_vector(off_gaba_indices);
off_gaba_tpoints = cell_tdata_vector(off_gaba_indices);
deflection_threshold = 5*std(cell_data_vector);
gaba_th = abs(gaba_points-median(gaba_points));
off_gaba_th = abs(off_gaba_points-median(off_gaba_points));
trimmed_gaba_points = gaba_points(find(gaba_th<deflection_threshold));
trimmed_off_gaba_points = off_gaba_points(find(off_gaba_th<deflection_threshold));

%g_low_range = min(gaba_points)-(max(gaba_points)-min(gaba_points));
%g_high_range = max(gaba_points)+(max(gaba_points)-min(gaba_points));
%og_low_range = min(off_gaba_points)-(max(off_gaba_points)-min(off_gaba_points));
%og_high_range = max(off_gaba_points)+(max(off_gaba_points)-min(off_gaba_points));
g_low_range = min(trimmed_gaba_points)-(max(trimmed_gaba_points)-min(trimmed_gaba_points));
g_high_range = max(trimmed_gaba_points)+(max(trimmed_gaba_points)-min(trimmed_gaba_points));
og_low_range = min(trimmed_off_gaba_points)-(max(trimmed_off_gaba_points)-min(trimmed_off_gaba_points));
og_high_range = max(trimmed_off_gaba_points)+(max(trimmed_off_gaba_points)-min(trimmed_off_gaba_points));
abs_low_range = min(g_low_range,og_low_range);
abs_high_range = max(g_high_range,og_high_range);
bin_step = 3e-12;
bin_edges = [abs_low_range:bin_step:abs_high_range];
%D_g = histc(gaba_points,bin_edges);
D_g = histc(trimmed_gaba_points,bin_edges);
%D_og = histc(off_gaba_points,bin_edges);
D_og = histc(trimmed_off_gaba_points,bin_edges);
bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2;
D_g = D_g(1:end-1);
D_og = D_og(1:end-1);
f = figure;
x_range_g = 4*std(gaba_points);
x_range_og = 4*std(off_gaba_points);
abs_x_range = max(x_range_g,x_range_og);
subplot(2,1,1);
bar(bin_centers,D_g,'r');
line([median(trimmed_gaba_points) median(trimmed_gaba_points)],[0 1.5e6]);
hold on;
ylabel('Count');
legend('GABA on','Location','East');
xlim([mean(trimmed_gaba_points)-abs_x_range mean(trimmed_gaba_points)+abs_x_range]);
subplot(2,1,2);
bar(bin_centers,D_og,'b');
line([median(trimmed_off_gaba_points) median(trimmed_off_gaba_points)],[0 1.5e6]);
ylabel('Count');
xlabel('Current values (Amps)');
xlim([mean(trimmed_gaba_points)-abs_x_range mean(trimmed_gaba_points)+abs_x_range]);
legend('GABA off','Location','East');  

%non_parametric
[H1,p1] = kstest2(D_g,D_og);
%parametric
[H2,p2] = ttest2(D_g,D_og);



