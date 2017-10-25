function [frequency_data,amplitude_data,duration_data,sampling_rate,is_E_test,acq_duration] = O_PSC_deconvolution( plot_it,save_it,bin_sizes,num_devs_E,num_devs_I,s_margin,windowed_baseline,th_follows_base,filter_window)
%Deconvolution method from Pernia-Andrade et al. 2012 (Biophysical Journal)
 
%First, setup output structures:
frequency_data = struct('binned_means',[],'binned_counts',[],'epoch_IEI',[],'alt_binned_means',[],'alt_binned_counts',[],'alt_epoch_IEI',[]);
amplitude_data = struct('deconv_trial_amplitudes',[],'deconv_trial_peaktimes',[],'alt_trial_amplitudes',[],'alt_trial_peaktimes',[]);
duration_data = struct('threshold_crossing_durations',[],'fwhm_durations',[],'risetimes',[]);

current_folder = pwd;
fullfilename = [current_folder filesep 'Clamp1_uncomp.ma'];
data = hdf5read(fullfilename,'/data');
tdata = hdf5read(fullfilename,'/info/1/values');
sampling_rate = hdf5read(fullfilename,'/info/2/DAQ/primary','rate');
acq_duration = length(tdata)/sampling_rate;
%in case of unintentional wrong sampling rates (i.e. anything other than
%10kHz)
if sampling_rate > 10000 && mod(sampling_rate,10000)==0,
    s = sampling_rate/10000;
    data_ds = data(1:s:end,2);
    tdata_ds = tdata(1:s:end);
    ds_flag = 1;
elseif sampling_rate > 10000 && mod(sampling_rate,10000)~=0,
    s_add = sampling_rate;
    while mod(s_add,10000)~=0,
        s_add = s_add+1;
    end
    sample_pad = zeros(s_add,1);
    data_tailed = [data(:,2);sample_pad];
    tdata_tailed = [tdata;sample_pad];
    data_ds = data_tailed(1:s_add,end,1);
    tdata_ds = tdata_tailed(1:s_add,end,1);
    ds_flag = 1;
elseif sampling_rate < 10000,
    disp('***ERROR: datafile returns unexpected/unsupported sampling rate.  Rate= ',num2str(sampling_rate),'***');
    return;
else
    ds_flag = 0;
end
        
test_type = hdf5read(fullfilename,'/info/2/DAQ/command','holding');

if test_type < -0.020,
    is_E_test = 1;
else
    is_E_test = 0;
end
%compute template -
%assumes double exponential conductance model 
if is_E_test,
    tau_rise = 0.9; %in ms, based on manual assessment
    tau_decay = 1.8; %in ms
    amp = -2e-11;
else
    tau_rise = 0.5;
    tau_decay = 5.5;
    amp = 7e-11;
end
div_amp = (tau_decay/tau_rise)^(tau_rise/(tau_rise-tau_decay));
if ds_flag == 1,
    rd = data_ds;
else
    rd = data(:,2);
end
ms_samp = sampling_rate*0.001;
predelay = tau_decay*1;
N = floor((tau_decay*5+predelay)*ms_samp);

for i = 1:N
    t = i*ms_samp;
    if(t >= predelay)
        tmp(i) = (amp/div_amp)*(-exp(-(t-predelay)/tau_rise))+exp(-(t-predelay)/tau_decay);
    else
        tmp(i) = 0;
    end;
end;
highpass = 600;  %default 0.1; possible these could be improved on or need adjustment for diff. datasets
if is_E_test,
    lowpass = 25;  %default 100
else
    lowpass = 25;  
end
H = fft(tmp,size(rd,1))';
R = fft(rd);
%deconvolution in frequency domain
D = R./H;
%filter
f = [0:size(rd,1)-1] * sampling_rate/size(rd,1);
smooth_transform = 1;
if smooth_transform == 0,
    D( f < lowpass | (highpass < f & f < sampling_rate-high_pass) | sampling_rate-lowpass < f) = 0;
else
    w = (1/sqrt(2*pi*highpass/sampling_rate)) * exp(-0.5*min([f;sampling_rate-f]/(highpass),[],1).^2);
    w( f < lowpass | sampling_rate - lowpass < f) = 0;
    D_filt = sampling_rate*w(:).*D;
end
%return to time domain
d = real(ifft(D_filt));

%determine threshold values 
step_range = abs(min(diff(d)));
[max_marker,m] = max(abs(d));
valence = sign(d(m));
%see alt. step options below
%if step_range < 1,
%    step_order = floor(log10(step_range));
%else
%    step_order = ceil(log10(step_range));
%end
%step = floor(10^(abs(step_order))*abs(step_range))*10^(step_order);
if max_marker < 1,
    max_order = floor(log10(max_marker));
else
    max_order = ceil(log10(max_marker));
end
max_range = ceil(10^(abs(max_order)+1)*abs(max_marker))*10^(max_order-1);
step = abs(max_range)/100;
edges = -(max_range):step:max_range;
hd = histc(d,edges);
% fit gaussian to amplitude histogram, use fitted sd coeffiencient to set 4 SDs. 
gauss = fittype('a*exp(-((x-b)/(c))^2)');
gauss_d = gauss;
fo = fitoptions(gauss_d);
[max_hd,location] = max(hd);
%position_test = hd(location,1);
width_test = 1*10^(max_order);
fo.StartPoint = [max_hd;0;width_test];
gauss_d = setoptions(gauss_d,fo);
[c,gof] = fit(edges',hd,gauss_d);
f = figure;
plot(edges,hd);
hold on;
plot(c);
if is_E_test,
    threshold = num_devs_E*c.c;  %recommended default: 2.5 SD
else
    threshold = num_devs_I*c.c;  %recommended default: 1.5 SD
end
line([threshold threshold],[0 max(hd)]);
line([-threshold -threshold],[0 max(hd)]);
if save_it == 1,
    saveas(gcf,[pwd filesep 'amp_hist_fit' '_filter_w_',num2str(filter_window), '.fig']);
    close(f);
else
end

%detect thresholded events in deconvolved trace
if windowed_baseline == 1
    win_d = 2001;
    %filter_shift = (win-1)/2;  
    baseline = medfilt1(data(:,2),win_d);
else
    baseline = median(data(:,2));
end

test_TH = threshold*valence;
d_pad = [d;+inf];
if valence == 1,
    pos = (d_pad(1:end-1) >= test_TH) & (d_pad(1:end-1) >  d_pad(2:end) ) & (d_pad(1:end-1) > [+inf;d_pad(1:end-2)]);
    ix  = (d_pad(1:end-1) >= test_TH) & (d_pad(1:end-1) == d_pad(2:end) ) & (d_pad(1:end-1) > [+inf;d_pad(1:end-2)]);
elseif valence == -1,
    pos = (d_pad(1:end-1) <= test_TH) & (d_pad(1:end-1) < d_pad(2:end) ) & (d_pad(1:end-1) < [+inf;d_pad(1:end-2)]);
    ix = (d_pad(1:end-1) <= test_TH) & (d_pad(1:end-1) == d_pad(2:end) ) & (d_pad(1:end-1) < [+inf;d_pad(1:end-2)]);
else
end

k   = 0; 
ix  = find(ix);
while (1) 
	k   = k+1;
	ix2 = ix + k;
	ix  = ix(d(ix2) <= d(ix));
	ix2 = ix + k;
	if (isempty(ix) || ~any(d(ix)==d(ix2))), 
        break; 
    end; 
end
pos = sort([find(pos); ix]);
%plot deconvolution detection 
f = figure;
plot(tdata,d);
hold on;
line([0 10],[test_TH test_TH]);
d_level = max(d)+(0.05*max(d));
d_lx = repmat(d_level,length(pos),1);
scatter(tdata(pos),d_lx,'x');
hold off;
if save_it == 1,
    saveas(gcf,[pwd filesep 'deconvolve_detection' '_filter_w_',num2str(filter_window), '.fig']);
    close(f);
else
end
%plot event detection
f = figure;
plot(tdata,data(:,2));
hold on;
if windowed_baseline == 1,
    base_x = (1/sampling_rate):(1/sampling_rate):10;
    h1 = plot(base_x,baseline);
else
    h1 = line([0 10],[baseline baseline]);
end
set(h1,'Color',[0 0 0],'LineWidth',2);
if is_E_test,
    level = min(data(:,2))-(0.05*min(data(:,2)));
else
    level = max(data(:,2))+(0.05*max(data(:,2)));
end
lx = repmat(level,length(pos),1);
scatter(tdata(pos),lx,'x');
xlabel('Time (sec.)');
ylabel('Syn. current (Amps)');
hold off;
if save_it == 1,
    saveas(gcf,[pwd filesep 'd_raw_I_thresholded_events' '_filter_w_',num2str(filter_window), '.fig']);
    close(f);
else
end

iei = diff(tdata(pos));
peak_times = tdata(pos);
total_count = length(pos);
for z = 1:length(bin_sizes),
    edges = [0:bin_sizes(1,z):(round(acq_duration/bin_sizes(1,z))*bin_sizes(1,z))];
    binned_means(z,1) = bin_sizes(1,z);
    binned_means(z,2) = mean(histc(peak_times(:,1),edges));
    binned_counts{z,1} = bin_sizes(1,z);
    binned_counts{z,2} = histc(peak_times(:,1),edges);
end
extra_counts = {acq_duration, total_count};
extra_means = [acq_duration, (total_count/acq_duration)];
frequency_data.binned_counts = cat(1,binned_counts,extra_counts);
frequency_data.binned_means = [binned_means; extra_means];
frequency_data.epoch_IEI = iei;


%**********find local extrema**********
if nargin < 5,
    s_margin = 0.002;  %default to 2 ms search window
else
end
for N = 1:length(pos),
    peak_search = (pos(N)-((s_margin/2)*sampling_rate)):(pos(N)+((s_margin/2)*sampling_rate));
    peak_search = nonzeros(peak_search.*(peak_search>=0));
    peak_search = nonzeros(peak_search.*(peak_search<=(length(data))));
    peak_search = reshape(peak_search,length(peak_search),1);
    if is_E_test,
        peak_series = (data(peak_search(1:end-1),2) > data(peak_search(2:end),2)) %& (data(peak_search(1:end-1),2) > [+inf;data(peak_search(1:end-2),2)]);
        peak_s_temp = peak_search(1:end-1,1);
        [peak_amp,ind] = max(data(peak_s_temp(peak_series),2));
    else
        peak_series = (data(peak_search(1:end-1),2) < data(peak_search(2:end),2)) %& (data(peak_search(1:end-1),2) < [+inf;data(peak_search(1:end-2),2)]);
        peak_s_temp = peak_search(1:end-1);
        [peak_amp,ind] = min(data(peak_s_temp(peak_series),2));
    end
    %if (peak_amp > 0) & (~isempty(peak_amp)),  ****NOTE: if conditional
    %statement is reinstated (lines 250-262), change N_peak_amps/N_peak_indices in this block to
    %pre_N_peak_amps/pre_N_peak_indices.
        N_peak_amps(N,1) = peak_amp;
        subset = peak_search(peak_series);
        N_peak_indices(N,1) = subset(ind);  %compute absolute index
    %else
        %pre_N_peak_amps(N,1) = NaN;
        %pre_N_peak_indices(N,1) = NaN;       
    %end
end
%N_peak_amps = ~isnan(pre_N_peak_amps);
%N_peak_indices = ~isnan(pre_N_peak_indices);
H_amps = abs(N_peak_amps);
new_step = max(H_amps)/100;
amp_bins = 0:new_step:(1.5*max(H_amps));
ac = cumsum(histc(H_amps,amp_bins));
f = figure;
plot(amp_bins,ac);
hold on;
xlabel('Peak amplitudes');
ylabel('Raw counts');
hold off;
if save_it == 1,
    saveas(gcf,[pwd filesep 'peak_amp_dist' '_filter_w_',num2str(filter_window), '.fig']);
    close(f);
else
end
amplitude_data.deconv_trial_amplitudes = N_peak_amps(:,1);
amplitude_data.deconv_trial_peaktimes = peak_times;

%**********find durations, rise and decay time constants**********
cutoff_f = 300;
%determine num of st. devs. to set threshold by first high-pass filtering
%signal (3rd-order Butterworth seems to keep ringing manageable, but will
%overestimate # of crossings except if used in combination with min_width cutoff ~7 ms:
[b,a] = butter(3,(cutoff_f/(sampling_rate/2)),'high');
data_hp = filter(b,a,data(:,2));
%hp_diff = std(data_hp);
hp_diff = max(abs(data_hp-median(data_hp)));
%take first deriv of the signal and smooth it using Savitsky-Golay filter
d_deriv = gradient(data(:,2))./(1/sampling_rate);
der_smooth = sgolayfilt(d_deriv,7,51);
slope_threshold = -1e-8;
%find peaks in slope of the signal
[der_peaks,inds] = findpeaks(der_smooth,'MinPeakProminence',abs(slope_threshold));

if windowed_baseline == 1,
    if ~(mod(filter_window,2)==1),       %obviously filter_window should be specified in *seconds*
        win = (filter_window*1000)+1;
    else
        win = filter_window*1000;
    end
    %filter_shift = (win-1)/2;    %these 3 lines provide shift correction
    %(which is generally not needed)
    %pad_end = data((end-filter_shift)+1:end,2);
    %corrected_data = [data(filter_shift+1:end,2);pad_end];
    filter_shift = (win-1)/2;     %ends of data trials are padded with mirrored data half the size of the sliding filter window
    pre_mirror_pad = flipud(data(1:filter_shift+1,2));
    post_mirror_pad = flipud(data((length(data)-filter_shift):end,2));
    corrected_data = [pre_mirror_pad;data(:,2);post_mirror_pad];
    baseline_overage = medfilt1(corrected_data,win);
    baseline = baseline_overage(length(pre_mirror_pad)+1:(length(corrected_data)-length(post_mirror_pad)));
else
    baseline = median(data(:,2));
end
flat_baseline = median(data(:,2));

if is_E_test,
    amp_threshold = baseline-(num_devs_E*median(abs(baseline-data(:,2)))/0.6745);
    adj_amp_threshold = baseline-(hp_diff/0.6745);
    if adj_amp_threshold > amp_threshold,
        amp_threshold = adj_amp_threshold;
    else
    end
else
    amp_threshold = baseline+(num_devs_I*median(abs(baseline-data(:,2)))/0.6745);
    adj_amp_threshold = baseline+(hp_diff/0.6745);
    if adj_amp_threshold < amp_threshold,
        amp_threshold = adj_amp_threshold;
    else
    end
end

format long;
data_smoothed = sgolayfilt(data(:,2),7,101);
dy = data_smoothed(:,1)-amp_threshold;
alt_dy = data_smoothed(:,1)-baseline;   %changed 9/10/15
dy_mult = dy(1:end-1).*dy(2:end);
alt_dy_mult = alt_dy(1:end-1).*alt_dy(2:end);   %9/10/15
crossings = find(dy_mult<0);
alt_crossings = find(alt_dy_mult<0);            %9/10/15
crossings_diffs = data(crossings(2:end),2)-data(crossings(1:end-1),2);
alt_crossings_diffs = data(alt_crossings(2:end),2)-data(alt_crossings(1:end-1),2);  %9/10/15
crossing_slopes = crossings_diffs.*sampling_rate;
alt_slopes = gradient(data_smoothed);           %9/11/15
crossing_alt_slopes = alt_slopes(crossings);    %9/11/15
alt_crossing_slopes = alt_crossings_diffs.*sampling_rate;  %9/10/15
f = figure;
plot(tdata,data(:,2));
hold on;
if windowed_baseline == 1,
    base_x = (1/sampling_rate):(1/sampling_rate):10;
    h1 = plot(base_x,baseline);
else
    h1 = line([0 10],[flat_baseline flat_baseline]);
end
if th_follows_base >= 1,
    h2 = plot(base_x,amp_threshold);
else
    h2 = line([0 10],[amp_threshold amp_threshold]);
end
set(h1,'Color',[0 0 0],'LineWidth',2);
set(h2,'Color',[1 0 0],'LineStyle','--');
plot(tdata,data_smoothed,'k');
scatter(tdata(crossings),data_smoothed(crossings),'o');
scatter(tdata(alt_crossings),data_smoothed(alt_crossings),'d'); %9/10/15
durations = [];
collected_peaks = zeros(length(crossings)-1,1);
collected_amps = zeros(length(crossings)-1,1);
for k=1:length(crossings)-1,
    if (~(is_E_test) & crossing_alt_slopes(k,1)>0)|(is_E_test & crossing_alt_slopes(k,1)<0),
        data_w = data_smoothed(crossings(k):crossings(k+1),1);
        durations(end+1,1) = tdata(crossings(k+1))-tdata(crossings(k));
        if is_E_test,
            [local_peak,I] = min(data_w);
        else
            [local_peak,I] = max(data_w);
        end
    else
        continue;
    end
    collected_peaks(k,1) = crossings(k)-1+I;
    collected_amps(k,1) = abs(local_peak);
end
[nondup_index_a,s_a] = find(collected_amps>0);
[nondup_index_b,s_b] = find(collected_peaks>0);
u_collected_amps = collected_amps(nondup_index_a,1);
u_collected_peaks = collected_peaks(nondup_index_b,1);
half_max = zeros(length(collected_amps),1);
pre_fwhm = zeros(length(collected_amps),1);
r_ten_values = 0.1.*collected_amps;          %10-90% rise time 9/28/16
r_ninety_values = 0.9.*collected_amps;       %10-90% rise time 9/28/16
for j=1:length(collected_amps),
    if collected_amps(j,1) == 0,
        continue
    else
        if is_E_test,
            half_max(j,1) = baseline(collected_peaks(j,1))-(abs(collected_amps(j,1))-abs(baseline(collected_peaks(j,1))))/2;
        else
            half_max(j,1) = baseline(collected_peaks(j,1))+(abs(collected_amps(j,1))-abs(baseline(collected_peaks(j,1))))/2;
        end
        rise_data_w = data_smoothed(crossings(j):collected_peaks(j,1),1);
        decay_data_w = data_smoothed(collected_peaks(j,1):crossings(j+1,1),1);
        point_check_rise = abs(rise_data_w-half_max(j,1));
        point_check_decay = abs(decay_data_w-half_max(j,1));
        [rise_value,rise_ind] = sort(point_check_rise);
        [decay_value,decay_ind] = sort(point_check_decay);
         [Rnew_ten_val,Rten_ind] = min(abs(rise_data_w-r_ten_values(j,1)));     %10-90% rise time 9/28/16
         [Rnew_ninety_val,Rninety_ind] = min(abs(rise_data_w-r_ninety_values(j,1)));        %10-90% rise time 9/28/16
         Rten_ninety_risetime(j,1) = (Rninety_ind-Rten_ind)*(1/sampling_rate);             %10-90% rise time 9/28/16
         d_crit_point = 1/(exp(1))*(abs(collected_amps-)
        pre_fwhm(j,1) = (length(rise_data_w)-rise_ind(1,1)+(decay_ind(1,1)-1))*(1/sampling_rate);
        if (rise_ind(1,1) == max(rise_ind))&&(decay_ind(1,1) == 1),
            sub_alt_crossings = nonzeros((alt_crossings<crossings(j)).*alt_crossings);
            base_distance = abs(sub_alt_crossings-crossings(j));
            [nearest_base,alt_ind] = sort(base_distance);
            lowrise_data_w = data_smoothed(alt_crossings(alt_ind(1,1)):crossings(j));
            adv = 0;
            while alt_crossings(alt_ind(1,1)+adv)<collected_peaks(j,1),
                adv = adv+1;
            end
            lowdecay_data_w = data_smoothed(collected_peaks(j,1):alt_crossings(alt_ind(1,1)+adv));
            new_point_check_rise = abs(lowrise_data_w-half_max(j,1));
            new_point_check_decay = abs(lowdecay_data_w-half_max(j,1));
            [new_rise_value,new_rise_ind] = sort(new_point_check_rise);
            [new_decay_value,new_decay_ind] = sort(new_point_check_decay);
            pre_fwhm(j,1) = (alt_crossings(alt_ind(1,1)+1)-(length(lowdecay_data_w)-new_decay_ind(1,1))-...
                (alt_crossings(alt_ind(1,1))+new_rise_ind(1,1)))*(1/sampling_rate);
        else
        end
    end
end
[nonz_index,s_c] = find(pre_fwhm>0);
fwhm = pre_fwhm(nonz_index,1);
scatter(tdata(u_collected_peaks),data(u_collected_peaks,2),'+');
hold off;
if save_it == 1,
    if windowed_baseline == 1,
        if th_follows_base >= 1,
            saveas (gcf,[pwd filesep 'c_adaptth_raw_I_thresholded_events' '_filter_w_',num2str(filter_window), '.fig']);
        else
            saveas (gcf,[pwd filesep 'c_adapt_raw_I_thresholded_events' '_filter_w_',num2str(filter_window), '.fig']);
        end
    else
        saveas(gcf,[pwd filesep 'c_raw_I_thresholded_events' '_filter_w_',num2str(filter_window), '.fig']);
    end
    close(f);
else
end
alt_iei = diff(tdata(u_collected_peaks));
alt_peak_times = tdata(u_collected_peaks);
alt_total_count = length(u_collected_peaks);
for z = 1:length(bin_sizes),
    alt_edges = [0:bin_sizes(1,z):(round(acq_duration/bin_sizes(1,z))*bin_sizes(1,z))];
    alt_binned_means(z,1) = bin_sizes(1,z);
    alt_binned_means(z,2) = mean(histc(alt_peak_times(:,1),edges));
    alt_binned_counts{z,1} = bin_sizes(1,z);
    alt_binned_counts{z,2} = histc(alt_peak_times(:,1),edges);
end
alt_extra_counts = {acq_duration, alt_total_count};
alt_extra_means = [acq_duration, (alt_total_count/acq_duration)];
frequency_data.alt_binned_counts = cat(1,alt_binned_counts,alt_extra_counts);
frequency_data.alt_binned_means = [alt_binned_means; alt_extra_means];
frequency_data.alt_epoch_IEI = alt_iei;

%durations histogram
edges = 0:0.0001:0.1;
bin_centers = (edges(1:end-1)+edges(2:end))/2;
dist = cumsum(histc(durations,bin_centers));
f1 = figure;
plot(bin_centers,(dist/max(dist)));
xlabel('Durations (s)');
ylabel('Cumulative proportion of all durations');
title('Threshold crossing-to-threshold crossing duration');
if save_it == 1,
    saveas(gcf,[pwd filesep 'durations_hist' '_filter_w_',num2str(filter_window), '.fig']);
    close(f1);
else
end
%fwhm histogram
edges_f = 0:0.0001:(max(fwhm)*1.5);
f_bin_centers = (edges_f(1:end-1)+edges_f(2:end))/2;
d_f = cumsum(histc(fwhm,f_bin_centers));
f2 = figure;
plot(f_bin_centers,(d_f/max(d_f)));
xlabel('FWHM (s)');
ylabel('Cumulative proportion of all durations');
title('FWHM duration');
if save_it == 1,
    saveas(gcf,[pwd filesep 'fwhm_hist' '_filter_w_',num2str(filter_window), '.fig']);
    close(f2);
else
end
%[nondup_index_a,s_a] = find(collected_amps~=0);
%[nondup_index_b,s_b] = find(collected_peaks~=0);
%u_collected_amps = abs(collected_amps(nondup_index_a,1));
%u_collected_peaks = abs(collected_peaks(nondup_index_b,1));
edges_a = 0:(max(u_collected_amps)/100):(max(u_collected_amps)*1.5);
a_bin_centers = (edges_a(1:end-1)+edges_a(2:end))/2;
d_a = cumsum(histc(u_collected_amps,a_bin_centers));
f3 = figure;
plot(a_bin_centers,(d_a/max(d_a)));
xlabel('Event amplitudes');
ylabel('Cumulative proportion');
title('Peak Amplitude distribution using first deriv method');
if save_it == 1,
    saveas(gcf,[pwd filesep 'event_amp_hist' '_filter_w_',num2str(filter_window), '.fig']);
    close(f3);
else
end
%risetime histogram     added 9/28/16
r_binsize = 2*(prctile(ten_ninety_risetime,75)-prctile(ten_ninety_risetime,25))*(1/(nthroot(length(ten_ninety_risetime),3)));
edges_r = 0:r_binsize:(max(ten_ninety_risetime)*1.5);
D_r = histc(ten_ninety_risetime,edges_r);
r_bin_centers = (edges_r(1:end-1)+edges_r(2:end))/2;
D_r = D_r(1:end-1);
f4 = figure;
bar(bin_centers,D_r*1000,'r');
hold on;
ylabel('Count');
xlabel('Rise time (msec)');
title('Distribution of PSC event 10-90% rise times');
if save_it == 1,
    saveas(gcf,[pwd filesep 'event_risetime' '_filter_w_',num2str(filter_window),'.fig']);
    close(f4);
else
end

amplitude_data.alt_trial_amplitudes = u_collected_amps;
amplitude_data.alt_trial_peaktimes = u_collected_peaks;
%duration_data = struct('threshold_crossing_durations',[],'fwhm_durations',[],'risetimes',[]);
duration_data.threshold_crossing_durations = durations;
duration_data.fwhm_durations = fwhm;
duration_data.risetimes = ten_ninety_risetime;

%ALT. METHOD: EMPLOYS MODEL FITS TO EACH EVENT****UNDER CONSTRUCTION****
%REQUIRES functions define_fit_EPSC.m, define_fit_IPSC.m, and
%lsfit_PSC_event (define_fit_EPSC and lsfit_PSC_event have been written)
%padsize = 100;
%prefix = +inf*ones(padsize,1);
%t_prefix = fliplr((-1:-1:-padsize)*(1/sampling_rate));
%suffix = +inf*ones(padsize,1);
%t_suffix = max(tdata)+((1:1:padsize)*(1/sampling_rate));
%data_pad = [prefix;data(:,2);suffix];
%tdata_pad = [t_prefix';tdata;t_suffix'];
%for j = 1:length(N_peak_indices),
%    trans_cross_pre = find(data_pad((N_peak_indices(j)+padsize)-(0.005*sampling_rate):(N_peak_indices(j)+padsize)) > baseline);
%    trans_cross_post = find(data_pad((N_peak_indices(j)+padsize):(N_peak_indices(j)+padsize)+(0.05*sampling_rate)) > baseline);
%    onset_pre = max(trans_cross_pre);
%    offset_post = min(trans_cross_post);
%    rel_absolute_onset(j) = (N_peak_indices(j)+padsize)-(0.005*sampling_rate)-onset_pre;
%    rel_absolute_offset(j) = (N_peak_indices(j)+padsize)+offset_post;
%    total_rise_time = tdata_pad(N_peak_indices(j)+padsize)-tdata_pad(rel_absolute_onset(j));
%    total_decay_time = tdata_pad(rel_absolute_offset(j))-tdata_pad(N_peak_indices(j)+padsize);
%    amp_base = abs(data_pad(N_peak_indices(j)+padsize)-data_pad(rel_absolute_onset(j)));
%    e_rise = amp_base*(1-1/exp(1));
%    e_decay = amp_base*(1/exp(1));
%    if is_E_test,
%        rise_point = data_pad(rel_absolute_onset(j))-e_rise;
%        decay_point = data_pad(N_peak_indices(j)+padsize)+e_decay;
%    else
%        rise_point = data_pad(rel_absolute_onset(j))+e_rise;
%        decay_point = data_pad(N_peak_indices(j)+padsize)-e_decay;
%    end
%    search_pre = abs(data_pad(rel_absolute_onset(j):N_peak_indices(j)+padsize)-rise_point);
%    search_post = abs(data_pad(N_peak_indices(j)+padsize:rel_absolute_offset(j))-decay_point);
%    [pre_nearest_point,pre_near_loc] = sort(search_pre);
%    [post_nearest_point,post_near_loc] = sort(search_post);
%    tau_rise(j) = pre_near_loc(1)*(1/sampling_rate);
%    tau_decay(j) = post_near_loc(1)*(1/sampling_rate);
%    absolute_onset(j) = rel_absolute_onset(j)-padsize;
%    absolute_offset(j) = rel_absolute_offset(j)-padsize;
%end
%use trans_pre_cross point to set up predelay for one-to-one event model
%fitting
%for k = 1:length(N_peak_indices),
%    start_index = absolute_onset(k);
%    plot_it = 1;
%    [hf,x] = lsfit_PSC_event(data,tdata,ms_samp,start_index,is_E_test,predelay,plot_it);
    

end
    