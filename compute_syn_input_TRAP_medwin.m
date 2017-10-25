h5function [ epoch_charge,sampling_rate,is_E_test,acq_duration ] = compute_syn_input_TRAP_medwin( plot_it,adapt_baseline,filter_win )
% COMPUTE_SYN_INPUT_TRAP - takes an identified voltage clamp file from
% spontaneous synaptic input protocol and returns the total charge
% accumulated in the given current trace.
%   ***INPUTS***
%       SPONTSYNI_FILE - filename of the current trace to be analyzed
%       IS_E_TEST - enter "1" if current trace file is test of excitatory
%           currents (holding potential at E_Cl); enter "0" if current trace
%           file is test of inhibitory currents (holding potential near 0 mV)
%       START_POINT - the beginning point, in seconds from the beginning of file, 
%           of the desired part of the recording to be analyzed; exclude if entire file is
%           desired.
%       END_POINT - the end point, in seconds from beginning of file, of
%           the desired part of the recording to be analyzed; again exclude
%           if want to analyze to end of file.
%   ***OUTPUT***
%       TOTAL_CHARGE - computed by setting baseline equal to the median of 
%           the total current trace, excluding all current deflections in the wrong 
%           direction (i.e. excitatory currrents are inward, so upward deflections are 
%           all set equal to the baseline), then computing the area under
%           the current trace with respect to the baseline.  The
%           composite Trapezoidal rule is the numeric integration method used.
%
%%%     NOTE -- needs correction; integral is defined with respect to zero,
%%%     not baseline -- need to include subtraction for this offset.
%%%     Probably explains results that are off by 6-7 orders of magnitude.
%%%     (get microCoulomb results rather than expected picroCoulombs)
%%%     UPDATE - the reason for this error had nothing to do with the
%%%     integration procedure -- previously failed to extract the 2nd
%%%     column (current amplitude) from the data structure at line 57.
%%%     Also corrected in other "simple y-summation" version.

current_folder = pwd;
fullfilename = [current_folder filesep 'Clamp1_uncomp.ma'];
data = hdf5read(fullfilename,'/data');
tdata = hdf5read(fullfilename,'/info/1/values');

test_type = hdf5read(fullfilename,'/info/2/Protocol/holding');
if test_type < -0.020,
    is_E_test = 1;
else 
    is_E_test = 0;
end

data_filtered = data;
%conditional block except for 'else' statement added for adaptive baseline
%variant
if adapt_baseline == 1,
    if ~(mod((filter_win*1000),2)==1),
        win = (filter_win*1000)+1;
    else
        win = filter_win*1000;
    end
    filter_shift = (win-1)/2;
    front_mirror_pad = flipud(data((1:filter_shift+1),2));
    back_mirror_pad = flipud(data((length(data)-filter_shift):end,2));
    extended_data = [front_mirror_pad;data(:,2);back_mirror_pad];
    baseline_overage = medfilt1(extended_data,win);
    baseline = baseline_overage(length(front_mirror_pad)+1:(length(extended_data)-length(back_mirror_pad)));
else
    baseline = median(data(:,2));
end
%count_threshold = 3*std(data(:,2),1,1);    %***original to non-adaptive
%baseline variant
if is_E_test == 1,
    for i = 1:length(tdata),
        if data(i,2) > baseline(i,1),
            data_filtered(i,2) = baseline(i,1);
            
        end
    end
    count_line = baseline-(1.5*median(abs(baseline(:,1)-data(:,2)))/0.6745);
    %count_line = baseline-count_threshold; %***
end
if ~(is_E_test == 1),
    for j = 1:length(tdata),
        if data(j,2) < baseline,
            data_filtered(j,2) = baseline(j,1);
        end
    end
    count_line = baseline+(1.5*median(abs(baseline(:,1)-data(:,2)))/0.6745);
    %count_line = baseline+count_threshold; %***
end
flat_baseline = median(data(:,2));

sampling_rate = hdf5read(fullfilename,'info/2/DAQ/primary/rate');
acq_duration = length(tdata)/sampling_rate;
y = data_filtered(:,2) - baseline;
x = tdata;
unit_integral = trapz(x,y);
count_index = find(abs(data_filtered(:,2)-count_line)<=10^-13);
%target_period_length = end_point - start_point;
%if target_period_length == acq_duration,
    epoch_charge = abs(unit_integral);  %mult. by dt not needed; x-component explicit calc. in trapz method
%end
%if target_period_length < acq_duration
%    fine_start_point = start_point*sampling_rate;
%    fine_end_point = end_point*sampling_rate;
%    unit_integral_woBEGIN = unit_integral(end) - unit_integral(fine_start_point);
%    unit_integral_woEND = unit_integral(end) - unit_integral(fine_end_point);
%    unit_integral_adjusted = unit_integral_woBEGIN - unit_integral_woEND;
%    total_charge = unit_integral_adjusted*1/sampling_rate;
%end
if plot_it == 1,
    f = figure;
    plot(tdata(:),data(:,2),'k');
    hold on;
    if adapt_baseline == 1,
        base_t = (1/sampling_rate):(1/sampling_rate):10;
        h1 = plot(base_t,baseline);
        h2 = plot(base_t,count_line);
    else
        h1 = line([0 10],[flat_baseline flat_baseline]);
        h2 = line([0 10],[count_line count_line]);
    end
    set(h1,'Color','r');
    set(h1,'LineStyle','-');
    set(h2,'LineStyle','--');
    %h.LineWeight = '2.0';
    xlabel('time (secs)');
    ylabel('Syn. current (Amps)');
    saveas(gcf,[pwd filesep 'raw_I' '.fig']);
    close(f);
else
end


end



