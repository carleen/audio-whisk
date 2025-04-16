% Determine the best candidates at time t_a fro the tempo and downbeat
% locations
function [beat_locations] = tempoDownbeatLocation(sum_flux, t_bins, t_a, ...
    time_window, N_D, show_plots)


    % Values for tempo, discrete
    N_R_bpm = (60:1:150);
    N_R_bps = N_R_bpm/60;
    
    % Create array of zeros for each of the candidate tempos, with N_D
    downbeat_locs = zeros(N_D, length(N_R_bpm));
    halfbeat_locs = zeros(N_D, length(N_R_bpm));
    quarterbeat_locs = zeros(N_D*2, length(N_R_bpm));
    
    for ii_bps=1:length(N_R_bps)
        cur_N_R_bps = N_R_bps(ii_bps);
        T_i = 1/cur_N_R_bps; % Candidate beat period
    
        candidate_downbeat_locs = zeros(N_D,1);
        for ii_j = 1:N_D
            t_i_j = t_a + T_i * ii_j/N_D;
            candidate_downbeat_locs(ii_j) = t_i_j;
        end
        
        downbeat_locs(:,ii_bps) = candidate_downbeat_locs;
        halfbeat_locs(:,ii_bps) = candidate_downbeat_locs + T_i/2;
    
        quarterbeat_locs(1:N_D, ii_bps) = downbeat_locs(:,ii_bps) + T_i/4;
        quarterbeat_locs(N_D+1:end, ii_bps) = halfbeat_locs(:,ii_bps) + T_i/4;
    
    end
    
    % For each tempo and each tempo downbeat location, get the cc_flux 
    
    cross_correlation_values = zeros(N_D, length(N_R_bpm));
    
    for ind_tempo = 1:length(N_R_bps)
    
        cur_downbeat_locs = downbeat_locs(:,ind_tempo);
        cur_tempo = N_R_bps(ind_tempo);
        
        
        cur_beat_rate = 1/cur_tempo;
        
        for ind_downbeat_candidate = 1:N_D
            cc_flux = 0;

            cur_start_downbeat = cur_downbeat_locs(ind_downbeat_candidate);
            
            cur_downbeat_times = (cur_start_downbeat:cur_beat_rate:cur_start_downbeat+time_window);
            start_time = cur_downbeat_times(1);
            end_time = cur_downbeat_times(end);
            
            % Get flux values within this window
            sum_flux_mask = t_bins < end_time & t_bins > start_time;
            sum_flux_at_window = sum_flux(sum_flux_mask);
            flux_time_values = t_bins(sum_flux_mask);
            
            % Combined time values
            all_downbeat_times = sort([cur_downbeat_times flux_time_values]);
            
            % Iterate over each of the downbeat times and get the nearest time
            % Mark off the downbeat candidates
            interpolated_flux = interp1(flux_time_values, sum_flux_at_window, all_downbeat_times, 'spline');
            % Normalized based on max value
            interpolated_flux = interpolated_flux/max(interpolated_flux);
            

            for ii_t = 1:length(cur_downbeat_times)
                ic = find(all_downbeat_times == cur_downbeat_times(ii_t));
                cc_flux = cc_flux + mean(interpolated_flux(ic));
            end
    
            % Do the same as above for the half beat times
            % Get flux values within this window
            cur_halfbeat_times = cur_downbeat_times + cur_beat_rate/2;
            start_time = cur_halfbeat_times(1);
            end_time = cur_halfbeat_times(end);
    
            sum_flux_mask = t_bins < end_time & t_bins > start_time;
            sum_flux_at_window = sum_flux(sum_flux_mask);
            flux_time_values = t_bins(sum_flux_mask);
            
            % Combined time values
            all_halfbeat_times = sort([cur_halfbeat_times flux_time_values]);
            
            % Iterate over each of the downbeat times and get the nearest time
            % Mark off the downbeat candidates
            interpolated_flux = interp1(flux_time_values, sum_flux_at_window, all_halfbeat_times, 'spline');
            % Normalized based on max value
            interpolated_flux = interpolated_flux/max(interpolated_flux);
    
            for ii_t = 1:length(cur_halfbeat_times)
                ic = find(all_halfbeat_times == cur_halfbeat_times(ii_t));
                cc_flux = cc_flux + mean(interpolated_flux(ic))*0.7;
            end

            % Do the same as above for the half beat times
            % Get flux values within this window
            cur_quarterbeat_times = [cur_downbeat_times + cur_beat_rate/4 cur_halfbeat_times + cur_beat_rate/4];
            cur_quarterbeat_times = sort(cur_quarterbeat_times);
            start_time = cur_quarterbeat_times(1);
            end_time = cur_quarterbeat_times(end);
    
            sum_flux_mask = t_bins < end_time & t_bins > start_time;
            sum_flux_at_window = sum_flux(sum_flux_mask);
            flux_time_values = t_bins(sum_flux_mask);
            
            % Combined time values
            all_quarterbeat_times = sort([cur_quarterbeat_times flux_time_values]);
            
            % Iterate over each of the downbeat times and get the nearest time
            % Mark off the downbeat candidates
            interpolated_flux = interp1(flux_time_values, sum_flux_at_window, all_quarterbeat_times, 'spline');
            % Normalized based on max value
            interpolated_flux = interpolated_flux/max(interpolated_flux);
    
            for ii_t = 1:length(cur_quarterbeat_times)
                ic = find(all_quarterbeat_times == cur_quarterbeat_times(ii_t));
                cc_flux = cc_flux + mean(interpolated_flux(ic))*0.1;
            end
    
            cross_correlation_values(ind_downbeat_candidate, ind_tempo) = cc_flux;
        end
    end
    
    % Find the 10 tempos with the highest values for the cross correlation
    cross_correlation_values_for_tempos = zeros(length(N_R_bpm),1);
    max_cc_start_time = zeros(length(N_R_bpm),1);
    
    for ii=1:length(cross_correlation_values_for_tempos)
        max_cc = max(cross_correlation_values(:,ii));
        cross_correlation_values_for_tempos(ii) = max_cc;
        indc = find(cross_correlation_values(:,ii) == max_cc);
        max_cc_start_time(ii) = downbeat_locs(indc, ii);
    end
    
    
    [sorted_cc, sorted_cc_indices] = sort(cross_correlation_values_for_tempos(:), 'descend');
    
    % Get the indices of the 10 greatest values
    top_indices = sorted_cc_indices(1:9);
    bpm_candidates = N_R_bpm(top_indices);
    start_time_for_candidates = max_cc_start_time(top_indices);

    if show_plots
    
        % For each of the candidates, get the start time and plot it with the flux
        f1=figure; hold on; grid on; 
        f1.Position = [680 254 1037 625];
        
        for ii=1:length(bpm_candidates)
            cur_bpm = bpm_candidates(ii);
            cur_start_time = start_time_for_candidates(ii);
            cur_beat_rate = 1/(cur_bpm/60);
        
            cur_downbeat_times = (cur_start_time:cur_beat_rate:cur_start_downbeat+time_window);
            cur_halfbeat_times = cur_downbeat_times + cur_beat_rate/2;
            cur_quarterbeat_times = [cur_halfbeat_times + cur_beat_rate/4 cur_downbeat_times + cur_beat_rate/4];
            cur_quarterbeat_times = sort(cur_quarterbeat_times);
            
        
            % Get flux values within this window
            sum_flux_mask = t_bins < cur_downbeat_times(end) & t_bins > cur_downbeat_times(1);
            sum_flux_at_window = sum_flux(sum_flux_mask);
            flux_time_values = t_bins(sum_flux_mask);
            max_flux = max(sum_flux_at_window);
        
            half_flux_mask = t_bins < cur_halfbeat_times(end) & t_bins > cur_halfbeat_times(1);
           
        
            hold on;
            subplot(3,3,ii)
            plot(flux_time_values, sum_flux_at_window, ...
                cur_downbeat_times, ones(length(cur_downbeat_times),1) * max_flux/2, 'x', ...
                cur_halfbeat_times, ones(length(cur_halfbeat_times),1) * max_flux/4, '*', ...
                cur_quarterbeat_times, ones(length(cur_quarterbeat_times),1) * max_flux/8, '^');
            xlim([t_a,t_a+time_window])
            title(sprintf('BPM: %i', cur_bpm))
            legend('Energy Flux', 'Downbeat', 'Half Beat')
        
        end

    end

    beat_locations = struct();
    beat_locations.bpm_candidates = bpm_candidates;
    beat_locations.start_time_for_candidates = start_time_for_candidates;
    beat_locations.cross_correlation_values = cross_correlation_values;
    beat_locations.downbeat_locs = downbeat_locs;
    beat_locations.top_indices = top_indices;

end



















