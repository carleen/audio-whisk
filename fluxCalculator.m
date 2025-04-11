function [spectrogram_struct, sum_flux] = fluxCalculator(spectrogram_struct, ...
    t_bins, f_bins, f_min, f_max, show_plots)
    % Input: 
    %   -spectrogram_struct: ch_num, s, f, t

    num_ch = length(spectrogram_struct);
    num_t_bins = length(t_bins);
    
    % Limit the frequency summed to f_min, f_max
    f_mask = f_bins > f_min & f_bins < f_max;

    % Initailize array for sum_flux
    sum_flux = zeros(num_t_bins, 1);

    for ii_ch=1:num_ch
        mag_ch = abs(spectrogram_struct(ii_ch).s);
        arcsin_ch = mag_ch.^(1/2);

        frame_flux_arr = zeros(num_t_bins, 1);
        
        for ii_t = 2:num_t_bins
            f_cur = arcsin_ch(f_mask, ii_t);
            f_prev = arcsin_ch(f_mask, ii_t-1);
        
            % Get the flux at each freq bin; mask out based on min/max freq
            f_flux = f_cur - f_prev;
        
            % Sum over the frequency bins
            frame_flux_arr(ii_t) = sum(f_flux);

        end

        flux_mask = frame_flux_arr < 0;
        frame_flux_arr(flux_mask) = 0;

        spectrogram_struct(ii_ch).flux = frame_flux_arr;
        sum_flux = spectrogram_struct(ii_ch).flux + sum_flux;

    end

    % Plot figures if wanted
    if show_plots
        f1=figure; hold on; grid on;
        f1.Position = [680 254 1037 625];
        ax=gca;
        ax.FontWeight='bold';
        ax.FontSize=14;
        title('Energy Flux')
        xlabel('Time (s)')
        ylabel('Energy Flux, E(n)')
        plot(t_bins, sum_flux)
    end
   
end
