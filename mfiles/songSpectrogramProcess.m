function [spectrogram_struct] = songSpectrogramProcess(filename, nfft, noverlap,...
    min_thresh, show_plots)

        % y is our sound data, Fs is our sample rate
        [y, Fs] = audioread(filename);
    
        samplenum = length(y);
        try
            % Get left/right channels
            y_ch1 = y(:,1);
            y_ch2 = y(:,2);
            
            % Filter out padded zeros on start/finish
            mask_ch1 = y_ch1 ~=0;
            y_ch1 = y_ch1(mask_ch1);
        
            mask_ch2 = y_ch2 ~=0;
            y_ch2 = y_ch2(mask_ch2);
        catch
            % Mono, single channel only
            y_ch1 = y(:,1);
            
            % Filter out padded zeros on start/finish
            mask_ch1 = y_ch1 ~=0;
            y_ch1 = y_ch1(mask_ch1);
        end
               
        %% Keep it simple with the spectrogram
        window = hanning(nfft);

        [pathstr, sname, ext] = fileparts(filename);
        
        if show_plots
            % Plot the spectrogram and save it to file. 
            f1=figure; hold on; grid on;
            f1.Position = [680 254 1037 625];

            subplot(2,1,1)
            spectrogram(y_ch1, window, noverlap, nfft, Fs, 'onesided','minthreshold', min_thresh, ...
                'power', 'yaxis');
            title(sprintf('CH 1: %s'), sname)
            ylim([0,22]); % This is in kHz
            ax=gca;
            ax.FontWeight='bold';
            ax.FontSize=14;
        
            subplot(2,1,2)
            spectrogram(y_ch2, window, noverlap, nfft, Fs, 'onesided','minthreshold', min_thresh, ...
            'power', 'yaxis');
            title(sprintf('CH 2: %s'), sname)
            ylim([0,22]); % This is in kHz
            ax=gca;
            ax.FontWeight='bold';
            ax.FontSize=14;


            fullfigname = sprintf('/Users/boyercm1/Documents/mfiles/spectra/%s', sname);
            savefig(fullfigname);
            close all;
        end



        [s_ch1, f_ch1, t_ch1] = spectrogram(y_ch1, window, noverlap, nfft, Fs, 'onesided','minthreshold', min_thresh, ...
            'power', 'yaxis');
        [s_ch2, f_ch2, t_ch2] = spectrogram(y_ch2, window, noverlap, nfft, Fs, 'onesided','minthreshold', min_thresh, ...
            'power', 'yaxis');
   


        %% Populate spectogram_struct with calculated values for each

        spectrogram_struct = struct();
        spectrogram_struct(1).ch_num = 1;
        spectrogram_struct(1).s = s_ch1;
        spectrogram_struct(1).f = f_ch1;
        spectrogram_struct(1).t = t_ch1;


        spectrogram_struct(2).ch_num = 2;
        spectrogram_struct(2).s = s_ch2;
        spectrogram_struct(2).f = f_ch2;
        spectrogram_struct(2).t = t_ch2;

    
    end
 
    
    
    
    
    
    

