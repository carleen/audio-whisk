clear; clc; close all;

%% User-defined variables
% Get the files of interests
songloc = '/Users/boyercm1/Documents/misc/tempo_processing/audio_files/sparhawk_white-roses/*.wav';
ii = 5; % Index of song to use
show_plots = 1;

nfft = 512; % Gives time resolution of ~ 10ms
noverlap = 0; % No overlap needed
min_thresh = -60;

f_min = 100; % Hz
f_max = 10e3; % Hz


%% File I/O
songs = dir(songloc);
filename = fullfile(songs(ii).folder, songs(ii).name);


%% Frequency-domain processing
[spectrogram_struct] = songSpectrogramProcess(filename, nfft, noverlap, ...
    min_thresh, show_plots);

%% Energy Flux Calculation
% Calculate the Energy Flux (Laroche, J. "Efficent Tempo and Beat Tracking
% in Audio Recordings")
t_bins = spectrogram_struct(1).t;
f_bins = spectrogram_struct(1).f;
[spectrogram_struct, sum_flux] = fluxCalculator(spectrogram_struct, ...
    t_bins, f_bins, f_min, f_max, show_plots);



