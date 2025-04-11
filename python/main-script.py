import os
import glob
import numpy as np
from song_spectrogram_process import song_spectrogram_process
from flux_calculator import flux_calculator

# User-defined variables
# Get the files of interest
songloc = os.path.join('/Users/carleen/Documents/misc_code/audio-whisk/audio_files/sparhawk_white-roses', '*.wav')
ii = 4  # Index of song to use (Python uses 0-based indexing, so 5-1)
show_plots = True

nfft = 512  # Gives time resolution of ~ 10ms
noverlap = 0  # No overlap needed
min_thresh = -60

f_min = 100  # Hz
f_max = 10e3  # Hz

# File I/O
songs = sorted(glob.glob(songloc))
if not songs:
    print(f"No .wav files found at: {songloc}")
    exit(1)
    
filename = songs[ii]

# Frequency-domain processing
spectrogram_struct = song_spectrogram_process(filename, nfft, noverlap, 
                                           min_thresh, show_plots)

# Energy Flux Calculation
# Calculate the Energy Flux (Laroche, J. "Efficient Tempo and Beat Tracking
# in Audio Recordings")
t_bins = spectrogram_struct[0]['t']
f_bins = spectrogram_struct[0]['f']
spectrogram_struct, sum_flux = flux_calculator(spectrogram_struct, 
                                            t_bins, f_bins, f_min, f_max, show_plots)

print(f"Processing complete for file: {os.path.basename(filename)}")
