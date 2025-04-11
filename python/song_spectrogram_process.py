import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import spectrogram
from scipy.signal.windows import hann
import os

def song_spectrogram_process(filename, nfft, noverlap, min_thresh, show_plots=False):
    """
    Process a .wav file to generate spectrograms for both audio channels.
    
    Parameters:
    ----------
    filename : str
        Path to the .wav file
    nfft : int
        Number of points for FFT
    noverlap : int
        Number of points to overlap between segments
    min_thresh : float
        Minimum threshold for spectrogram
    show_plots : bool, optional
        Whether to display and save plots (default False)
        
    Returns:
    -------
    spectrogram_struct : list of dict
        List containing spectrogram data for each channel
    """
    # Read the audio file
    Fs, y = wavfile.read(filename)
    
    # Check if data is mono or stereo
    if len(y.shape) == 1:
        # Mono - single channel
        y_ch1 = y
        
        # Filter out padded zeros on start/finish
        mask_ch1 = y_ch1 != 0
        y_ch1 = y_ch1[mask_ch1]
        
        # Create a zero array for y_ch2 as a placeholder
        y_ch2 = np.zeros_like(y_ch1)
    else:
        # Stereo - left/right channels
        y_ch1 = y[:, 0]
        y_ch2 = y[:, 1]
        
        # Filter out padded zeros on start/finish
        mask_ch1 = y_ch1 != 0
        y_ch1 = y_ch1[mask_ch1]
        
        mask_ch2 = y_ch2 != 0
        y_ch2 = y_ch2[mask_ch2]
    
    # Create Hanning window
    window = hann(nfft)
    
    # Get the filename without extension for plot titles
    basename = os.path.basename(filename)
    sname = os.path.splitext(basename)[0]
    
    if show_plots:
        plt.figure(figsize=(12, 8))
        
        # Plot channel 1
        plt.subplot(2, 1, 1)
        f_ch1, t_ch1, Sxx_ch1 = spectrogram(y_ch1, fs=Fs, window=window, nperseg=nfft, 
                                         noverlap=noverlap, scaling='density', mode='psd')
        # Apply min threshold
        Sxx_ch1[Sxx_ch1 < min_thresh] = min_thresh
        
        # Convert to dB
        Sxx_ch1_db = 10 * np.log10(Sxx_ch1)
        
        plt.pcolormesh(t_ch1, f_ch1/1000, Sxx_ch1_db, shading='gouraud')
        plt.title(f'CH 1: {sname}', fontweight='bold', fontsize=14)
        plt.ylabel('Frequency [kHz]', fontweight='bold', fontsize=14)
        plt.ylim(0, 22)
        plt.grid(True)
        
        # Plot channel 2
        plt.subplot(2, 1, 2)
        f_ch2, t_ch2, Sxx_ch2 = spectrogram(y_ch2, fs=Fs, window=window, nperseg=nfft, 
                                         noverlap=noverlap, scaling='density', mode='psd')
        # Apply min threshold
        Sxx_ch2[Sxx_ch2 < min_thresh] = min_thresh
        
        # Convert to dB
        Sxx_ch2_db = 10 * np.log10(Sxx_ch2)
        
        plt.pcolormesh(t_ch2, f_ch2/1000, Sxx_ch2_db, shading='gouraud')
        plt.title(f'CH 2: {sname}', fontweight='bold', fontsize=14)
        plt.ylabel('Frequency [kHz]', fontweight='bold', fontsize=14)
        plt.xlabel('Time [sec]', fontweight='bold', fontsize=14)
        plt.ylim(0, 22)
        plt.grid(True)
        
        plt.tight_layout()
        
        # Save the figure
        # Note: Changed the path to a more generic one for platform independence
        # You should modify this to your desired path
        os.makedirs('spectra', exist_ok=True)
        fullfigname = os.path.join('spectra', f'{sname}.png')
        plt.savefig(fullfigname)
        plt.close()
    
    # Calculate spectrograms
    f_ch1, t_ch1, s_ch1 = spectrogram(y_ch1, fs=Fs, window=window, nperseg=nfft, 
                                    noverlap=noverlap, scaling='density', mode='psd')
    
    f_ch2, t_ch2, s_ch2 = spectrogram(y_ch2, fs=Fs, window=window, nperseg=nfft, 
                                    noverlap=noverlap, scaling='density', mode='psd')
    
    # Apply min threshold
    s_ch1[s_ch1 < min_thresh] = min_thresh
    s_ch2[s_ch2 < min_thresh] = min_thresh
    
    # Create spectrogram structure
    spectrogram_struct = [
        {
            'ch_num': 1,
            's': s_ch1,
            'f': f_ch1,
            't': t_ch1
        },
        {
            'ch_num': 2,
            's': s_ch2,
            'f': f_ch2,
            't': t_ch2
        }
    ]
    
    return spectrogram_struct
