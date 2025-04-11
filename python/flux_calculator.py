import numpy as np
import matplotlib.pyplot as plt

def flux_calculator(spectrogram_struct, t_bins, f_bins, f_min, f_max, show_plots=False):
    """
    Calculate spectral flux from spectrogram data across frequency bins.
    
    Parameters:
    ----------
    spectrogram_struct : list of dict
        List containing spectrogram data for each channel
    t_bins : array
        Time bins from spectrogram
    f_bins : array
        Frequency bins from spectrogram
    f_min : float
        Minimum frequency to include in flux calculation
    f_max : float
        Maximum frequency to include in flux calculation
    show_plots : bool, optional
        Whether to display and save plots (default False)
        
    Returns:
    -------
    spectrogram_struct : list of dict
        Updated spectrogram data with flux for each channel
    sum_flux : array
        Sum of flux across all channels
    """
    num_ch = len(spectrogram_struct)
    num_t_bins = len(t_bins)
    
    # Limit the frequency summed to f_min, f_max
    f_mask = (f_bins > f_min) & (f_bins < f_max)
    
    # Initialize array for sum_flux
    sum_flux = np.zeros(num_t_bins)
    
    for ii_ch in range(num_ch):
        # Get magnitude of complex spectrogram
        mag_ch = np.abs(spectrogram_struct[ii_ch]['s'])
        
        # Apply arcsin-like power transformation
        arcsin_ch = np.power(mag_ch, 1/2)
        
        frame_flux_arr = np.zeros(num_t_bins)
        
        # Calculate flux for each time bin (starting from second bin)
        for ii_t in range(1, num_t_bins):
            f_cur = arcsin_ch[f_mask, ii_t]
            f_prev = arcsin_ch[f_mask, ii_t-1]
            
            # Get the flux at each freq bin
            f_flux = f_cur - f_prev
            
            # Sum over the frequency bins
            frame_flux_arr[ii_t] = np.sum(f_flux)
        
        # Set negative flux values to zero
        flux_mask = frame_flux_arr < 0
        frame_flux_arr[flux_mask] = 0
        
        # Store flux in structure and update sum
        spectrogram_struct[ii_ch]['flux'] = frame_flux_arr
        sum_flux = spectrogram_struct[ii_ch]['flux'] + sum_flux
    
    # Plot figures if wanted
    if show_plots:
        plt.figure(figsize=(12, 8))
        plt.grid(True)
        plt.title('Energy Flux', fontweight='bold', fontsize=14)
        plt.xlabel('Time (s)', fontweight='bold', fontsize=14)
        plt.ylabel('Energy Flux, E(n)', fontweight='bold', fontsize=14)
        plt.plot(t_bins, sum_flux)
        plt.tight_layout()
        plt.show()
    
    return spectrogram_struct, sum_flux
