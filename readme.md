# Frequency Identification using MUSIC and ESPRIT Algorithms

This project provides a set of MATLAB functions for high-resolution frequency estimation of multiple narrow-band disturbances using Digital Signal Processing (DSP) techniques. The core algorithms implemented are MUSIC (Multiple Signal Classification) and ESPRIT (Estimation of Signal Parameters via Rotational Invariance Techniques).

These methods offer significantly more accurate frequency estimates from smaller data samples compared to traditional Fast Fourier Transform (FFT) based methods, especially for signals buried in white noise.

## Getting Started

### Prerequisites

- MATLAB

### Usage

The main script to run the frequency identification is `music_freq_ID.m`. It uses the sample data provided in `y1.mat`.

To run the analysis, simply execute the main script in your MATLAB environment:

```matlab
music_freq_ID
```

This will process the sample data and identify three narrow-band frequencies using the implemented algorithms.

You should be able to identify three narrow band frequencies using the algorithms.

## File Descriptions

- **`music_freq_ID.m`**: The main script that demonstrates the use of the frequency identification algorithms with the sample data.
- **`m_music.m`**: An implementation of the MUSIC algorithm for frequency estimation.
- **`m_esprit.m`**: An implementation of the ESPRIT algorithm for frequency estimation.
- **`Denoising.m`**: A function for denoising the input signal. For more advanced cases, pre-filtering using MATLAB's Signal Processing Toolbox is recommended.
- **`AF_with_denoising.m`**: Calculates the autocorrelation function with denoising.
- **`specCal.m`**: Calculates the signal spectrum.
- **`convm.m`**: Generates a convolution matrix.
- **`extrema.m`**: A utility function to find the local extrema in a data set.
- **`y1.mat`**: A sample dataset containing a signal with three narrow-band frequencies.
- **`m_AF.m`**: Calculates the autocorrelation function.

## Applicability
The algorithms here can be well applied to handle narrow band signals buried in white noise. For other cases, pre-filtering (denoising) is important. I suggest using MATLAB's signal processing toolbox to design low-pass or band-pass filters for that purpose.

## Author

Xu Chen
