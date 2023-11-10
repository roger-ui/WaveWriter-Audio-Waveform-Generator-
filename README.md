# Wave Audio File Generator

## Description
This C++ program generates various types of sound waves and writes them to a WAV file. It supports wave types like square, sawtooth, triangle, sine, white noise, brown noise, and chirp. The program is also capable of creating WAV files shorter than 1 second. 

Additional features include the application of low-pass filtering, amplitude modulation, and amplitude envelopes.

## Installation

1. Clone the repository:
   ```
   git clone [repository URL]
   ```
2. Navigate to the cloned directory:
   ```
   cd WaveFileGenerator
   ```
3. Compile the source code:
   ```
   g++ -o WaveFileGenerator WriteWaveFile.cpp
   ```


## Usage

To use the WaveFileGenerator, run the executable from the command line with the desired parameters.

### Command Line Arguments

- `--type` : Type of wave (square, sawtooth, triangle, sine, whitenoise, brownnoise, chirp)
- `--freq` : Frequency of the wave (for tonal waves)
- `--duration` : Duration of the audio in seconds
- `--output` : Output file name
- `--use-filter` : Apply a low pass filter
- `--cutoff` : Cutoff frequency for the low pass filter
- `--amplitude` : Amplitude factor of the wave
- `--start-freq` : Start frequency for chirp waves
- `--end-freq` : End frequency for chirp waves
- `--apply-envelope` : Apply an amplitude envelope to the wave

### Examples

- Generate a 5-second sine wave at 440 Hz:
  ```
  ./WaveFileGenerator.exe --type sine --freq 440 --duration 5 --output my_sine_wave.wav
  ```
- Generate a 10-second chirp from 200 Hz to 2000 Hz:
  ```
  ./WaveFileGenerator.exe --type chirp --start-freq 200 --end-freq 2000 --duration 10 --output chirp_wave.wav
  ```

### Licence 

WaveFileGenerator is licensed under the MIT License. The software as a whole is provided under this license. 

See LICENSE in this repository for details.