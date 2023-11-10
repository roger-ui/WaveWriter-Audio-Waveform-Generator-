// WriteWaveFile.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

constexpr int SAMPLE_RATE = 44100;
//constexpr int NUM_SECONDS = 5;
//constexpr double FREQUENCY = 200.0; // Frequency set to A note (440 Hz)

class LowPassFilter {
private:
    double alpha;
    double prevY = 0.0;

public:
    LowPassFilter(double cutoffFrequency) {
        double RC = 1.0 / (cutoffFrequency * 2 * M_PI);
        double dt = 1.0 / SAMPLE_RATE;
        alpha = dt / (RC + dt);
    }

    short filter(short sample) {
        double x = sample / 32767.0;
        double y = alpha * x + (1.0 - alpha) * prevY;
        prevY = y;
        return static_cast<short>(y * 32767);
    }
};

void writeWavFile(const std::string& filename, const std::vector<short>& samples)
{
    std::ofstream file(filename, std::ios::binary);

    //WAV header 
    file.write("RIFF", 4);
    int fileLength = 36 + samples.size() * sizeof(short);
    file.write(reinterpret_cast<const char*>(&fileLength), 4);
    file.write("WAVE", 4);
    file.write("fmt ", 4);

    int fmtChunkSize = 16;
    short audioFormat = 1;
    short numChannels = 1;
    int byteRate = SAMPLE_RATE * numChannels * sizeof(short);
    short blockAlign = numChannels * sizeof(short);
    short bitsPerSample = 8 * sizeof(short);

    file.write(reinterpret_cast<const char*>(&fmtChunkSize), 4);
    file.write(reinterpret_cast<const char*>(&audioFormat), 2);
    file.write(reinterpret_cast<const char*>(&numChannels), 2);
    file.write(reinterpret_cast<const char*>(&SAMPLE_RATE), 4);
    file.write(reinterpret_cast<const char*>(&byteRate), 4);
    file.write(reinterpret_cast<const char*>(&blockAlign), 2);
    file.write(reinterpret_cast<const char*>(&bitsPerSample), 2);

    file.write("data", 4);
    int dataChunkSize = samples.size() * sizeof(short);
    file.write(reinterpret_cast<const char*>(&dataChunkSize), 4);
    file.write(reinterpret_cast<const char*>(samples.data()), samples.size() * sizeof(short));

}

void applyLowPassFilter(std::vector<short>& samples, double cutoffFrequency) 
{
    
    LowPassFilter filter(cutoffFrequency);
    
    for (int i = 0; i < samples.size(); ++i) 
    {
        samples[i] = filter.filter(samples[i]);
    }
}

void generateSquareWave(std::vector<short>& samples) 
{
    
    bool high = true;

    for (int i = 0; i < samples.size(); ++i) 
    {
        samples[i] = high ? 32767 : -32768;
        high = !high;
    }
}

void generateSoftenedSawtoothWave(std::vector<short>& samples, double frequency) 
{
   
    for (int i = 0; i < samples.size(); ++i) 
    {
        double t = static_cast<double>(i) / SAMPLE_RATE;
        double value = fmod(t * frequency, 1.0);
        value = value * value; // Introduce some non-linearity
        samples[i] = static_cast<short>((value * 2 - 1) * 32767);
    }
}

void generateTriangleWave(std::vector<short>& samples, double frequency) 
{
    double period = SAMPLE_RATE / frequency;
    
    for (int i = 0; i < samples.size(); ++i) 
    {
        double phase = fmod(static_cast<double>(i), period) / period; // Normalize to [0, 1]
        double value;
        if (phase < 0.5) {
            value = 2.0 * phase; // Ramp up
        }
        else {
            value = 2.0 * (1.0 - phase); // Ramp down
        }
        samples[i] = static_cast<short>((value * 2 - 1) * 32767);
    }
}

void generateSineWave(std::vector<short>& samples, double frequency) 
{
    
    for (int i = 0; i < samples.size(); ++i) 
    {
        // Sine wave formula: A * sin(2 * PI * frequency * t + phase)
        double t = static_cast<double>(i) / SAMPLE_RATE; // time in seconds
        double value = sin(2.0 * M_PI * frequency * t);
        samples[i] = static_cast<short>(value * 32767); // A is assumed to be 1 (max amplitude)
    }
}

void generateWhiteNoise(std::vector<short>& samples)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<short> dis(-32768, 32767);

    for (int i = 0; i < samples.size(); ++i) {
        samples[i] = dis(gen);
    }
}


void generateBrownNoise(std::vector<short>& samples) 
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dis(0.0, 1.0);

    double brown = 0.0;
    for (int i = 0; i < samples.size(); ++i) {
        double white = dis(gen);

        // Integrate the white noise to get brown noise
        brown += white;
        
        // Clip the values to prevent overflow and ensure it stays within the 16-bit range
        if (brown < -32768) brown = -32768;
        if (brown > 32767) brown = 32767;

        samples[i] = static_cast<short>(brown);

        // Apply a very simple filter to shape the spectrum for brown noise characteristics
        brown *= 0.998; 
    }
}

void generateChirp(std::vector<short>& samples, double startFrequency, double endFrequency, double durationSeconds) 
{
    double frequencyRange = endFrequency - startFrequency;
    double chirpRate = frequencyRange / durationSeconds; // frequency change per second

    for (size_t i = 0; i < samples.size(); ++i) {
        double t = static_cast<double>(i) / SAMPLE_RATE; // time in seconds
        
        // Linearly change frequency with time
        double instantaneousFrequency = startFrequency + chirpRate * t;
        double value = sin(2.0 * M_PI * instantaneousFrequency * t);
        samples[i] = static_cast<short>(value * 32767); // Amplitude is assumed to be 1 (max amplitude)
    }
}

void modulateAmplitude(std::vector<short>& samples, double amplitudeFactor) {
    for (int i = 0; i < samples.size(); ++i) {
        samples[i] = static_cast<short>(samples[i] * amplitudeFactor);
    }
}

void applyEnvelope(std::vector<short>& samples, double numSeconds) 
{
    int attack_samples = static_cast<int>(SAMPLE_RATE * numSeconds * 0.1); // 10% of the total duration
    int release_samples = static_cast<int>(SAMPLE_RATE * numSeconds * 0.1);

    for (int i = 0; i < attack_samples; ++i) {
        double multiplier = static_cast<double>(i) / attack_samples;
        samples[i] *= multiplier;
    }
    for (int i = samples.size() - release_samples; i < samples.size(); ++i) {
        double multiplier = static_cast<double>(samples.size() - i) / release_samples;
        samples[i] *= multiplier;
    }
}

void printHelp() {
    std::cout << "Usage: WriteWaveFile [OPTION]...\n"
        << "Generate a wave file with specified options.\n\n"
        << "--type        Set the type of wave to generate (square, sawtooth, triangle, sine, whitenoise, brownnoise, chirp).\n"
        << "              Example: --type sine\n"
        << "--freq        Set the frequency of the wave (applies to tonal wave types only).\n"
        << "              Example: --freq 440\n"
        << "--duration    Set the duration of the audio in seconds. Fractional seconds can be specified.\n"
        << "              Example: --duration 2.5\n"
        << "--output      Set the output file name.\n"
        << "              Example: --output my_wave.wav\n"
        << "--use-filter  Apply a low pass filter to the generated wave. (No additional value needed)\n"
        << "              Example: --use-filter\n"
        << "--cutoff      Set the cutoff frequency for the low pass filter.\n"
        << "              Example: --cutoff 1500\n"
        << "--amplitude   Set the amplitude factor (0.0 to 1.0) for the generated wave.\n"
        << "              Example: --amplitude 0.5\n"
        << "--start-freq  Set the start frequency for chirp waves.\n"
        << "              Example: --start-freq 200\n"
        << "--end-freq    Set the end frequency for chirp waves.\n"
        << "              Example: --end-freq 2000\n"
        << "--apply-envelope Apply an amplitude envelope to the wave (attack and release phases).\n"
        << "              Example: --apply-envelope\n"
        << "--help        Display this help and exit.\n\n"
        << "Full Example:\n"
        << "WriteWaveFile --type sine --freq 440 --duration 5 --output my_sine_wave.wav --amplitude 0.5\n"
        << "This will generate a 5-second sine wave with a frequency of 440 Hz, at half the maximum amplitude,\n"
        << "and save it to 'my_sine_wave.wav'.\n\n"
        << "WriteWaveFile --type chirp --start-freq 200 --end-freq 2000 --duration 10 --output chirp_wave.wav\n"
        << "This will generate a 10-second chirp wave starting at 200 Hz and ending at 2000 Hz,\n"
        << "and save it to 'chirp_wave.wav'.\n"
        << "Fractional durations are supported, e.g., --duration 2.5 will generate a 2.5-second audio file.\n";
}



int main(int argc, char* argv[]) 
{
    double numSeconds = 5; // Default duration
    double frequency = 440; // Default frequency
    double cutoffFrequency = 2000.0; // Default cutoff frequency for low pass filter
    bool useLowPassFilter = false; // Flag to apply low pass filter
    bool useEnvelope = false; 
    std::string waveType = "square"; // Default wave type
    std::string outputFilename = "output.wav"; // Default output file name
    std::vector<short> samples;
    double amplitudeFactor = 1.0; // Default amplitude factor
    bool amplitudeModulation = false; // Flag to indicate if amplitude modulation is requested
    double startFrequency = 200.0;
    double endFrequency = 2000.0;


    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--type") == 0 && i + 1 < argc) {
            waveType = argv[++i];
        }
        else if (strcmp(argv[i], "--freq") == 0 && i + 1 < argc) {
            frequency = atof(argv[++i]); // atof can convert a string to a double
        }
        else if (strcmp(argv[i], "--duration") == 0 && i + 1 < argc) {
            numSeconds = atof(argv[++i]); 
        }
        else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            outputFilename = argv[++i];
        }
        else if (strcmp(argv[i], "--use-filter") == 0) {
            useLowPassFilter = true;
        }
        else if (strcmp(argv[i], "--cutoff") == 0 && i + 1 < argc) {
            cutoffFrequency = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--amplitude") == 0 && i + 1 < argc) {
            amplitudeFactor = atof(argv[++i]);
            amplitudeModulation = true; 
        }
        else if (strcmp(argv[i], "--apply-envelope") == 0)
        {
            useEnvelope = true; 
        }
        else if (strcmp(argv[i], "--start-freq") == 0 && i + 1 < argc)
        {
            startFrequency = atof(argv[++i]);   
        }
        else if (strcmp(argv[i], "--end-freq") == 0 && i + 1 < argc)
        {
            endFrequency = atof(argv[++i]); 
        }

        else if (strcmp(argv[i], "--help") == 0) {
			printHelp();
			return 0;
		}
    }

    samples.resize(static_cast<size_t>(ceil(SAMPLE_RATE * numSeconds)));

    // Generate the appropriate wave based on the waveType
    if (waveType == "square") {
        generateSquareWave(samples);
    }
    else if (waveType == "sawtooth") {
        generateSoftenedSawtoothWave(samples, frequency);
    }
    else if (waveType == "triangle") {
        generateTriangleWave(samples, frequency);
    }
    else if (waveType == "sine") {
        generateSineWave(samples, frequency);
    }
    else if (waveType == "whitenoise") {
		generateWhiteNoise(samples);
	}
    else if (waveType == "brownnoise") {
		generateBrownNoise(samples);
	}
    else if (waveType == "chirp")
    {
        generateChirp(samples, startFrequency, endFrequency, numSeconds);
    }
    else {
		std::cout << "Invalid wave type: " << waveType << std::endl;
		return 1;
	}

    // Apply amplitude modulation if requested
    if (amplitudeModulation) {
        modulateAmplitude(samples, amplitudeFactor);
    }
    
    // Apply envelope if requested
    if (useEnvelope)
    {
        applyEnvelope(samples, numSeconds);
    }
    
    // Apply low pass filter if requested
    if (useLowPassFilter) 
    {
        applyLowPassFilter(samples, cutoffFrequency);
    }

    writeWavFile(outputFilename, samples);

    std::cout << "WAV file " << outputFilename << " has been generated!" << std::endl;

    return 0;
}


