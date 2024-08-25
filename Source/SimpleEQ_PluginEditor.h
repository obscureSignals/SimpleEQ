
#pragma once

#include "../../Gui/PresetPanel.h"
#include "../../PluginProcessor.h"
#include "SimpleEQ_PluginProcessor.h"
#include <gin/gin.h>
#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>
#include <juce_gui_basics/juce_gui_basics.h>
#include <sound_meter/sound_meter.h>

// TODO: numAves should be a constructor argument?
template <typename BlockType>
struct FFTDataGenerator
{
    /**
     Produces the FFT data from an audio buffer
     audioData is incoming time domain data
     */

    explicit FFTDataGenerator (const float negativeInfinity) : fftOrder (0), fftSize (0), negativeInfinity (negativeInfinity)
    {
    }
    void produceFFTDataForRendering (const juce::AudioBuffer<float>& audioData, float newFFTOrder)
    {
        if (newFFTOrder != fftOrder)
        {
            changeOrder (newFFTOrder);
        }

        // Copy contents of audioData into fftData
        auto* readIndex = audioData.getReadPointer (0);
        std::copy (readIndex, readIndex + fftSize, fftData.begin());

        // Window time domain data
        window->multiplyWithWindowingTable (fftData.data(), fftSize);

        // Perform FFT for just magnitude values
        // First half of fftData is now magnitude spectrum
        forwardFFT->performFrequencyOnlyForwardTransform (fftData.data());

        // Normalize the FFT values
        int numBins = (int) fftSize / 2;
        for (int i = 0; i < numBins; ++i)
        {
            auto v = fftData[i];

            if (!std::isinf (v) && !std::isnan (v))
            {
                v /= static_cast<float> (numBins);
            }
            else
            {
                v = 0.f; // Set infs & NaNs to zero
            }
            fftData[i] = v;
        }

        // Convert magnitude values to dB
        for (int i = 0; i < numBins; ++i)
        {
            fftData[i] = juce::Decibels::gainToDecibels (fftData[i], negativeInfinity);
        }

        // Add fftData to fftDataVector
        std::copy (fftData.begin(), fftData.begin() + fftSize - 1, fftDataVector[aveCounter].begin());

        aveCounter++;
        if (aveCounter % numAves == 0)
        {
            aveCounter = 0;
        }

        // Find average of fftDataVector and copy into fftDataAve
        std::fill (fftDataAve.begin(), fftDataAve.end(), 0.f);
        for (int aveIdx = 0; aveIdx < numAves; aveIdx++)
        {
            for (int dataIdx = 0; dataIdx < fftSize; dataIdx++)
            {
                fftDataAve[dataIdx] = fftDataAve[dataIdx] + fftDataVector[aveIdx][dataIdx] / static_cast<float> (numAves);
            }
        }

        // Push averaged data to FIFO
        fftDataFifo.push (fftDataAve);
    }

    void changeOrder (int newOrder)
    {
        // When you change order, recreate the window, forwardFFT, fifo, fftData
        // Also reset the fifoIndex
        // Things that need recreating should be created on the heap via std::make_unique<>

        fftOrder = newOrder;
        fftSize = getFFTSize();

        forwardFFT = std::make_unique<juce::dsp::FFT> (fftOrder); // Make FFT object
        window = std::make_unique<juce::dsp::WindowingFunction<float>> (fftSize, juce::dsp::WindowingFunction<float>::hann);

        fftData.clear();
        fftData.resize (fftSize * 2, 0); // FFT input/output array must be twice as long as the incoming time domain data

        fftDataVector.resize (numAves);
        for (int aveIdx = 0; aveIdx < numAves; aveIdx++)
        {
            fftDataVector[aveIdx].clear();
            fftDataVector[aveIdx].resize (fftSize, negativeInfinity);
        }
        fftDataAve.clear();
        fftDataAve.resize (fftSize, negativeInfinity);

        fftDataFifo.prepare (fftData.size());
    }
    [[nodiscard]] int getFFTSize() const { return 1 << fftOrder; } // FFT size is 2^fftOrder
    [[nodiscard]] int getNumAvailableFFTDataBlocks() const { return fftDataFifo.getNumAvailableForReading(); }
    bool getFFTData (BlockType& fftData) { return fftDataFifo.pull (fftData); }

private:
    int fftOrder;
    int fftSize; // Size of time domain data that FFT will work on is 2^fftOrder
    BlockType fftData; // Raw, preaveraged FFT data
    std::vector<BlockType> fftDataVector; // Array of raw, preaveraged FFT data blocks
    BlockType fftDataAve; // Averaged FFT data that is returned

    const int numAves = 4; // Moving average in freq domain, each frame is average of last numAves frames - this could be weighted?
    int aveCounter = 0;
    std::unique_ptr<juce::dsp::FFT> forwardFFT;
    std::unique_ptr<juce::dsp::WindowingFunction<float>> window;

    Fifo<BlockType> fftDataFifo;

    const float negativeInfinity; // Lowest value in dB  - zeros (or negative values) in the magnitude response will be converted to negativeInfinity dB
};

//======================================================

// AnalyzerPathGenerator creates a juce::Path for the spectrum analyzer trace
template <typename PathType>
struct AnalyzerPathGenerator
{
    /*
     Converts 'renderData' into a juce::Path
     */

    explicit AnalyzerPathGenerator (const float negativeInfinity) : negativeInfinity (negativeInfinity)
    {
    }

    //    renderData is a magnitude spectrum in dB
    //    fftBounds is a rectangle representing the bounds of the path
    //       - left side is 0Hz
    //       - right side is 20kHz
    //       - top is 0dB
    //       - bottom (+ 6) is negative infintity dB
    //    binWidth is width of FFT bins in Hz = sampleRate/fftSize
    void generatePath (const std::vector<float>& renderData,
        juce::Rectangle<float> fftBounds,
        int fftSize,
        float binWidth)
    {
        auto top = fftBounds.getY();
        auto bottom = fftBounds.getHeight();
        auto width = fftBounds.getWidth();

        int numBins = ceil (20000.f / binWidth); // Only build the path up to 20kHz

        PathType pSpline; // Spline path

        // These point arrays hold coordinates in pixel space.
        // There are two because gin spline wants gin::Point<double>, but path wants juce::Point<float>
        juce::Array<juce::Point<float>> points;
        juce::Array<gin::Point<double>> dpoints;
        pSpline.preallocateSpace (3 * static_cast<int> (fftBounds.getWidth()));

        // create horizontal map from dB value to pixel space
        const auto negInf = negativeInfinity;
        auto map = [bottom, top, negInf] (float v) {
            return juce::jmap (v, negInf, 0.f, (bottom + 6.f), top);
        };

        // Height of bin
        float prevY = NAN; // Height of previous bin

        // A 'run' is when more than one frequency is mapping to the same pixel in the horizontal axis
        // This will happen at higher frequencies where the bins are close together because of the log scale
        // A typical case would be > 3k frequencies and < 1k pixels
        // We run a 1-point delay at all times when assigning coordinates to point arrays in order to catch runs
        // If we are in a run, we wait until the end of the run to assign the value to the point array

        float prevBinX = NAN; // Previous bin's horizontal coordinate
        int runIdx = 0; // How many bins in the run so far
        bool runFlag = false; // Are we in a run
        bool lastFlag = false; // Last bin
        float minRunY = bottom + 6; // Minimum of values in run

        bool firstFlag = true; // We don't add the first prevY to our point array becuase it doesn't exist yet

        // Fill point arrays
        // Start at idx 1, because 0 is 0Hz which maps to -inf because of log scale
        for (int binNum = 1; binNum <= numBins; binNum++)
        {
            const float y = map (renderData[binNum]); // Set horizontal coordinate for current value

            if (binNum == numBins)
            {
                lastFlag = true;
            }

            if (!std::isnan (y) && !std::isinf (y)) // There should not be NaNs or infs
            {
                auto binFreq = binNum * binWidth; // Frequency of bin in Hz
                auto normalizedBinX = juce::mapFromLog10 (binFreq, 20.f, 20000.f); // Normalized value (0-1) for horizontal coordinate
                float binX = std::floor (normalizedBinX * width); // Horizontal coordinate in pixels
                if (binX != prevBinX)
                { // If the horizontal coordinate is not the same as it was for the last bin
                    if (runFlag or lastFlag)
                    { // If run flag or the last flag are tripped, this must be the end of a run
                        prevY = minRunY;

                        // Reset Run
                        minRunY = bottom + 6.f;
                        runIdx = 0;
                        runFlag = false;
                    }
                    if (!firstFlag)
                    { // Add the previous points to our point arrays
                        points.add ({ prevBinX, prevY });
                        dpoints.add ({ static_cast<double> (prevBinX), static_cast<double> (prevY) });
                    }
                }
                else
                { // If the horizontal coordinate is same as it was for the last bin
                    if (!runFlag)
                    { // If the run flag was not already tripped, this is a new run
                        minRunY = juce::jmin (prevY, y); // Find min of this y and the last y - smaller values are higher up in gui
                        runFlag = true; // We are running
                    }
                    else
                    { // If the run flag was already tripped, we are already in a run
                        minRunY = juce::jmin (minRunY, y); // Keep track of min value - we only want to retain the highest value for pixel
                        runIdx++;
                    }
                }

                // Assign current value as prev values
                prevBinX = binX;
                prevY = y;

                if (binNum == 1)
                {
                    firstFlag = false;
                }
            }
        }

        // Add last bin
        points.add ({ prevBinX, prevY });
        dpoints.add ({ static_cast<double> (prevBinX), static_cast<double> (prevY) });

        // Spline interpolate points
        gin::Spline spline (dpoints);
        pSpline.startNewSubPath (points.getFirst().toFloat());
        for (auto x1 = static_cast<int> (points.getFirst().getX()); x1 < static_cast<int> (points.getLast().getX()); x1++)
        {
            auto y1 = spline.interpolate (x1);
            pSpline.lineTo (x1, y1);
        }

        // Enclose the path so we can fill it
        pSpline.lineTo (points.getLast().getX(), bottom + 10.f);
        pSpline.lineTo (points.getFirst().getX(), bottom + 10.f);
        pSpline.lineTo (points.getFirst().toFloat());

        // Push path to path Fifo
        splinePathFifo.push (pSpline);
    }

    [[nodiscard]] int getNumPathsAvailable() const
    {
        return splinePathFifo.getNumAvailableForReading();
    }

    bool getPath (PathType& path)
    {
        return splinePathFifo.pull (path);
    }

private:
    Fifo<PathType> splinePathFifo;
    const float negativeInfinity;
};

//=================================================

struct LookAndFeel : juce::LookAndFeel_V4
{
    void drawRotarySlider (juce::Graphics&,
        int x,
        int y,
        int width,
        int height,
        float sliderPosProportional,
        float rotaryStartAngle,
        float rotaryEndAngle,
        juce::Slider&) override;

private:
    Gui::colorPalette colors;
};

//==============================================

struct RotarySliderWithLabels : juce::Slider
{
    explicit RotarySliderWithLabels (const juce::String& unitSuffix) : juce::Slider (juce::Slider::SliderStyle::Rotary,
                                                                           juce::Slider::TextEntryBoxPosition::TextBoxBelow)
    {
        setLookAndFeel (&lnf);
        setTextValueSuffix (unitSuffix); // set unit suffixes for value display
        setColour (juce::Slider::textBoxOutlineColourId, juce::Colours::transparentWhite); // make border of value display inivisible
        setColour (juce::Slider::textBoxTextColourId, colors.textColor);
        setColour (juce::Slider::thumbColourId, colors.cloud.withBrightness (0.9f));
        // setSliderSnapsToMousePosition(false); // not available for rotary sliders :(
    }

    ~RotarySliderWithLabels() override
    {
        setLookAndFeel (nullptr);
    }

    // Tic labels and postions
    struct LabelPos
    {
        float pos;
        juce::String label;
    };
    juce::Array<LabelPos> labels;

    juce::String title;

    void paint (juce::Graphics& g) override;
    juce::Rectangle<int> getSliderBounds() const;
    static int getTextHeight() { return 14; }

private:
    Gui::colorPalette colors;
    LookAndFeel lnf;
};

// ===================================================================

struct PathProducer
{
    PathProducer (StereoSummingSampleFifo<PlayBackEQAudioProcessor::BlockType>& sssf, PlayBackEQAudioProcessor& p, const float negativeInfinity) : Fifo (&sssf), audioProcessor (p), FFTDataGenerator (negativeInfinity), pathGenerator (negativeInfinity), negativeInfinity (negativeInfinity)
    {
        double sampleRate = audioProcessor.getSampleRate();
        // Set up the FFTDataGenerator to operate on at ~ 1/6 of a second of time domain data at a time
        int fftOrder = ceil (log2 (sampleRate / 6));
        FFTDataGenerator.changeOrder (fftOrder);
        monoBuffer.setSize (1, FFTDataGenerator.getFFTSize());
    }

    void process (juce::Rectangle<float> fftBounds, double sampleRate);
    juce::Path getPath() { return FFTPath; }

private:
    // This will be spectrumAnalyzerSampleFifo which is owned by the processor and pushed to in processBlock
    StereoSummingSampleFifo<PlayBackEQAudioProcessor::BlockType>* Fifo;

    PlayBackEQAudioProcessor& audioProcessor;

    juce::AudioBuffer<float> monoBuffer; // Will hold time domain data from Fifo to pass to FFTDataGenerator

    FFTDataGenerator<std::vector<float>> FFTDataGenerator;

    AnalyzerPathGenerator<juce::Path> pathGenerator;

    const float negativeInfinity = -80; // Lowest value in dB  - zeros (or negative values) in the magnitude response will be converted to negativeInfinity dB

    juce::Path FFTPath;

    std::vector<float> fftData;

    int fftCount = 0;
};

// ===================================================================

class SpectrumDisplay final : public juce::Component
{
public:
    explicit SpectrumDisplay (PlayBackEQAudioProcessor&);
    ~SpectrumDisplay() override;
    void paint (juce::Graphics& g) override;
    void pathProducerProcess (juce::Rectangle<float> fftBounds);
    void resized() override;
    void updateResponseCurve();
    void updateFiltersGUI();
    void setDisplayColor (juce::Colour newDisplayColor)
    {
        displayColor = newDisplayColor;
    }

private:
    Gui::colorPalette colors;
    juce::Colour displayColor = colors.blue;
    PlayBackEQAudioProcessor& audioProcessor;
    PathProducer leftPathProducer;
    juce::dsp::ProcessorChain<Filter> leftBassChain, leftTrebleChain, leftRIAAChain;
    juce::dsp::ProcessorChain<Filter, Filter, Filter> leftHPChain;
    juce::dsp::ProcessorChain<Filter, Filter> leftLPChain;
    juce::Path responseCurve;
};

//=============================================================

struct ResponseCurveComponent : juce::Component,
                                juce::AudioProcessorParameter::Listener,
                                juce::Timer
{
    explicit ResponseCurveComponent (PlayBackEQAudioProcessor&);
    ~ResponseCurveComponent() override;

    void parameterValueChanged (int parameterIndex, float newValue) override;

    void parameterGestureChanged (int parameterIndex, bool gestureIsStarting) override {}

    void timerCallback() override;

    void paint (juce::Graphics& g) override;

    void paintOverChildren (juce::Graphics& g) override;

    void mouseMove (const juce::MouseEvent& event) override;

    void resized() override;

    void updateReponse()
    {
        spectrumDisplay.resized();
        repaint();
    }

    void setSpectrumDisplayColor (juce::Colour newSpectrumDisplayColor)
    {
        spectrumDisplay.setDisplayColor (newSpectrumDisplayColor);
    }

    void mouseEnter (const juce::MouseEvent& event) override;
    void mouseExit (const juce::MouseEvent& event) override;

private:
    Gui::colorPalette colors;

    PlayBackEQAudioProcessor& audioProcessor;

    SpectrumDisplay spectrumDisplay;

    juce::Atomic<bool> parametersChanged { true };

    void drawBackgroundGrid (juce::Graphics& g) const;
    void drawTextLabels (juce::Graphics& g) const;

    static std::vector<float> getFrequencies();
    static std::vector<float> getGains();
    static std::vector<float> getXs (const std::vector<float>& freqs, float left, float width);

    juce::Rectangle<int> getAnalysisArea() const;

    bool isMouseOverAnalysisArea {};
    float freq {}, mag {};
    int mouseX {}, mouseY {};

    juce::Component coordinateComponent; // Offset so that coordinate display corrosponds to very point of mouse

    // For overlaying multiple paths with some weighted transparency - might be useful at some point
    //    std::vector<juce::Path> FFTpaths;
    //    int numPaths = 10;
    //    int pathCount = 0;
};