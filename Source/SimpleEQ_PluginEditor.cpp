
#include "SimpleEQ_PluginEditor.h"
#include "../../Gui/GetMaxStringWidth.h"

//===========================================// Spectrum Display

SpectrumDisplay::SpectrumDisplay (PlayBackEQAudioProcessor& p) : audioProcessor (p)
// leftPathProducer (audioProcessor.spectrumAnalyzerSampleFifo, p, -80)
{
    // Create log-scale frequency array from 20Hz to 20kHz for frequency response plot
    frequencies.resize (NUM_FREQUENCIES);
    for (size_t i = 0; i < frequencies.size(); ++i)
    {
        frequencies[i] = juce::mapToLog10 (static_cast<double> (i) / static_cast<double> (frequencies.size()), 20.0, 20000.0);
    }

    overallMagnitudeResponse.resize (frequencies.size());
    std::ranges::fill (overallMagnitudeResponse, 1.0);

    for (auto& response : filterMagnitudeResponses)
    {
        response.fill (1.0);
    }

    // Get inverse RIAA filter coefficients
    const auto RIAAcoeffs = getRIAACoeffs (audioProcessor.getUpSampleRate());
    // Convert array to coefficients pointer
    coeffs.riaaCoeffsGui = new juce::dsp::IIR::Coefficients<float> (Coefficients (RIAAcoeffs));

    // Set RIAA magnitude response
    coeffs.riaaCoeffsGui->getMagnitudeForFrequencyArray (frequencies.data(),
        filterMagnitudeResponses[static_cast<size_t> (filters::RIAA)].data(),
        frequencies.size(),
        audioProcessor.getUpSampleRate());

    setInterceptsMouseClicks (false, false); // Passes mouse events through to component behind it
}

SpectrumDisplay::~SpectrumDisplay() = default;

// void SpectrumDisplay::pathProducerProcess (const juce::Rectangle<float> fftBounds)
// {
//     leftPathProducer.process (fftBounds, audioProcessor.getSampleRate());
// }

void SpectrumDisplay::paint (juce::Graphics& g)
{
    // TRACE_COMPONENT();

    using namespace juce;

    // Get FFT path, fill  it and outline it
    // const auto FFTpath = leftPathProducer.getPath();
    // g.setColour (displayColor);
    // g.setOpacity (0.6);
    // g.fillPath (FFTpath);
    // g.strokePath (FFTpath, PathStrokeType (1.f));

    audioProcessor.createAnalyserPlot (analyserPath, getLocalBounds(), 20.0f);
    g.setColour (displayColor);
    g.setOpacity (0.6);
    g.fillPath (analyserPath);
    g.strokePath (analyserPath, PathStrokeType (1.f));

    // Draw response curve
    g.setColour (Colours::white.withBrightness (0.9));

    g.strokePath (magnitudeResponse, PathStrokeType (2.f));
}

void SpectrumDisplay::updateCoeffsGUI() // Updates the filter coefficients owned by the editor with the current control values
{
    const double sampleRate = audioProcessor.getSampleRate();

    // Calculate bass filter coefficients
    const auto bassCoeffs = getBassCoeffs (audioProcessor.getCurrentSettings(), sampleRate);
    // Convert array to coefficients pointer
    coeffs.bassCoeffsGui = new juce::dsp::IIR::Coefficients<float> (Coefficients (bassCoeffs));

    // Calculate treble filter coefficients
    const auto trebleCoeffs = getTrebleAttenCoeffs (audioProcessor.getCurrentSettings(), audioProcessor.getUpSampleRate());
    // Convert array to coefficients pointer
    coeffs.trebleCoeffsGui = new juce::dsp::IIR::Coefficients<float> (Coefficients (trebleCoeffs));

    // Get coefficients for high-pass filter
    coeffs.hpCoeffsGui = getHighPassCoeffs (audioProcessor.getCurrentSettings(), sampleRate);

    // Get low pass coefficients
    coeffs.lpCoeffsGui = getLowPassCoeffs (audioProcessor.getCurrentSettings(), audioProcessor.getSampleRate());

    // Set coefficients for special curve filters
    juce::ReferenceCountedArray<juce::dsp::IIR::Coefficients<float>> specialCurveCoeffs;
    if (audioProcessor.getCurrentSettings().specialCurve == SpecialCurve::BBC)
    {
        coeffs.specialCoeffsGui = getBBCcoeffs (audioProcessor.getSampleRate());
    }
    else if (audioProcessor.getCurrentSettings().specialCurve == SpecialCurve::BBC_XTR_311)
    {
        coeffs.specialCoeffsGui = getXTRcoeffs (audioProcessor.getSampleRate());
    }
}

void SpectrumDisplay::createFrequencyPlot (juce::Path& p, const std::vector<double>& mags, const juce::Rectangle<int> bounds, const float maxGainDB) const
{
    const auto y1 = juce::jmap (-1.f * juce::Decibels::gainToDecibels (static_cast<float> (mags[0])), -1 * maxGainDB, maxGainDB, static_cast<float> (bounds.getY()), static_cast<float> (bounds.getBottom()));

    p.startNewSubPath (static_cast<float> (bounds.getX()), y1);
    const auto xFactor = static_cast<double> (bounds.getWidth()) / frequencies.size();
    for (size_t i = 1; i < frequencies.size(); ++i)
    {
        const auto y = juce::jmap (-1.f * juce::Decibels::gainToDecibels (static_cast<float> (mags[i])), -1 * maxGainDB, maxGainDB, static_cast<float> (bounds.getY()), static_cast<float> (bounds.getBottom()));
        p.lineTo (static_cast<float> (bounds.getX() + static_cast<double> (i) * xFactor), y);
    }
}

void SpectrumDisplay::updateOverallMagResponse()
{
    // Set overallMagnitudeResponse to all ones
    std::ranges::fill (overallMagnitudeResponse, 1.0);

    // Multiply overallMagnitudeResponse by response of each filter that is not bypassed
    if (audioProcessor.getCurrentSettings().bassBoostEnable.load())
        juce::FloatVectorOperations::multiply (overallMagnitudeResponse.data(), filterMagnitudeResponses[static_cast<size_t> (filters::bass)].data(), NUM_FREQUENCIES);
    if (audioProcessor.getCurrentSettings().trebleCutEnable.load())
        juce::FloatVectorOperations::multiply (overallMagnitudeResponse.data(), filterMagnitudeResponses[static_cast<size_t> (filters::treble)].data(), NUM_FREQUENCIES);
    if (audioProcessor.getCurrentSettings().RIAAEnable.load())
        juce::FloatVectorOperations::multiply (overallMagnitudeResponse.data(), filterMagnitudeResponses[static_cast<size_t> (filters::RIAA)].data(), NUM_FREQUENCIES);
    if (audioProcessor.getCurrentSettings().specialCurveEnable.load())
        juce::FloatVectorOperations::multiply (overallMagnitudeResponse.data(), filterMagnitudeResponses[static_cast<size_t> (filters::special)].data(), NUM_FREQUENCIES);
    if (audioProcessor.getCurrentSettings().highPassSlope.load() != Slope::Off)
        juce::FloatVectorOperations::multiply (overallMagnitudeResponse.data(), filterMagnitudeResponses[static_cast<size_t> (filters::highPass)].data(), NUM_FREQUENCIES);
    if (audioProcessor.getCurrentSettings().lowPassSlope.load() != Slope::Off)
        juce::FloatVectorOperations::multiply (overallMagnitudeResponse.data(), filterMagnitudeResponses[static_cast<size_t> (filters::lowPass)].data(), NUM_FREQUENCIES);
}

void SpectrumDisplay::updateMagnitudeResponses()
{
    if (coeffs.bassCoeffsGui != nullptr)
    {
        coeffs.bassCoeffsGui->getMagnitudeForFrequencyArray (frequencies.data(),
            filterMagnitudeResponses[static_cast<size_t> (filters::bass)].data(),
            frequencies.size(),
            audioProcessor.getSampleRate());
    }
    if (coeffs.trebleCoeffsGui != nullptr)
    {
        coeffs.trebleCoeffsGui->getMagnitudeForFrequencyArray (frequencies.data(),
            filterMagnitudeResponses[static_cast<size_t> (filters::treble)].data(),
            frequencies.size(),
            audioProcessor.getUpSampleRate());
    }

    // Fill highPass magnitude response with ones
    filterMagnitudeResponses[static_cast<size_t> (filters::highPass)].fill (1.0);

    switch (audioProcessor.getCurrentSettings().highPassSlope.load())
    { // Based on the slope, enable certain sections. Notice that the cases have no break, therefore the steeper slopes will enable all butterworth filters from the cases with gentler slopes as well.
        case Slope_24:
        {
            // When you get the high pass coeffients, that method knows the slope, so it will make one or two meaningful butterworth filters based on the current slope. The first filter (<1>) only affects slopes with 18 or 24 dB rolloff. The second filter (<0>) either works with the first filter to achieve 24/18 dB rolloff, or is used in isolation to affect slopes with 12 or 6 dB.
        }
        case Slope_18:
        {
            // Multiply highPass frequency response with linear gain values for this filter
            std::array<double, NUM_FREQUENCIES> magnitudeResponse {};
            if (coeffs.hpCoeffsGui[1] != nullptr)
            {
                coeffs.hpCoeffsGui[1]->getMagnitudeForFrequencyArray (frequencies.data(),
                    magnitudeResponse.data(),
                    frequencies.size(),
                    audioProcessor.getSampleRate());
                juce::FloatVectorOperations::multiply (filterMagnitudeResponses[static_cast<size_t> (filters::highPass)].data(), magnitudeResponse.data(), frequencies.size());
            }
        }
        case Slope_12:
        {
        }
        case Slope_6: // This will not be Butterworth
        {
            // Multiply highPass frequency response with linear gain values for this filter
            std::array<double, NUM_FREQUENCIES> magnitudeResponse {};
            if (coeffs.hpCoeffsGui[0] != nullptr)
            {
                coeffs.hpCoeffsGui[0]->getMagnitudeForFrequencyArray (frequencies.data(),
                    magnitudeResponse.data(),
                    frequencies.size(),
                    audioProcessor.getSampleRate());
                juce::FloatVectorOperations::multiply (filterMagnitudeResponses[static_cast<size_t> (filters::highPass)].data(), magnitudeResponse.data(), frequencies.size());
            }
        }
        case Off:
        {
        }
    }

    // Fill lowPass magnitude response with ones
    filterMagnitudeResponses[static_cast<size_t> (filters::lowPass)].fill (1.0);

    switch (audioProcessor.getCurrentSettings().lowPassSlope.load())
    { // Based on the slope, enable certain sections. Notice that the cases have no break, therefore the steeper slopes will enable all butterworth filters from the cases with gentler slopes as well.
        case Slope_24:
        { // When you get the low pass coeffients, that method knows the slope, so it will make one or two meaningful butterworth filters based on the current slope. The first filter (<1>) only affects slopes with 18 or 24 dB rolloff. The second filter (<0>) either works with the first filter to achieve 24/18 dB rolloff, or is used in isolation to affect slopes with 12 or 6 dB.
        }
        case Slope_18:
        {
            // Multiply lowPass frequency response with linear gain values for this filter
            std::array<double, NUM_FREQUENCIES> magnitudeResponse {};
            if (coeffs.lpCoeffsGui[1] != nullptr)
            {
                coeffs.lpCoeffsGui[1]->getMagnitudeForFrequencyArray (frequencies.data(),
                    magnitudeResponse.data(),
                    frequencies.size(),
                    audioProcessor.getSampleRate());
                juce::FloatVectorOperations::multiply (filterMagnitudeResponses[static_cast<size_t> (filters::lowPass)].data(), magnitudeResponse.data(), frequencies.size());
            }
        }
        case Slope_12:
        {
        }
        case Slope_6:
        {
            // Multiply lowPass frequency response with linear gain values for this filter
            std::array<double, NUM_FREQUENCIES> magnitudeResponse {};
            if (coeffs.lpCoeffsGui[0] != nullptr)
            {
                coeffs.lpCoeffsGui[0]->getMagnitudeForFrequencyArray (frequencies.data(),
                    magnitudeResponse.data(),
                    frequencies.size(),
                    audioProcessor.getSampleRate());
                juce::FloatVectorOperations::multiply (filterMagnitudeResponses[static_cast<size_t> (filters::lowPass)].data(), magnitudeResponse.data(), frequencies.size());
            }
        }
        case Off:
        {
        }
    }

    // Fill special curve magnitude response with ones
    filterMagnitudeResponses[static_cast<size_t> (filters::special)].fill (1.0);

    for (auto&& specialCurveCoeff : coeffs.specialCoeffsGui)
    {
        std::array<double, NUM_FREQUENCIES> magnitudeResponse {};
        specialCurveCoeff->getMagnitudeForFrequencyArray (frequencies.data(),
            magnitudeResponse.data(),
            frequencies.size(),
            audioProcessor.getSampleRate());
        juce::FloatVectorOperations::multiply (filterMagnitudeResponses[static_cast<size_t> (filters::special)].data(), magnitudeResponse.data(), frequencies.size());
    }
}

void SpectrumDisplay::resized()
{
}

// ============================== from frequalizer ================================

void SpectrumDisplay::updateMagnitudeResponsePath (juce::Rectangle<int> plotFrame)
{
    magnitudeResponse.clear();
    createFrequencyPlot (magnitudeResponse, getOverallMagnitudeResponse(), plotFrame, 30.f);
}

//===========================================// Response Curve

ResponseCurveComponent::ResponseCurveComponent (PlayBackEQAudioProcessor& p, Gui::TopPanel& pp) : audioProcessor (p), presetPanel (pp), spectrumDisplay (audioProcessor)
{

    // Set this component up as a listener for every parameter so that we can redraw the response curve whenver a parameter is changed
    const auto& params = audioProcessor.getParameters();
    for (const auto param : params)
    {
        param->addListener (this);
    }

    startTimerHz (refreshRate);

    addAndMakeVisible (spectrumDisplay);
    addAndMakeVisible (coordinateComponent);
    coordinateComponent.addMouseListener (this, true);
    coordinateComponent.setMouseCursor (juce::MouseCursor::CrosshairCursor);
}

ResponseCurveComponent::~ResponseCurveComponent()
{
    const auto& params = audioProcessor.getParameters();
    for (const auto param : params)
    {
        param->removeListener (this);
    }
}

void ResponseCurveComponent::parameterValueChanged (int parameterIndex, float newValue)
{
    parametersChanged.set (true);
}

void ResponseCurveComponent::timerCallback()
{
    if (checkRefreshB4Repaint)
    {
        if (presetPanel.getRefreshResponse() == false)
            return;
    }

    // const auto fftBounds = getAnalysisArea().toFloat();
    // spectrumDisplay.pathProducerProcess (fftBounds);

    spectrumDisplay.updateCoeffsGUI();
    spectrumDisplay.updateMagnitudeResponses();
    spectrumDisplay.updateOverallMagResponse();
    spectrumDisplay.updateMagnitudeResponsePath (spectrumDisplay.getLocalBounds());
    spectrumDisplay.repaint();
}

void ResponseCurveComponent::paint (juce::Graphics& g)
{
    // TRACE_COMPONENT();

    // g.fillAll (juce::Colours::black);
    drawBackgroundGrid (g);
    g.setColour (colors.textColor);
    drawTextLabels (g);
}

void ResponseCurveComponent::paintOverChildren (juce::Graphics& g)
{
    // TRACE_COMPONENT();

    if (isMouseOverAnalysisArea)
    {
        juce::String coord = "(";
        if (freq < 100)
        {
            coord += juce::String (freq, 1);
            coord += "Hz, ";
        }
        else if (freq < 1000)
        {
            coord += juce::String (round (freq));
            coord += "Hz, ";
        }
        else if (freq < 10000)
        {
            coord += juce::String (freq / 1000, 2);
            coord += "kHz, ";
        }
        else
        {
            coord += juce::String (freq / 1000, 1);
            coord += "kHz, ";
        }

        if (mag > 0)
        {
            coord += "+";
        }
        coord += juce::String (mag, 1);
        coord += "dB)";

        constexpr int fontHeight = 12;
        g.setFont (fontHeight);

        int xLoc = mouseX + 7;
        if (xLoc + juce::GlyphArrangement::getStringWidthInt (g.getCurrentFont(), coord) > getAnalysisArea().getRight())
        {
            xLoc = mouseX - juce::GlyphArrangement::getStringWidthInt (g.getCurrentFont(), coord) - 6;
        }

        const int yLoc = mouseY - 18;

        const auto coordRect = juce::Rectangle<int> (xLoc, yLoc, juce::GlyphArrangement::getStringWidthInt (g.getCurrentFont(), coord), fontHeight);
        g.setColour (juce::Colours::black.withAlpha (0.5f));
        g.fillRect (coordRect);
        g.setColour (colors.textColor);
        g.drawFittedText (coord, coordRect, juce::Justification::centred, 1);
    }
}

void ResponseCurveComponent::mouseEnter (const juce::MouseEvent& event)
{
    if (event.eventComponent == &coordinateComponent)
    {
        isMouseOverAnalysisArea = true;
    }
}

void ResponseCurveComponent::mouseExit (const juce::MouseEvent& event)
{
    if (event.eventComponent == &coordinateComponent)
    {
        isMouseOverAnalysisArea = false;
    }
}

void ResponseCurveComponent::mouseMove (const juce::MouseEvent& event)
{
    if (event.eventComponent == &coordinateComponent)
    {
        const auto responseArea = coordinateComponent.getBounds();
        const auto w = responseArea.getWidth();
        freq = juce::mapToLog10 (event.position.x / static_cast<float> (w), 20.f, 20000.f);
        const auto outputMin = static_cast<float> (responseArea.getHeight());
        mag = juce::jmap (event.position.y, outputMin, 0.f, -30.f, 30.f);
        const auto relEvent = event.getEventRelativeTo (this);
        mouseX = relEvent.x;
        mouseY = relEvent.y;
    }
}

juce::Rectangle<int> ResponseCurveComponent::getAnalysisArea() const
{
    auto bounds = getLocalBounds();

    bounds.removeFromTop (22);
    bounds.removeFromBottom (4);
    bounds.removeFromLeft (32);
    bounds.removeFromRight (7);

    return bounds;
}

std::vector<float> ResponseCurveComponent::getFrequencies()
{
    return std::vector<float> {
        20, 100, 1000, 10000, 20000
    };
}

std::vector<float> ResponseCurveComponent::getGains()
{
    return std::vector<float> {
        -30, -20, -10, 0, 10, 20, 30
    };
}

std::vector<float> ResponseCurveComponent::getXs (const std::vector<float>& freqs, const float left, const float width)
{
    // Get horizontal positons for each frequency tic

    std::vector<float> xs;
    for (const auto f : freqs)
    {
        const auto normX = juce::mapFromLog10 (f, 20.f, 20000.f);
        xs.push_back (left + width * normX);
    }

    return xs;
}

void ResponseCurveComponent::drawBackgroundGrid (juce::Graphics& g) const
{
    using namespace juce;
    const auto freqs = getFrequencies();

    const auto renderArea = getAnalysisArea();
    const auto left = renderArea.getX();
    const auto right = renderArea.getRight();
    const auto top = renderArea.getY();
    const auto bottom = renderArea.getBottom();
    const auto width = renderArea.getWidth();

    const auto xs = getXs (freqs, left, width);

    // Draw vertical lines at each frequency tic
    g.setColour (Colour (60u, 60u, 60u));
    for (const auto x : xs)
    {
        g.drawVerticalLine (x, top, bottom);
    }

    // Draw horizontal lines at each level tic
    const auto gain = getGains();
    for (const auto gDb : gain)
    {
        const auto y = jmap (gDb, -30.f, 30.f, static_cast<float> (bottom), static_cast<float> (top));
        g.drawHorizontalLine (y, left, right);
    }
}

void ResponseCurveComponent::drawTextLabels (juce::Graphics& g) const
{
    using namespace juce;
    g.setColour (Colours::lightgrey);
    constexpr int fontHeight = 12;
    g.setFont (fontHeight);

    const auto renderArea = getAnalysisArea();
    // const auto left = renderArea.getX();

    const auto top = renderArea.getY();
    const auto bottom = renderArea.getBottom();
    const auto gain = getGains();
    for (const auto gDb : gain)
    {
        const auto y = jmap (gDb, -30.f, 30.f, static_cast<float> (bottom), static_cast<float> (top));

        String str;
        if (gDb > 0)
            str << "+";
        str << gDb;

        if (gDb == 0)
            str << "dB";

        const auto textWidth = GlyphArrangement::getStringWidthInt (g.getCurrentFont(), str);

        Rectangle<int> r;
        r.setSize (textWidth, fontHeight);
        r.setCentre (16, y);

        g.setColour (colors.textColor);
        g.drawFittedText (str, r, juce::Justification::centred, 1);
    }
}

void ResponseCurveComponent::resized()
{
    auto coordinateBounds = getAnalysisArea();
    coordinateBounds.translate (2, 0);
    coordinateComponent.setBounds (coordinateBounds);
    spectrumDisplay.setBounds (getAnalysisArea());

    //    responseCurve.preallocateSpace(getWidth() * 3);
}

//=================================================================
// void PathProducer::process (const juce::Rectangle<float> fftBounds, const double sampleRate)
// {
//     /** This function will:
//      - Query the processor FIFO to see if any buffers are available
//         - If so shift the oldest data out of the beginning of monoBuffer to make room for the new buffer
//         - Send the monoBuffer to the FFTDataGenerator
//         - (If the size of the FIFO buffers are smaller than the aumount of data required by the FFT, then this has the effect of a moving average)
//      - Query the FFTDataGenerator to see if there are processed FFT data blocks availabe
//         - If so send to pathGenerator
//      - Query the pathGenerator to see if there are paths available
//         - If so, set FFTpath to the new path - FFTpath will SpectrumDisplay::paint to draw the spectrum display
//     */
//
//     // Calculate the appropriate FFT order to pass to FFTDataGenerator.
//     // FFTDataGenerator will change its FFT order if they are different
//     const int newFFTOrder = ceil (log2 (sampleRate / 6));
//
//     juce::AudioBuffer<float> tempIncomingBuffer;
//
//     while (Fifo->getNumCompleteBuffersAvailable() > 0)
//     {
//         if (Fifo->getAudioBuffer (tempIncomingBuffer))
//         {
//             const auto size = tempIncomingBuffer.getNumSamples();
//
//             // Shift every value in monoBuffer starting at index size to the left by size positions
//             juce::FloatVectorOperations::copy (monoBuffer.getWritePointer (0, 0), // Destination
//                 monoBuffer.getReadPointer (0, size), // Source
//                 monoBuffer.getNumSamples() - size); // Number of values
//
//             // Copy the new buffer from Fifo into the end of monoBuffer
//             juce::FloatVectorOperations::copy (monoBuffer.getWritePointer (0, monoBuffer.getNumSamples() - size), // Destination
//                 tempIncomingBuffer.getReadPointer (0, 0), // Source
//                 size); // Number of values
//
//             fftDataGenerator.produceFFTDataForRendering (monoBuffer, newFFTOrder);
//         }
//     }
//
//     const auto fftSize = fftDataGenerator.getFFTSize();
//     const auto binWidth = sampleRate / static_cast<double> (fftSize);
//
//     while (fftDataGenerator.getNumAvailableFFTDataBlocks() > 0)
//     {
//         if (fftDataGenerator.getFFTData (fftData))
//         {
//             pathGenerator.generatePath (fftData, fftBounds, fftSize, binWidth);
//         }
//     }
//
//     while (pathGenerator.getNumPathsAvailable() > 0)
//     {
//         pathGenerator.getPath (FFTPath);
//     }
// }
