
#include "SimpleEQ_PluginEditor.h"

void LookAndFeel::drawRotarySlider (juce::Graphics& g, int x, int y, int width, int height, float sliderPos, const float rotaryStartAngle, const float rotaryEndAngle, juce::Slider& slider)
{
    using namespace juce;

    auto outline = slider.findColour (Slider::rotarySliderOutlineColourId);
    auto fill = slider.findColour (Slider::rotarySliderFillColourId);

    auto bounds = Rectangle<int> (x, y, width, height).toFloat().reduced (16);

    auto radius = jmin (bounds.getWidth(), bounds.getHeight()) / 2.0f;
    auto toAngle = rotaryStartAngle + sliderPos * (rotaryEndAngle - rotaryStartAngle); // Angle corrosponding to the current value of the control
    auto lineW = 1.6f * (jmin (8.0f, radius * 0.5f));
    auto arcRadius = radius - lineW * 0.5f;

    // This arc spans the angles corrosponding the full range of possible values
    Path backgroundArc;
    backgroundArc.addCentredArc (bounds.getCentreX(),
        bounds.getCentreY(),
        arcRadius,
        arcRadius,
        0.0f,
        rotaryStartAngle,
        rotaryEndAngle,
        true);

    g.setColour (outline);
    g.strokePath (backgroundArc, PathStrokeType (lineW, PathStrokeType::curved, PathStrokeType::butt));

    // This arc is the portion of the arc that is below the current value angle
    Path valueArc;
    valueArc.addCentredArc (bounds.getCentreX(),
        bounds.getCentreY(),
        arcRadius,
        arcRadius,
        0.0f,
        rotaryStartAngle,
        toAngle,
        true);
    g.setColour (fill);
    g.strokePath (valueArc, PathStrokeType (lineW, PathStrokeType::curved, PathStrokeType::butt));

    // The thumb is a section of the background arc
    Path thumbArc;
    thumbArc.addCentredArc (bounds.getCentreX(),
        bounds.getCentreY(),
        arcRadius,
        arcRadius,
        0.0f,
        toAngle - 0.18f,
        toAngle + 0.18f,
        true);

    // If control is bypassed, make the thumb color dim
    if (slider.isEnabled())
    {
        g.setColour (slider.findColour (Slider::thumbColourId));
    }
    else
    {
        g.setColour (slider.findColour (Slider::thumbColourId).withBrightness (0.4));
    }

    // Draw the thumb
    g.strokePath (thumbArc, PathStrokeType (lineW, PathStrokeType::curved, PathStrokeType::butt));
}

void RotarySliderWithLabels::paint (juce::Graphics& g)
{
    using namespace juce;

    auto startAng = degreesToRadians (180.f + 45.f);
    auto endAng = degreesToRadians (180.f - 45.f) + MathConstants<float>::twoPi;

    auto range = getRange();

    auto sliderBounds = getSliderBounds();

    getLookAndFeel().drawRotarySlider (g,
        sliderBounds.getX(),
        sliderBounds.getY(),
        sliderBounds.getWidth(),
        sliderBounds.getHeight(),
        jmap (getValue(), range.getStart(), range.getEnd(), 0.0, 1.0),
        startAng,
        endAng,
        *this);

    auto center = sliderBounds.toFloat().getCentre();
    auto radius = sliderBounds.getWidth() * 0.4f;

    // The text will be dark grey if the control is bypassed
    if (isEnabled())
    {
        g.setColour (colors.textColor);
    }
    else
    {
        g.setColour (juce::Colours::darkgrey);
    }

    g.setFont (getTextHeight() * 1.f);

    // Tic labels - min and max values for a continuous control, each value for a rotary switch
    auto numChoices = labels.size();
    for (int i = 0; i < numChoices; ++i)
    {
        auto pos = labels[i].pos;
        jassert (0.f <= pos);

        auto ang = jmap (pos, 0.f, static_cast<float> (numChoices) - 1.f, startAng, endAng);

        auto c = center.getPointOnCircumference (radius + getTextHeight() * 0.5f + 1.f, ang);

        // Some height corrections for certain controls - this might be unnecessary with some adjustments to how c is calculated and/or how the text is justified?
        if (numChoices == 2)
        {
            c.y = c.y - 8;
        }

        if (title == "HP\nSlope\n(dB/oct.)" || title == "LP\nSlope\n(dB/oct.)")
        {
            if (i == 0 or i == 4)
            {
                c.y = c.y - 8;
            }
            if (i == 1)
            {
                c.y = c.y - 15;
            }
            if (i == 2)
            {
                c.y = c.y - 12;
            }
            if (i == 3)
            {
                c.y = c.y - 15;
            }
        }

        // Create rectangle for text bounds
        Rectangle<float> r;
        auto str = labels[i].label;
        r.setSize (g.getCurrentFont().getStringWidth (str), getTextHeight());
        r.setCentre (c);
        r.setY (r.getY() + getTextHeight());

        // Draw tic labels
        g.drawFittedText (str, r.toNearestInt(), juce::Justification::centred, 1);
    }

    // Name of control
    Rectangle<float> r;
    r.setSize (g.getCurrentFont().getStringWidth (title), getTextHeight());
    r.setCentre (center);
    g.drawFittedText (title, r.toNearestInt(), juce::Justification::centred, 1);
}

juce::Rectangle<int> RotarySliderWithLabels::getSliderBounds() const
{
    return getLocalBounds();
}

//===========================================// Spectrum Display

SpectrumDisplay::SpectrumDisplay (PlayBackEQAudioProcessor& p) : audioProcessor (p),
                                                                 leftPathProducer (audioProcessor.spectrumAnalyzerSampleFifo, p, -80)
{
    setInterceptsMouseClicks (false, false); // Passes mouse events through to component behind it
}

SpectrumDisplay::~SpectrumDisplay() = default;

void SpectrumDisplay::pathProducerProcess (juce::Rectangle<float> fftBounds)
{
    leftPathProducer.process (fftBounds, audioProcessor.getSampleRate());
}

void SpectrumDisplay::paint (juce::Graphics& g)
{
    using namespace juce;

    // Get FFT path, fill  it and outline it
    auto FFTpath = leftPathProducer.getPath();
    g.setColour (displayColor);
    //    float alphaVal = 0.7;
    //    g.setColour (juce::Colours::whitesmoke.withMultipliedAlpha(alphaVal));

    g.setOpacity (0.6);
    g.fillPath (FFTpath);
    g.strokePath (FFTpath, PathStrokeType (1.f));

    // For overlaying multiple paths with some weighted transparency - might be useful at some point
    //    float opaque = 0.5f;
    //    for (int pathNum = 0; pathNum < numPaths; pathNum++) {
    //        int idx = (pathCount+pathNum)%numPaths;
    //        g.setOpacity(opaque);
    //        g.fillPath(FFTpaths[idx]);
    //        g.strokePath(FFTpaths[idx], PathStrokeType(1.f));
    //        opaque = opaque/2.f;
    //    }
    //    pathCount++;
    //    if (pathCount%numPaths == 0){
    //        pathCount = 0;
    //    }

    // Draw response curve
    g.setColour (Colours::white.withBrightness (0.9));
    g.strokePath (responseCurve, PathStrokeType (2.f));
}

void SpectrumDisplay::updateFiltersGUI()
{
    // Updates the filters owned by the editor with the current control values

    auto settings = audioProcessor.getSettings (audioProcessor.apvts);
    double sampleRate = audioProcessor.getSampleRate();

    // Calculate bass filter coefficients
    auto bassCoeffs = getBassCoeffs (settings, sampleRate);
    // Convert array to coefficients pointer
    juce::dsp::IIR::Coefficients<float>::Ptr bassCoeffsPtr (new juce::dsp::IIR::Coefficients<float> (Coefficients (bassCoeffs)));
    // Apply coefficients to filters
    updateCoefficients (leftBassChain.get<0>().coefficients, bassCoeffsPtr);
    // Apply bass boost bypass
    leftBassChain.setBypassed<0> (!settings.bassBoostEnable);

    // Calculate treble filter coefficients
    auto trebleCoeffs = getTrebleAttenCoeffs (settings, audioProcessor.getUpSampleRate());
    // Convert array to coefficients pointer
    juce::dsp::IIR::Coefficients<float>::Ptr trebleCoeffsPtr (new juce::dsp::IIR::Coefficients<float> (Coefficients (trebleCoeffs)));
    // Apply coefficients to filters
    updateCoefficients (leftTrebleChain.get<0>().coefficients, trebleCoeffsPtr);
    // Apply treble bypass
    leftTrebleChain.setBypassed<0> (!settings.trebleCutEnable);

    // Get inverse RIAA filter coefficients
    auto RIAAcoeffs = getRIAACoeffs (settings, audioProcessor.getUpSampleRate());
    // Convert array to coefficients pointer
    juce::dsp::IIR::Coefficients<float>::Ptr RIAACoeffsPtr (new juce::dsp::IIR::Coefficients<float> (Coefficients (RIAAcoeffs)));
    // Apply coefficients to filters
    updateCoefficients (leftRIAAChain.get<0>().coefficients, RIAACoeffsPtr);
    // Apply RIAA bypass
    leftRIAAChain.setBypassed<0> (!settings.RIAAEnable);

    // Get high pass coefficients
    auto hpfcoeffs = getHighPasscoeffs (settings, sampleRate);
    // Convert array to coefficients pointer
    juce::dsp::IIR::Coefficients<float>::Ptr hpfcoeffsPtr (new juce::dsp::IIR::Coefficients<float> (Coefficients (hpfcoeffs)));
    // apply coefficients to filters depending on slope
    leftHPChain.setBypassed<0> (true);
    leftHPChain.setBypassed<1> (true);
    leftHPChain.setBypassed<2> (true);

    // Get coefficients for butterworth filters (any order over 1st)
    auto hpfcoeffsPtrButter = getHighPassCoeffsButter (settings, sampleRate);

    switch (settings.highPassSlope)
    { // Based on the slope, enable certain sections. Notice that the cases have no break, therefore the steeper slopes will enable all of the SOS from the cases with gentler slopes as well.
        case Slope_24:
        {
        }
        case Slope_18:
        {
            updateCoefficients (leftHPChain.get<2>().coefficients, hpfcoeffsPtrButter[1]);
            leftHPChain.setBypassed<2> (false);
        }
        case Slope_12:
        {
            updateCoefficients (leftHPChain.get<1>().coefficients, hpfcoeffsPtrButter[0]);
            leftHPChain.setBypassed<1> (false);
        }
        break;
        case Slope_6: // This will not be Butterworth
        {
            updateCoefficients (leftHPChain.get<0>().coefficients, hpfcoeffsPtr);
            leftHPChain.setBypassed<0> (false);
        }
        case Off:
        {
        }
    }

    // get low pass coefficients
    auto lpfcoeffsPtr = getLowPasscoeffs (settings, audioProcessor.getSampleRate());
    // apply coefficients to filters depending on slope
    leftLPChain.setBypassed<0> (true);
    leftLPChain.setBypassed<1> (true);
    switch (settings.lowPassSlope)
    {
        case Slope_24:
        {
        }
        case Slope_18:
        {
            updateCoefficients (leftLPChain.get<1>().coefficients, lpfcoeffsPtr[1]);
            leftLPChain.setBypassed<1> (false);
        }
        case Slope_12:
        {
        }
        case Slope_6:
        {
            updateCoefficients (leftLPChain.get<0>().coefficients, lpfcoeffsPtr[0]);
            leftLPChain.setBypassed<0> (false);
        }
        case Off:
        {
        }
    }
}

void SpectrumDisplay::updateResponseCurve()
{
    using namespace juce;

    auto responseArea = getLocalBounds();
    auto w = responseArea.getWidth();

    auto userPreferences = audioProcessor.getUserPreferences();
    auto sampleRate = audioProcessor.getSampleRate();
    auto upSampleRate = audioProcessor.getUpSampleRate();

    // This will be a vector of magnitude values for each pixel/frequency that we will conbvert to dB and then turn into a path
    std::vector<double> mags;

    mags.resize (w);

    responseCurve.clear();

    // Create map to convert values in dB to vertical position in responseArea
    const double outputMin = responseArea.getBottom();
    const double outputMax = responseArea.getY();
    auto map = [outputMin, outputMax] (double input) {
        return jmap (input, -30.0, 30.0, outputMin, outputMax);
    };

    for (int pix = 0; pix < w; ++pix)
    {
        // Magnitude for every frequency is 1 if no filters are modifying it
        double mag = 1.f;

        // Find the frequency asscoiated with this pixel in horizontal dimension of responseArea
        auto freq = mapToLog10 (static_cast<double> (pix) / static_cast<double> (w), 20.0, 20000.0);

        // Get the magnitude repsonse at this frequendcy from each filter and multiply them together
        if (!leftTrebleChain.isBypassed<0>())
        {
            mag *= leftTrebleChain.get<0>().coefficients->getMagnitudeForFrequency (freq, upSampleRate);
        }

        if (!leftBassChain.isBypassed<0>())
        {
            mag *= leftBassChain.get<0>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }

        if (!leftHPChain.isBypassed<0>())
        {
            mag *= leftHPChain.get<0>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }
        if (!leftHPChain.isBypassed<1>())
        {
            mag *= leftHPChain.get<1>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }
        if (!leftHPChain.isBypassed<2>())
        {
            mag *= leftHPChain.get<2>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }

        if (!leftLPChain.isBypassed<0>())
        {
            mag *= leftLPChain.get<0>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }
        if (!leftLPChain.isBypassed<1>())
        {
            mag *= leftLPChain.get<1>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }

        if (!leftRIAAChain.isBypassed<0>() && userPreferences.showRIAA)
        {
            mag *= leftRIAAChain.get<0>().coefficients->getMagnitudeForFrequency (freq, upSampleRate);
        }

        // Magnitude to dB for this pixel/frequency
        mags[pix] = Decibels::gainToDecibels (mag);

        auto y = map (mags[pix]);
        if (pix == 0)
        {
            // Start path
            responseCurve.startNewSubPath (responseArea.getX(), y);
        }
        else
        {
            // Build path
            responseCurve.lineTo (responseArea.getX() + pix, y);
        }
    }
}

void SpectrumDisplay::resized()
{
    updateResponseCurve();
}

//===========================================// Response Curve

ResponseCurveComponent::ResponseCurveComponent (PlayBackEQAudioProcessor& p) : audioProcessor (p),
                                                                               spectrumDisplay (audioProcessor)
{
    // For overlaying multiple paths with some weighted transparency - might be useful at some point
    //    FFTpaths.clear();
    //    FFTpaths.resize(numPaths, path);

    // Set this component up as a listener for every parameter so that we can redraw the response curve whenver a parameter is changed
    const auto& params = audioProcessor.getParameters();
    for (auto param : params)
    {
        param->addListener (this);
    }

    startTimerHz (60);

    addAndMakeVisible (spectrumDisplay);
    addAndMakeVisible (coordinateComponent);
    coordinateComponent.addMouseListener (this, true);
}

ResponseCurveComponent::~ResponseCurveComponent()
{
    const auto& params = audioProcessor.getParameters();
    for (auto param : params)
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
    auto fftBounds = getAnalysisArea().toFloat();

    spectrumDisplay.setBounds (getAnalysisArea());

    spectrumDisplay.pathProducerProcess (fftBounds);

    if (parametersChanged.compareAndSetBool (false, true))
    {
        spectrumDisplay.updateFiltersGUI();
        spectrumDisplay.updateResponseCurve();
    }

    repaint();
}

void ResponseCurveComponent::paint (juce::Graphics& g)
{
    drawBackgroundGrid (g);
    g.setColour (colors.textColor);
    drawTextLabels (g);
}

void ResponseCurveComponent::paintOverChildren (juce::Graphics& g)
{
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

        g.setColour (colors.textColor);

        constexpr int fontHeight = 12;
        g.setFont (fontHeight);

        int xLoc = mouseX + 7;
        if (xLoc + g.getCurrentFont().getStringWidth (coord) > getAnalysisArea().getRight())
        {
            xLoc = mouseX - g.getCurrentFont().getStringWidth (coord) - 6;
        }

        int yLoc = mouseY - 12;

        auto coordRect = juce::Rectangle<int> (xLoc, yLoc, g.getCurrentFont().getStringWidth (coord), fontHeight);
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
        freq = juce::mapToLog10 (static_cast<float> (event.x) / static_cast<float> (w), 20.f, 20000.f);
        const auto outputMin = responseArea.getHeight();
        mag = juce::jmap (static_cast<float> (event.y), static_cast<float> (outputMin), 0.f, -30.f, 30.f);
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
    bounds.removeFromRight (20);

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

std::vector<float> ResponseCurveComponent::getXs (const std::vector<float>& freqs, float left, float width)
{
    // Get horizontal positons for each frequency tic

    std::vector<float> xs;
    for (auto f : freqs)
    {
        auto normX = juce::mapFromLog10 (f, 20.f, 20000.f);
        xs.push_back (left + width * normX);
    }

    return xs;
}

void ResponseCurveComponent::drawBackgroundGrid (juce::Graphics& g) const
{
    using namespace juce;
    auto freqs = getFrequencies();

    auto renderArea = getAnalysisArea();
    auto left = renderArea.getX();
    auto right = renderArea.getRight();
    auto top = renderArea.getY();
    auto bottom = renderArea.getBottom();
    auto width = renderArea.getWidth();

    auto xs = getXs (freqs, left, width);

    // Draw vertical lines at each frequency tic
    g.setColour (Colour (60u, 60u, 60u));
    for (auto x : xs)
    {
        g.drawVerticalLine (x, top, bottom);
    }

    // Draw horizontal lines at each level tic
    auto gain = getGains();
    for (auto gDb : gain)
    {
        auto y = jmap (gDb, -30.f, 30.f, static_cast<float> (bottom), static_cast<float> (top));
        g.drawHorizontalLine (y, left, right);
    }
}

void ResponseCurveComponent::drawTextLabels (juce::Graphics& g) const
{
    using namespace juce;
    g.setColour (Colours::lightgrey);
    constexpr int fontHeight = 12;
    g.setFont (fontHeight);

    auto renderArea = getAnalysisArea();
    auto left = renderArea.getX();

    auto top = renderArea.getY();
    auto bottom = renderArea.getBottom();
    auto width = renderArea.getWidth();

    // Draw frequency tic labels

    auto freqs = getFrequencies();
    auto xs = getXs (freqs, left, width);

    for (int i = 0; i < freqs.size(); ++i)
    {
        auto f = freqs[i];
        auto x = xs[i];

        bool addK = false;
        String str;
        if (f > 999.f)
        {
            addK = true;
            f /= 1000.f;
        }

        str << f;
        if (addK)
            str << "k";
        str << "Hz";

        auto textWidth = g.getCurrentFont().getStringWidth (str);

        Rectangle<int> r;

        r.setSize (textWidth, fontHeight);
        if (f == 20.f && !addK)
        {
            r.setCentre (x + textWidth / 2.f, 0);
        }
        else
        {
            r.setCentre (x, 0);
        }
        r.setY (fontHeight / 2);

        g.drawFittedText (str, r, juce::Justification::centred, 1);
    }

    // Draw dB tic labels

    auto gain = getGains();
    for (auto gDb : gain)
    {
        auto y = jmap (gDb, -30.f, 30.f, static_cast<float> (bottom), static_cast<float> (top));

        String str;
        if (gDb > 0)
            str << "+";
        str << gDb;

        if (gDb == 0)
            str << "dB";

        auto textWidth = g.getCurrentFont().getStringWidth (str);

        Rectangle<int> r;
        r.setSize (textWidth, fontHeight);
        r.setCentre (16, y);

        g.setColour (colors.textColor);
        g.drawFittedText (str, r, juce::Justification::centred, 1);

        // For tick labels on right for spectrum analyzer
        //        str.clear();
        //        str << (gDb - 30.f);
        //
        //        r.setX(1);
        //        textWidth = g.getCurrentFont().getStringWidth(str);
        //        r.setSize(textWidth, fontHeight);
        //        g.setColour(Colours::lightgrey);
        //        g.drawFittedText(str, r, juce::Justification::centred, 1);
    }
}

void ResponseCurveComponent::resized()
{
    auto coordinateBounds = getAnalysisArea();
    coordinateBounds.translate (1, 3);
    coordinateComponent.setBounds (coordinateBounds);
    //    responseCurve.preallocateSpace(getWidth() * 3);
}

//=================================================================
void PathProducer::process (juce::Rectangle<float> fftBounds, double sampleRate)
{
    /** This function will:
     - Query the processor FIFO to see if any buffers are available
        - If so shift the oldest data out of the beginning of monoBuffer to make room for the new buffer
        - Send the monoBuffer to the FFTDataGenerator
        - (If the size of the FIFO buffers are smaller than the aumount of data required by the FFT, then this has the effect of a moving average)
     - Query the FFTDataGenerator to see if there are processed FFT data blocks availabe
        - If so send to pathGenerator
     - Query the pathGenerator to see if there are paths available
        - If so, set FFTpath to the new path - FFTpath will SpectrumDisplay::paint to draw the spectrum display
    */

    // Calculate the appropriate FFT order to pass to FFTDataGenerator.
    // FFTDataGenerator will change its FFT order if they are different
    const int newFFTOrder = ceil (log2 (sampleRate / 6));

    juce::AudioBuffer<float> tempIncomingBuffer;

    while (Fifo->getNumCompleteBuffersAvailable() > 0)
    {
        if (Fifo->getAudioBuffer (tempIncomingBuffer))
        {
            auto size = tempIncomingBuffer.getNumSamples();

            // Shift every value in monoBuffer starting at index size to the left by size positions
            juce::FloatVectorOperations::copy (monoBuffer.getWritePointer (0, 0), // Destination
                monoBuffer.getReadPointer (0, size), // Source
                monoBuffer.getNumSamples() - size); // Number of values

            // Copy the new buffer from Fifo into the end of monoBuffer
            juce::FloatVectorOperations::copy (monoBuffer.getWritePointer (0, monoBuffer.getNumSamples() - size), // Destination
                tempIncomingBuffer.getReadPointer (0, 0), // Source
                size); // Number of values

            FFTDataGenerator.produceFFTDataForRendering (monoBuffer, newFFTOrder);
        }
    }

    const auto fftSize = FFTDataGenerator.getFFTSize();
    const auto binWidth = sampleRate / static_cast<double> (fftSize);

    while (FFTDataGenerator.getNumAvailableFFTDataBlocks() > 0)
    {
        if (FFTDataGenerator.getFFTData (fftData))
        {
            pathGenerator.generatePath (fftData, fftBounds, fftSize, binWidth);
        }
    }

    while (pathGenerator.getNumPathsAvailable() > 0)
    {
        pathGenerator.getPath (FFTPath);
    }
}
