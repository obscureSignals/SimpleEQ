
#include "SimpleEQ_PluginEditor.h"

void LookAndFeelShort::drawRotarySlider (juce::Graphics& g, const int x, const int y, const int width, const int height, const float sliderPos, const float rotaryStartAngle, const float rotaryEndAngle, juce::Slider& slider)
{
    using namespace juce;

    const auto outline = slider.findColour (Slider::rotarySliderOutlineColourId);
    const auto fill = slider.findColour (Slider::rotarySliderFillColourId);

    auto bounds = Rectangle<int> (x, y, width, height).toFloat().reduced (16);
    const auto radius = jmin (bounds.getWidth(), bounds.getHeight()) / 2.f;
    const auto toAngle = rotaryStartAngle + sliderPos * (rotaryEndAngle - rotaryStartAngle); // Angle corrosponding to the current value of the control
    const auto lineW = 1.6f * (jmin (8.0f, radius * 0.5f));
    const auto arcRadius = radius - lineW * 0.5f;

    const auto yCenter = bounds.getCentreY() - 13;

    // This arc spans the angles corrosponding the full range of possible values
    Path backgroundArc;
    backgroundArc.addCentredArc (bounds.getCentreX(),
        yCenter,
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
        yCenter,
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
        yCenter,
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

//===================================================

void RotarySliderWithLabelsShort::paint (juce::Graphics& g)
{
    // TRACE_COMPONENT();

    using namespace juce;

    constexpr auto startAng = degreesToRadians (180.f + 45.f);
    constexpr auto endAng = degreesToRadians (180.f - 45.f) + MathConstants<float>::twoPi;

    const auto sliderBounds = getSliderBounds();

    getLookAndFeel().drawRotarySlider (g,
        sliderBounds.getX(),
        sliderBounds.getY(),
        sliderBounds.getWidth(),
        sliderBounds.getHeight(),
        getNormalisableRange().convertTo0to1 (getValue()),
        startAng,
        endAng,
        *this);

    auto center = sliderBounds.toFloat().getCentre();
    center.setY (center.getY() - 13);

    const auto radius = jmin (getSliderBounds().reduced (16).getWidth(), getSliderBounds().reduced (16).getHeight()) / 2.f;


    if (isEnabled())
    {
        g.setColour (colors.textColor);
    }
    else
    {
        g.setColour (Colours::darkgrey);
    }

    // g.setFont (static_cast<float> (getTextHeight()));

    g.setFont (13);
    const auto numChoices = labels.size();
    for (int i = 0; i < numChoices; ++i)
    {
        const auto pos = labels[i].pos;
        jassert (0.f <= pos);

        const auto ang = jmap (pos, 0.f, static_cast<float> (numChoices) - 1.f, startAng, endAng);

        auto c = center.getPointOnCircumference (radius + getTextHeight() * 0.5f + 2.f, ang);

        // Create rectangle for text bounds
        Rectangle<float> r;
        auto str = labels[i].label;
        r.setSize (GlyphArrangement::getStringWidth (g.getCurrentFont(), str), getTextHeight());
        c.setY (c.getY() + 4);
        r.setCentre (c);

        // Draw tic labels
        auto test = g.getCurrentFont();
        g.drawFittedText (str, r.toNearestInt(), juce::Justification::centred, 1);
    }

    // Name of control
    Rectangle<float> r;
    // g.setFont (14);
    r.setSize (getMaxStringWidth (g.getCurrentFont(), title), getTextHeight());
    r.setCentre (center);
    g.drawFittedText (title, r.toNearestInt(), juce::Justification::centred, 1);
}

juce::Rectangle<int> RotarySliderWithLabelsShort::getSliderBounds() const
{
    return getLocalBounds();
}

double RotarySliderWithLabelsShort::snapValue (double attemptedValue, DragMode dragMode)
{
    if (snap && snapValues != nullptr && !getCommandControlPressed())
    {
        double minVal = INFINITY;
        int idxMin = 0;
        for (auto idx = 0; idx < snapValues->size(); idx++)
        {
            const auto val = std::abs ((*snapValues)[idx] - attemptedValue);
            if (val < minVal)
            {
                minVal = val;
                idxMin = idx;
            }
        }

        if (isMouseDrag || isMouseDown) // and command/ctrl is NOT pressed and snap is set
        {
            return (*snapValues)[idxMin];
        }

        if (isMouseWheelMove) // and command/ctrl is NOT pressed and snap is set
        {
            // If this is coming within ~160ms of the last time we processed a mouse wheel move, ignore it. This keeps a small mouse wheel move from moving the slider all the way to the other side.
            auto elapsedTime = juce::Time().currentTimeMillis() - getLastMouseWheelMove();
            if (elapsedTime < 160)
            {
                return getValue();
            }
            setLastMouseWheelMove (juce::Time().currentTimeMillis());

            // A single mouse wheel move is usually not enough to move the slider value far enough so that it snaps to the next value. Instead, we look to see which way the slider is being moved and jump to the next snap value in that direction.
            if (juce::approximatelyEqual ((*snapValues)[idxMin], static_cast<float>(getValue())), 2)
            {
                if (attemptedValue < (*snapValues)[idxMin])
                {
                    return (*snapValues)[idxMin - 1];
                }
                if (attemptedValue > (*snapValues)[idxMin])
                {
                    return (*snapValues)[idxMin + 1];
                }
            }
            return (*snapValues)[idxMin];
        }
    }
    return attemptedValue; // command/ctrl pressed, value entered in text box, or no snap values assigned
}

void RotarySliderWithLabelsShort::mouseDrag (const juce::MouseEvent& event)
{
    setIsMouseDrag (true);
    if (event.mods == 24 or event.mods == 18)
    {
        setCommandControlPressed (true);
    }
    Slider::mouseDrag (event);
    setCommandControlPressed (false);
    setIsMouseDrag (false);
}

void RotarySliderWithLabelsShort::mouseWheelMove (const juce::MouseEvent& e, const juce::MouseWheelDetails& wheel)
{
    setIsMouseWheelMove (true);
    if (e.mods == 8 or e.mods == 2)
    {
        setCommandControlPressed (true);
    }
    Slider::mouseWheelMove (e, wheel);
    setCommandControlPressed (false);
    setIsMouseWheelMove (false);
}

void RotarySliderWithLabelsShort::mouseDown (const juce::MouseEvent& event)
{
    setIsMouseDown (true);
    if (event.mods == 24 or event.mods == 18)
    {
        setCommandControlPressed (true);
    }
    Slider::mouseDown (event);
    setCommandControlPressed (false);
    setIsMouseDown (false);
}

void RotarySliderWithLabelsShort::mouseEnter (const juce::MouseEvent& event)
{
    Slider::mouseEnter (event);

    if (legend != nullptr && isEnabled())
    {
        // setAppearing (false);
        // startTimerHz (20);
        legend->setVisible (true);
        legend->setAppearing (true);
        legend->startTimerHz (20);
    }
}
void RotarySliderWithLabelsShort::mouseExit (const juce::MouseEvent& event)
{
    Slider::mouseExit (event);

    if (legend != nullptr && isEnabled())
    {
        // setAppearing (true);
        // startTimerHz (20);
        legend->setAppearing (false);
        legend->startTimerHz (20);
    }
}

void RotarySliderWithLabelsShort::timerCallback()
{
    // if (appearing)
    // {
    //     currentAlpha = currentAlpha + 0.1f;
    //     if (currentAlpha >= 1.0f)
    //     {
    //         currentAlpha = 1.0f;
    //         stopTimer();
    //     }
    // }
    // else
    // {
    //     currentAlpha = currentAlpha - 0.1f;
    //     if (currentAlpha <= 0)
    //     {
    //         currentAlpha = 0.f;
    //         stopTimer();
    //         return;
    //     }
    // }
}

//=========================================================

void RotarySliderLegendComponent::paint (juce::Graphics& g)
{
    using namespace juce;

    constexpr auto startAng = degreesToRadians (180.f + 45.f);
    constexpr auto endAng = degreesToRadians (180.f - 45.f) + MathConstants<float>::twoPi;

    const auto bounds = getLocalBounds();

    auto center = bounds.toFloat().getCentre();
    center.setY (center.getY() - 13);

    const auto radius = jmin (bounds.reduced (16).getWidth(), bounds.reduced (16).getHeight()) / 2.f;

    // The text will be dark grey if the control is bypassed
    if (isEnabled())
    {
        g.setColour (colors.textColor.withAlpha (currentAlpha));
    }
    else
    {
        g.setColour (Colours::darkgrey.withAlpha (currentAlpha));
    }

    g.setFont (13.f);

    for (int i = 1; i < labels.size() - 1; ++i)
    {
        const auto pos = labels[i].pos;
        jassert (0.f <= pos);

        const auto ang = jmap (pos, 0.f, static_cast<float> (labels.size()) - 1.f, startAng, endAng);

        auto c = center.getPointOnCircumference (radius - 2, radius - 6, ang);

        // Create rectangle for text bounds
        Rectangle<float> r;
        auto str = labels[i].label;
        r.setSize (GlyphArrangement::getStringWidth (g.getCurrentFont(), str), 14);
        r.setCentre (c);

        // Draw tic labels
        g.drawFittedText (str, r.toNearestInt(), juce::Justification::centred, 1);
    }
}

//=========================================================

void RotarySliderLegendComponent::timerCallback()
{
    if (appearing)
    {
        currentAlpha = currentAlpha + 0.1f;
        if (currentAlpha >= 1.0f)
        {
            currentAlpha = 1.0f;
            stopTimer();
        }
    }
    else
    {
        currentAlpha = currentAlpha - 0.1f;
        if (currentAlpha <= 0)
        {
            currentAlpha = 0.f;
            stopTimer();
            setVisible (false);
            return;
        }
    }
    repaint();
}

void LookAndFeel::drawRotarySlider (juce::Graphics& g, const int x, const int y, const int width, const int height, const float sliderPos, const float rotaryStartAngle, const float rotaryEndAngle, juce::Slider& slider)
{
    using namespace juce;

    const auto outline = slider.findColour (Slider::rotarySliderOutlineColourId);
    const auto fill = slider.findColour (Slider::rotarySliderFillColourId);

    const auto bounds = Rectangle<int> (x, y, width, height).toFloat().reduced (16);

    const auto radius = jmin (bounds.getWidth(), bounds.getHeight()) / 2.0f;
    const auto toAngle = rotaryStartAngle + sliderPos * (rotaryEndAngle - rotaryStartAngle); // Angle corrosponding to the current value of the control
    const auto lineW = 1.6f * (jmin (8.0f, radius * 0.5f));
    const auto arcRadius = radius - lineW * 0.5f;

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
    // TRACE_COMPONENT();

    using namespace juce;

    constexpr auto startAng = degreesToRadians (180.f + 45.f);
    constexpr auto endAng = degreesToRadians (180.f - 45.f) + MathConstants<float>::twoPi;

    const auto range = getRange();

    const auto sliderBounds = getSliderBounds();

    getLookAndFeel().drawRotarySlider (g,
        sliderBounds.getX(),
        sliderBounds.getY(),
        sliderBounds.getWidth(),
        sliderBounds.getHeight(),
        getNormalisableRange().convertTo0to1 (getValue()),
        startAng,
        endAng,
        *this);

    const auto center = sliderBounds.toFloat().getCentre();
    const auto radius = sliderBounds.getWidth() * 0.4f;

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
    const auto numChoices = labels.size();
    for (int i = 0; i < numChoices; ++i)
    {
        const auto pos = labels[i].pos;
        jassert (0.f <= pos);

        const auto ang = jmap (pos, 0.f, static_cast<float> (numChoices) - 1.f, startAng, endAng);

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
        r.setSize (GlyphArrangement::getStringWidth (g.getCurrentFont(), str), getTextHeight());
        r.setCentre (c);
        r.setY (r.getY() + getTextHeight());

        // Draw tic labels
        g.drawFittedText (str, r.toNearestInt(), juce::Justification::centred, 1);
    }

    // Name of control
    Rectangle<float> r;
    r.setSize (getMaxStringWidth (g.getCurrentFont(), title), getTextHeight());
    r.setCentre (center);
    g.drawFittedText (title, r.toNearestInt(), juce::Justification::centred, 1);
}

juce::Rectangle<int> RotarySliderWithLabels::getSliderBounds() const
{
    return getLocalBounds();
}

void RotarySliderWithLabels::mouseWheelMove (const juce::MouseEvent& e, const juce::MouseWheelDetails& wheel)
{
    // If this is coming within ~160ms of the last time we processed a mouse wheel move, ignore it. This keeps a small mouse wheel move from moving the slider all the way to the other side for controls with stops
    auto elapsedTime = juce::Time().currentTimeMillis() - getLastMouseWheelMove();
    if (elapsedTime < 160)
    {
        return;
    }
    setLastMouseWheelMove (juce::Time().currentTimeMillis());
    Slider::mouseWheelMove (e, wheel);
}

//===========================================// Spectrum Display

SpectrumDisplay::SpectrumDisplay (PlayBackEQAudioProcessor& p) : audioProcessor (p),
                                                                 leftPathProducer (audioProcessor.spectrumAnalyzerSampleFifo, p, -80)
{
    setInterceptsMouseClicks (false, false); // Passes mouse events through to component behind it
}

SpectrumDisplay::~SpectrumDisplay() = default;

void SpectrumDisplay::pathProducerProcess (const juce::Rectangle<float> fftBounds)
{
    leftPathProducer.process (fftBounds, audioProcessor.getSampleRate());
}

void SpectrumDisplay::paint (juce::Graphics& g)
{
    // TRACE_COMPONENT();

    using namespace juce;

    // Get FFT path, fill  it and outline it
    const auto FFTpath = leftPathProducer.getPath();
    g.setColour (displayColor);
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

    const auto settings = audioProcessor.getSettings (audioProcessor.apvts);
    const double sampleRate = audioProcessor.getSampleRate();

    // Calculate bass filter coefficients
    const auto bassCoeffs = getBassCoeffs (settings, sampleRate);
    // Convert array to coefficients pointer
    const juce::dsp::IIR::Coefficients<float>::Ptr bassCoeffsPtr (new juce::dsp::IIR::Coefficients<float> (Coefficients (bassCoeffs)));
    // Apply coefficients to filters
    updateCoefficients (bassChain_graphic.get<0>().coefficients, bassCoeffsPtr);
    // Apply bass boost bypass
    bassChain_graphic.setBypassed<0> (!settings.bassBoostEnable);

    // Calculate treble filter coefficients
    const auto trebleCoeffs = getTrebleAttenCoeffs (settings, audioProcessor.getUpSampleRate());
    // Convert array to coefficients pointer
    const juce::dsp::IIR::Coefficients<float>::Ptr trebleCoeffsPtr (new juce::dsp::IIR::Coefficients<float> (Coefficients (trebleCoeffs)));
    // Apply coefficients to filters
    updateCoefficients (trebleChain_graphic.get<0>().coefficients, trebleCoeffsPtr);
    // Apply treble bypass
    trebleChain_graphic.setBypassed<0> (!settings.trebleCutEnable);

    // Get inverse RIAA filter coefficients
    const auto RIAAcoeffs = getRIAACoeffs (audioProcessor.getUpSampleRate());
    // Convert array to coefficients pointer
    const juce::dsp::IIR::Coefficients<float>::Ptr RIAACoeffsPtr (new juce::dsp::IIR::Coefficients<float> (Coefficients (RIAAcoeffs)));
    // Apply coefficients to filters
    updateCoefficients (RIAAChain_graphic.get<0>().coefficients, RIAACoeffsPtr);
    // Apply RIAA bypass
    RIAAChain_graphic.setBypassed<0> (!settings.RIAAEnable);

    // Get coefficients for high-pass filter
    const auto hpfcoeffsPtr = getHighPassCoeffs (settings, sampleRate);

    // apply coefficients to filters depending on slope
    HPChain_graphic.setBypassed<0> (true);
    HPChain_graphic.setBypassed<1> (true);

    switch (settings.highPassSlope)
    {
        case Slope_24:
        {
        }
        case Slope_18:
        {
            updateCoefficients (HPChain_graphic.get<1>().coefficients, hpfcoeffsPtr[1]);
            HPChain_graphic.setBypassed<1> (false);
        }
        case Slope_12:
        {
        }
        case Slope_6: // This will not be Butterworth
        {
            updateCoefficients (HPChain_graphic.get<0>().coefficients, hpfcoeffsPtr[0]);
            HPChain_graphic.setBypassed<0> (false);
        }
        case Off:
        {
        }
    }

    // get low pass coefficients
    const auto lpfcoeffsPtr = getLowPassCoeffs (settings, audioProcessor.getSampleRate());

    // apply coefficients to filters depending on slope
    LPChain_graphic.setBypassed<0> (true);
    LPChain_graphic.setBypassed<1> (true);

    switch (settings.lowPassSlope)
    {
        case Slope_24:
        {
        }
        case Slope_18:
        {
            updateCoefficients (LPChain_graphic.get<1>().coefficients, lpfcoeffsPtr[1]);
            LPChain_graphic.setBypassed<1> (false);
        }
        case Slope_12:
        {
        }
        case Slope_6:
        {
            updateCoefficients (LPChain_graphic.get<0>().coefficients, lpfcoeffsPtr[0]);
            LPChain_graphic.setBypassed<0> (false);
        }
        case Off:
        {
        }
    }

    specialCurveChain_graphic.setBypassed<0> (!settings.specialCurveEnable);
    specialCurveChain_graphic.setBypassed<1> (!settings.specialCurveEnable);
    specialCurveChain_graphic.setBypassed<2> (!settings.specialCurveEnable);
    specialCurveChain_graphic.setBypassed<3> (!settings.specialCurveEnable);

    // Set coefficients for special curve filters
    juce::ReferenceCountedArray<juce::dsp::IIR::Coefficients<float>> specialCurveCoeffs;
    if (settings.specialCurve == SpecialCurve::BBC)
    {
        specialCurveCoeffs = getBBCcoeffs (audioProcessor.getSampleRate());
    }
    else if (settings.specialCurve == SpecialCurve::BBC_XTR_311)
    {
        specialCurveCoeffs = getXTRcoeffs (audioProcessor.getSampleRate());
    }
    updateCoefficients (specialCurveChain_graphic.get<0>().coefficients, specialCurveCoeffs[0]);
    updateCoefficients (specialCurveChain_graphic.get<1>().coefficients, specialCurveCoeffs[1]);
    updateCoefficients (specialCurveChain_graphic.get<2>().coefficients, specialCurveCoeffs[2]);
    updateCoefficients (specialCurveChain_graphic.get<3>().coefficients, specialCurveCoeffs[3]);
}

void SpectrumDisplay::updateResponseCurve()
{
    using namespace juce;

    const auto responseArea = getLocalBounds();
    const auto w = responseArea.getWidth();

    const auto userPreferences = audioProcessor.getUserPreferences();
    const auto sampleRate = audioProcessor.getSampleRate();
    const auto upSampleRate = audioProcessor.getUpSampleRate();

    // This will be a vector of magnitude values for each pixel/frequency that we will conbvert to dB and then turn into a path
    std::vector<double> mags;

    mags.resize (w);

    responseCurve.clear();

    // Create map to convert values in dB to vertical position in responseArea
    const double outputMin = responseArea.getBottom();
    const double outputMax = responseArea.getY();
    auto map = [outputMin, outputMax] (const double input) {
        return jmap (input, -30.0, 30.0, outputMin, outputMax);
    };

    for (int pix = 0; pix < w; ++pix)
    {
        // Magnitude for every frequency is 1 if no filters are modifying it
        double mag = 1.f;

        // Find the frequency asscoiated with this pixel in horizontal dimension of responseArea
        const auto freq = mapToLog10 (static_cast<double> (pix) / static_cast<double> (w), 20.0, 20000.0);

        // Get the magnitude repsonse at this frequency from each filter and multiply them together
        if (!trebleChain_graphic.isBypassed<0>())
        {
            mag *= trebleChain_graphic.get<0>().coefficients->getMagnitudeForFrequency (freq, upSampleRate);
        }

        if (!bassChain_graphic.isBypassed<0>())
        {
            mag *= bassChain_graphic.get<0>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }

        if (!HPChain_graphic.isBypassed<0>())
        {
            mag *= HPChain_graphic.get<0>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }
        if (!HPChain_graphic.isBypassed<1>())
        {
            mag *= HPChain_graphic.get<1>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }

        if (!LPChain_graphic.isBypassed<0>())
        {
            mag *= LPChain_graphic.get<0>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }
        if (!LPChain_graphic.isBypassed<1>())
        {
            mag *= LPChain_graphic.get<1>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }

        if (!RIAAChain_graphic.isBypassed<0>() && userPreferences.showRIAA)
        {
            mag *= RIAAChain_graphic.get<0>().coefficients->getMagnitudeForFrequency (freq, upSampleRate);
        }

        if (!specialCurveChain_graphic.isBypassed<0>())
        {
            mag *= specialCurveChain_graphic.get<0>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }
        if (!specialCurveChain_graphic.isBypassed<1>())
        {
            mag *= specialCurveChain_graphic.get<1>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }
        if (!specialCurveChain_graphic.isBypassed<2>())
        {
            mag *= specialCurveChain_graphic.get<2>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }
        if (!specialCurveChain_graphic.isBypassed<3>())
        {
            mag *= specialCurveChain_graphic.get<3>().coefficients->getMagnitudeForFrequency (freq, sampleRate);
        }

        // Magnitude to dB for this pixel/frequency
        mags[pix] = Decibels::gainToDecibels (mag);

        const auto y = map (mags[pix]);
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

ResponseCurveComponent::ResponseCurveComponent (PlayBackEQAudioProcessor& p, Gui::PresetPanel& pp) : audioProcessor (p), presetPanel (pp), spectrumDisplay (audioProcessor)
{
    // setOpaque (true);
    // For overlaying multiple paths with some weighted transparency - might be useful at some point
    //    FFTpaths.clear();
    //    FFTpaths.resize(numPaths, path);

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

    const auto fftBounds = getAnalysisArea().toFloat();

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
    // const auto width = renderArea.getWidth();

    // Draw frequency tic labels

    // const auto freqs = getFrequencies();
    // const auto xs = getXs (freqs, left, width);
    //
    // for (int i = 0; i < freqs.size(); ++i)
    // {
    //     auto f = freqs[i];
    //     const auto x = xs[i];
    //
    //     bool addK = false;
    //     String str;
    //     if (f > 999.f)
    //     {
    //         addK = true;
    //         f /= 1000.f;
    //     }
    //
    //     str << f;
    //     if (addK)
    //         str << "k";
    //     str << "Hz";
    //
    //     const auto textWidth = g.getCurrentFont().getStringWidth (str);
    //
    //     Rectangle<int> r;
    //
    //     r.setSize (textWidth, fontHeight);
    //     if (f == 20.f && !addK)
    //     {
    //         r.setCentre (x + textWidth / 2.f, 0);
    //     }
    //     else
    //     {
    //         r.setCentre (x, 0);
    //     }
    //     r.setY (fontHeight / 2);
    //
    //     g.drawFittedText (str, r, juce::Justification::centred, 1);
    // }

    // Draw dB tic labels

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
    coordinateBounds.translate (2, 0);
    coordinateComponent.setBounds (coordinateBounds);
    spectrumDisplay.setBounds (getAnalysisArea());

    //    responseCurve.preallocateSpace(getWidth() * 3);
}

//=================================================================
void PathProducer::process (const juce::Rectangle<float> fftBounds, const double sampleRate)
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
            const auto size = tempIncomingBuffer.getNumSamples();

            // Shift every value in monoBuffer starting at index size to the left by size positions
            juce::FloatVectorOperations::copy (monoBuffer.getWritePointer (0, 0), // Destination
                monoBuffer.getReadPointer (0, size), // Source
                monoBuffer.getNumSamples() - size); // Number of values

            // Copy the new buffer from Fifo into the end of monoBuffer
            juce::FloatVectorOperations::copy (monoBuffer.getWritePointer (0, monoBuffer.getNumSamples() - size), // Destination
                tempIncomingBuffer.getReadPointer (0, 0), // Source
                size); // Number of values

            fftDataGenerator.produceFFTDataForRendering (monoBuffer, newFFTOrder);
        }
    }

    const auto fftSize = fftDataGenerator.getFFTSize();
    const auto binWidth = sampleRate / static_cast<double> (fftSize);

    while (fftDataGenerator.getNumAvailableFFTDataBlocks() > 0)
    {
        if (fftDataGenerator.getFFTData (fftData))
        {
            pathGenerator.generatePath (fftData, fftBounds, fftSize, binWidth);
        }
    }

    while (pathGenerator.getNumPathsAvailable() > 0)
    {
        pathGenerator.getPath (FFTPath);
    }
}
