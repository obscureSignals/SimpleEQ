/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"


PlayBackEQAudioProcessorEditor::PlayBackEQAudioProcessorEditor (PlayBackEQAudioProcessor& p)
    : AudioProcessorEditor (&p), audioProcessor (p),

gainSlider("dB"),
balanceSlider(""),
bassShelfSlider("Hz"),
trebleSlider("kHz"),
bassTurnOverSlider("Hz"),
highPassFreqSlider("Hz"),
lowPassFreqSlider("kHz"),
hpSlopeSelect("dB/Oct"),
lpSlopeSelect("dB/Oct"),

responseCurveComponent(audioProcessor),
gainSliderAttachment(audioProcessor.apvts, "Gain", gainSlider),
balanceSliderAttachment(audioProcessor.apvts, "Balance", balanceSlider),
bassShelfSliderAttachment(audioProcessor.apvts, "Bass Shelf Frequency", bassShelfSlider),
trebleSliderAttachment(audioProcessor.apvts, "Treble Turnover Frequency", trebleSlider),
bassTurnOverSliderAttachment(audioProcessor.apvts, "Bass Turnover Frequency", bassTurnOverSlider),
highPassFreqSliderAttachment(audioProcessor.apvts, "High-Pass Frequency", highPassFreqSlider),
lowPassFreqSliderAttachment(audioProcessor.apvts, "Low-Pass Frequency", lowPassFreqSlider),
hpSlopeSelectAttachment(audioProcessor.apvts, "High-Pass Slope", hpSlopeSelect),
lpSlopeSelectAttachment(audioProcessor.apvts, "Low-Pass Slope", lpSlopeSelect),
bassShelfEnableButtonAttachment(audioProcessor.apvts, "Bass Shelf Enabled", bassShelfEnableButton),
bassBoostEnableButtonAttachment(audioProcessor.apvts, "Bass Boost Enabled", bassBoostEnableButton),
trebleAttenEnableButtonAttachment(audioProcessor.apvts, "Treble Attenuation Enabled", trebleAttenEnableButton),

mask(),

saveBoxComp("untitled"),

unlockBoxComp(),

messageBoxComp(),

preferencesBoxComp(p.getUserPreferences()),

presetPanel(p.getPresetManager(), &mask, &saveBoxComp, &preferencesBoxComp, &unlockBoxComp, &messageBoxComp, p),

saveButton("testImButton")

{
    auto logoImage = juce::ImageCache::getFromMemory(BinaryData::OSlogo_png, BinaryData::OSlogo_pngSize);
    logo.setImage(logoImage);
    
    auto logoGreyImage = juce::ImageCache::getFromMemory(BinaryData::OSlogoGrey_png, BinaryData::OSlogoGrey_pngSize);
    logoGrey.setImage(logoGreyImage);

    hpSlopeSelect.labels.add({0.f, "Off"});
    hpSlopeSelect.labels.add({1.f, "6"});
    hpSlopeSelect.labels.add({2.f, "12"});
    hpSlopeSelect.labels.add({3.f, "18"});
    hpSlopeSelect.labels.add({4.f, "24"});
    hpSlopeSelect.title = "HP\nSlope\n(dB/oct.)";
    hpSlopeSelect.setTextBoxStyle(juce::Slider::TextEntryBoxPosition::TextBoxBelow, 1, 0, 0);
    hpSlopeSelect.setColour(juce::Slider::textBoxTextColourId, colors.textColor);
    hpSlopeSelect.setColour(juce::Slider::thumbColourId, colors.cloud.withBrightness(0.9f));
    
    gainSlider.labels.add({0.f, "-30"});
    gainSlider.labels.add({1.f, "+10"});
    gainSlider.title = "Input\nGain";
    gainSlider.setColour(juce::Slider::textBoxTextColourId, colors.textColor);
    gainSlider.setColour(juce::Slider::thumbColourId, colors.cloud.withBrightness(0.9f));
    
    balanceSlider.labels.add({0.f, "Left"});
    balanceSlider.labels.add({1.f, "Right"});
    balanceSlider.title = "Balance";
    balanceSlider.textFromValueFunction = [](double value)
    {
        if (juce::isWithin<double>(value,0.f,0.001f)) {
            juce::String ctrString = "Center";
            return ctrString;
        }
        else {
            return juce::String(value,2);
        }
    };
    balanceSlider.setColour(juce::Slider::textBoxTextColourId, colors.textColor);
    balanceSlider.setColour(juce::Slider::thumbColourId, colors.cloud.withBrightness(0.9f));
        
    lowPassFreqSlider.labels.add({0.f, "2kHz"});
    lowPassFreqSlider.labels.add({1.f, "22kHz"});
    lowPassFreqSlider.title = "Low\nPass\nFreq.";
    lowPassFreqSlider.setColour(juce::Slider::textBoxTextColourId, colors.textColor);
    lowPassFreqSlider.setColour(juce::Slider::thumbColourId, colors.cloud.withBrightness(0.9f));
    
    lpSlopeSelect.labels.add({0.f, "Off"});
    lpSlopeSelect.labels.add({1.f, "6"});
    lpSlopeSelect.labels.add({2.f, "12"});
    lpSlopeSelect.labels.add({3.f, "18"});
    lpSlopeSelect.labels.add({4.f, "24"});
    lpSlopeSelect.title = "LP\nSlope\n(dB/oct.)";
    lpSlopeSelect.setTextBoxStyle(juce::Slider::TextEntryBoxPosition::TextBoxBelow, 1, 0, 0);
    lpSlopeSelect.setColour(juce::Slider::textBoxTextColourId, colors.textColor);
    lpSlopeSelect.setColour(juce::Slider::thumbColourId, colors.cloud.withBrightness(0.9f));
    
    setupSliders();
    
    for( auto* comp : getComps() )
    {
        addAndMakeVisible(comp);
    }
    
    addChildComponent(trialMessage);

    addChildComponent(saveBoxComp);
    addChildComponent(unlockBoxComp);
    addChildComponent(messageBoxComp);
    addChildComponent(preferencesBoxComp);
    
    addChildComponent(mask);
    
    addChildComponent(logo);
    addChildComponent(logoGrey);

    bassShelfEnableButton.setLookAndFeel(&lnf);
    
    bassBoostEnableButton.setLookAndFeel(&lnf);

    trebleAttenEnableButton.setLookAndFeel(&lnf);
    
    // meters
    
    const int refreshRate_hz = 30;

    // Set the meter OPTONS ...
    sd::SoundMeter::Options meterOptions;
    meterOptions.refreshRate           = refreshRate_hz;  // Frequency of the meter updates (when using the internal timer).
    meterOptions.showPeakHoldIndicator = true;            // Show the peak hold indicator
    meterOptions.tickMarks             = { 0.f, -3.0f, -6.0f, -12.0f};  // Positions (in decibels) of the tick-marks.
    meterOptions.decayTime_ms          = 4000.0f;         // Meter will decay to zero in this this many ms.
    meterOptions.showClipIndicator     = true;
    meterOptions.peakDecayTime_ms      = 3000.0f;         // Peak bars will hold this long in ms.
    meterOptions.clipIndicatorColor    = colors.red;
    
    // Set the meter's segment options ...
    inputMeters.setOptions (meterOptions);
    outputMeters.setOptions (meterOptions);
    
    std::vector<sd::SoundMeter::SegmentOptions> segmentOptions = { { { -60.0f, -18.0f }, { 0.0f, 0.5f }, colors.blue, colors.yellowOrange },
                                                                   { { -18.0f, -3.0f }, { 0.5f, 0.90f }, colors.yellowOrange, colors.orange },
                                                                   { { -3.0f, 0.0f }, { 0.90f, 1.0f }, colors.orange, colors.red } };
    
    inputMeters.setMeterSegments (segmentOptions);
    outputMeters.setMeterSegments (segmentOptions);
    
    auto saveIcon =  juce::ImageCache::getFromMemory(BinaryData::saveIcon_png, BinaryData::saveIcon_pngSize);
    saveButton.setImages(false, true, true, saveIcon,
                         1.0f, {}, saveIcon,
                         1.0f, {}, saveIcon,
                         1.0f, {});
    saveButton.setMouseCursor(juce::MouseCursor::PointingHandCursor);
    
    setSize (651, 597);
    
    // Start the 'polling' timer, that checks the audio processor for peak levels...
    startTimerHz (refreshRate_hz);
}

PlayBackEQAudioProcessorEditor::~PlayBackEQAudioProcessorEditor()
{
}

// The 'polling' timer.
void PlayBackEQAudioProcessorEditor::timerCallback()
{
    if (audioProcessor.getTrialMode() && firstTrialFlag){
        firstTrialFlag = false;
        mask.setVisible(true);
        trialMessage.setBannerText("Welcome to your playBackEQ trial!");
        juce::String trialMessageStr = "Thanks for trying out playBackEQ! You have ";
        int dl;
        if (audioProcessor.isTrialExpired()){
            dl = 0;
        }
        else {
            dl = round(14-audioProcessor.getEven());
        }
        
        trialMessageStr += juce::String(dl);
        trialMessageStr += " day(s) remaining in your trial. If you like it, please consider purchasing a license at obscuresignals.com.";
        trialMessage.setMessageText(trialMessageStr);
        trialMessage.toFront(true);
        // make weak references to components we are passing to lambda - this way we can check if they are null pointers before we use them and avoid trying trying to access deleted objects, should that situation arise
        auto trialMessageWeakReference = juce::WeakReference<Component> (&trialMessage);
        auto maskWeakReference = juce::WeakReference<Component> (&mask);
        trialMessage.enterModalState (false, juce::ModalCallbackFunction::create ([trialMessageWeakReference, maskWeakReference] (int result) {
            if (trialMessageWeakReference != nullptr)
                trialMessageWeakReference->setVisible(false);
            if (maskWeakReference != nullptr)
                maskWeakReference->setVisible(false);
        }));
    }
    
    // Loop through all meters (channels)...
    for (int meterIndex = 0; meterIndex < 2; ++meterIndex)
    {
       // Get the level, of the specified meter (channel), from the audio processor...
        inputMeters.setInputLevel (meterIndex, audioProcessor.getPeakLevelInput (meterIndex));
        outputMeters.setInputLevel (meterIndex, audioProcessor.getPeakLevelOutput (meterIndex));
    }
    inputMeters.refresh(true); // forceRefresh = true makes peak hold dissapear after meter has sat at -inf
    outputMeters.refresh(true); // forceRefresh = true makes peak hold dissapear after meter has sat at -inf
    
    auto settings = audioProcessor.getSettings(audioProcessor.apvts);
    
    if(firstMeterFlag){
        prevSettings.bypassEnable = !(settings.bypassEnable);
        firstMeterFlag = false;
    }
    
    if (settings.bypassEnable != prevSettings.bypassEnable){
        if (settings.bypassEnable){
            responseCurveComponent.setSpectrumDisplayColor(colors.blue);
                        
            const int refreshRate_hz = 30;

            // Set the meter OPTONS ...
            sd::SoundMeter::Options meterOptions;
            meterOptions.refreshRate           = refreshRate_hz;  // Frequency of the meter updates (when using the internal timer).
            meterOptions.showPeakHoldIndicator = true;            // Show the peak hold indicator
            meterOptions.tickMarks             = { 0.f, -3.0f, -6.0f, -12.0f};  // Positions (in decibels) of the tick-marks.
            meterOptions.decayTime_ms          = 4000.0f;         // Meter will decay to zero in this this many ms.
            meterOptions.showClipIndicator     = true;
            meterOptions.peakDecayTime_ms      = 3000.0f;         // Peak bars will hold this long in ms.
            meterOptions.clipIndicatorColor    = colors.red;
            
            // Set the meter's segment options ...
            inputMeters.setOptions (meterOptions);
            outputMeters.setOptions (meterOptions);
            
            std::vector<sd::SoundMeter::SegmentOptions> segmentOptions =
                { { { -60.0f, -18.0f }, { 0.0f, 0.5f }, colors.blue, colors.yellowOrange },
                { { -18.0f, -3.0f }, { 0.5f, 0.90f }, colors.yellowOrange, colors.orange },
                { { -3.0f, 0.0f }, { 0.90f, 1.0f }, colors.orange, colors.red } };
            inputMeters.setMeterSegments (segmentOptions);
            outputMeters.setMeterSegments (segmentOptions);
        }
        else {
            responseCurveComponent.setSpectrumDisplayColor(juce::Colours::whitesmoke.withMultipliedAlpha(0.5));
                        
            const int refreshRate_hz = 30;
            
            auto meterBypassColor = juce::Colours::whitesmoke.withMultipliedAlpha(0.8);

            // Set the meter OPTONS ...
            sd::SoundMeter::Options meterOptions;
            meterOptions.refreshRate           = refreshRate_hz;  // Frequency of the meter updates (when using the internal timer).
            meterOptions.showPeakHoldIndicator = true;            // Show the peak hold indicator
            meterOptions.tickMarks             = { 0.f, -3.0f, -6.0f, -12.0f};  // Positions (in decibels) of the tick-marks.
            meterOptions.decayTime_ms          = 4000.0f;         // Meter will decay to zero in this this many ms.
            meterOptions.showClipIndicator     = true;
            meterOptions.peakDecayTime_ms      = 3000.0f;         // Peak bars will hold this long in ms.
            meterOptions.clipIndicatorColor    = meterBypassColor;
            
            // Set the meter's segment options ...
            inputMeters.setOptions (meterOptions);
            outputMeters.setOptions (meterOptions);
                        
            std::vector<sd::SoundMeter::SegmentOptions> segmentOptions =
            { { { -60.0f, -18.0f }, { 0.0f, 0.5f },meterBypassColor,meterBypassColor},
                { { -18.0f, -3.0f }, { 0.5f, 0.90f },meterBypassColor,meterBypassColor},
                { { -3.0f, 0.0f }, { 0.90f, 1.0f },meterBypassColor,meterBypassColor}};
            inputMeters.setMeterSegments (segmentOptions);
            outputMeters.setMeterSegments (segmentOptions);
            logo.setVisible(false);
            logoGrey.setVisible(true);
        }
    }
    
    logo.setVisible(settings.bypassEnable);
    logoGrey.setVisible(!(settings.bypassEnable));
    
    hpSlopeSelect.setValue(settings.highPassSlope, juce::dontSendNotification);
    hpSlopeSelect.setEnabled(settings.eqEnable);
    if (settings.highPassSlope == Slope::Off){
        highPassFreqSlider.setEnabled(false);
    }
    else{
        highPassFreqSlider.setEnabled(true);
    }
    
    lpSlopeSelect.setEnabled(settings.eqEnable);
    lpSlopeSelect.setValue(settings.lowPassSlope, juce::dontSendNotification);
    if (settings.lowPassSlope == Slope::Off){
        lowPassFreqSlider.setEnabled(false);
    }
    else{
        lowPassFreqSlider.setEnabled(true);
    }

    bassShelfSlider.setEnabled(settings.bassShelfEnable);
    bassShelfEnableButton.setToggleState(settings.bassShelfEnable, juce::dontSendNotification);
    bassShelfEnableButton.setEnabled(settings.eqEnable);
    
    bassTurnOverSlider.setEnabled(settings.bassBoostEnable);
    bassBoostEnableButton.setToggleState(settings.bassBoostEnable, juce::dontSendNotification);
    bassBoostEnableButton.setEnabled(settings.eqEnable);
    
    trebleSlider.setEnabled(settings.trebleCutEnable);
    trebleAttenEnableButton.setToggleState(settings.trebleCutEnable, juce::dontSendNotification);
    trebleAttenEnableButton.setEnabled(settings.eqEnable);

    gainSlider.setEnabled(settings.bypassEnable);
    
    balanceSlider.setEnabled(settings.bypassEnable);
    
    presetPanel.RIAAEnableButton.setToggleState(settings.RIAAEnable, juce::dontSendNotification);
    presetPanel.RIAAEnableButton.setEnabled(settings.bypassEnable);

    presetPanel.delayEnableButton.setToggleState(settings.delayEnable, juce::dontSendNotification);
    presetPanel.delayEnableButton.setEnabled(settings.bypassEnable);
    
    presetPanel.eqBypassEnableButton.setToggleState(settings.eqEnable, juce::dontSendNotification);
    presetPanel.eqBypassEnableButton.setEnabled(settings.bypassEnable);

    presetPanel.latVertButton.setToggleState(settings.vertical, juce::dontSendNotification);
    presetPanel.latVertButton.setEnabled(settings.bypassEnable);
    
    presetPanel.learnDelayButton.setEnabled(settings.bypassEnable);
    
    presetPanel.stereoMonoButton.setToggleState(settings.monoSum, juce::dontSendNotification);
    presetPanel.stereoMonoButton.setEnabled(settings.bypassEnable);
    
    presetPanel.bypassEnableButton.setEnabled(audioProcessor.getUnlocked());

    prevSettings = settings;
    
    // Check if any sliders have had their units changed and need to be updated
    if (preferencesBoxComp.slidersDirty) {
        setupSliders();
        
        trebleSlider.repaint();
        trebleSlider.updateText();
        
//          This is to try and force the host to call the parameter string from value function so that the value in the host changes when user preferences were updated, but it seems like you have to actually change the value to do that - you can't just reset it to what it was before, or nothing happens.
//        audioProcessor.apvts.getParameter("Treble Turnover Frequency")->beginChangeGesture();
//        audioProcessor.apvts.getParameter("Treble Turnover Frequency")->endChangeGesture();
//        audioProcessor.apvts.getParameter("Treble Turnover Frequency")->setValueNotifyingHost(audioProcessor.apvts.getParameter("Treble Turnover Frequency")->getValue());
        
        bassTurnOverSlider.repaint();
        bassTurnOverSlider.updateText();
        
        bassShelfSlider.repaint();
        bassShelfSlider.updateText();

        highPassFreqSlider.repaint();
        highPassFreqSlider.updateText();

        preferencesBoxComp.slidersDirty = false;
    }
    
    if (preferencesBoxComp.responseDirty) {
        responseCurveComponent.updateReponse();
        preferencesBoxComp.responseDirty = false;
    }
    
}

void PlayBackEQAudioProcessorEditor::paint (juce::Graphics& g)
{
    using namespace juce;
    g.fillAll (Colours::black);
    
    g.setFont(13);
    g.setColour(colors.textColor);

    String inStr = "In";
    Rectangle<int> rIn;

    auto h = g.getCurrentFont().getHeight();
    rIn.setSize(inputMeters.getWidth(), h);
    rIn.setCentre((int)(inputMeters.getX() + inputMeters.getWidth()/2), (int)(inputMeters.getBottom()+h/2+1));

    g.drawFittedText(inStr, rIn, juce::Justification::centred, 1);
    
    String outStr = "Out";
    Rectangle<int> rOut;

    rOut.setSize(outputMeters.getWidth(), h);
    rOut.setCentre((int)(outputMeters.getX() + outputMeters.getWidth()/2), (int)(outputMeters.getBottom()+h/2+1));
        
    g.drawFittedText(outStr, rOut, juce::Justification::centred, 1);
    
    g.setColour(juce::Colours::silver);
    
}

void PlayBackEQAudioProcessorEditor::resized()
{
    auto totalWidth = getLocalBounds().getWidth();

    presetPanel.setBounds(juce::Rectangle<int>(0, 0, totalWidth, presetPanelHeight));
    
    mask.setBounds(getLocalBounds()); // Covers entire window
    
    saveBoxComp.setSize(250, 150);
    saveBoxComp.setCentrePosition(getLocalBounds().getCentre());
    
    unlockBoxComp.setSize(340, 300);
    unlockBoxComp.setCentrePosition(getLocalBounds().getCentre());
    
    messageBoxComp.setSize(250, 160);
    messageBoxComp.setCentrePosition(getLocalBounds().getCentre());
    
    trialMessage.setSize(250, 220);
    trialMessage.setCentrePosition(getLocalBounds().getCentre());

    preferencesBoxComp.setSize(350, 260);
    preferencesBoxComp.setCentrePosition(getLocalBounds().getCentre());

    responseCurveComponent.setBounds(juce::Rectangle<int>(33, presetPanelHeight, 631-33, 212));
    inputMeters.setBounds(juce::Rectangle<int>(6, 22+presetPanelHeight, 27, 186));
    outputMeters.setBounds(juce::Rectangle<int>(618, 22+presetPanelHeight, 27, 186));

    int yTopRowSlider = 198+presetPanelHeight;
    int wSlider = 130;
    int hSlider = 156;
    
    // top row controls
    highPassFreqSlider.setBounds(juce::Rectangle<int>(0, yTopRowSlider, wSlider, hSlider));
    bassShelfSlider.setBounds(juce::Rectangle<int>(130, yTopRowSlider, wSlider, hSlider));
    bassTurnOverSlider.setBounds(juce::Rectangle<int>(260, yTopRowSlider, wSlider, hSlider));
    
    // treble controls
    trebleSlider.setBounds(juce::Rectangle<int>(390, yTopRowSlider, wSlider, hSlider));
        lowPassFreqSlider.setBounds(juce::Rectangle<int>(520, yTopRowSlider, wSlider, hSlider));
    
    int yBottomRowSlider = 364+presetPanelHeight;
    // bottom row controls
    hpSlopeSelect.setBounds(juce::Rectangle<int>(0, yBottomRowSlider, wSlider, hSlider));
    gainSlider.setBounds(juce::Rectangle<int>(130, yBottomRowSlider, wSlider, hSlider));
    balanceSlider.setBounds(juce::Rectangle<int>(260, yBottomRowSlider, wSlider, hSlider));
    
    int logoDim = 108;
    int logoCenterX = 455;
    int logoCenterY = (int)(yBottomRowSlider+hSlider/2)+10;
    juce::Rectangle<int> rLogo;
    rLogo.setHeight(logoDim);
    rLogo.setWidth(logoDim);
    rLogo.setCentre(logoCenterX, logoCenterY);
    logo.setBounds(rLogo);
    logoGrey.setBounds(rLogo);

    lpSlopeSelect.setBounds(juce::Rectangle<int>(520, yBottomRowSlider, wSlider, hSlider));
    
    int wButton = 36;
    int hButton = 18;
    int yButton = 354+presetPanelHeight;
    
    // buttons
    bassShelfEnableButton.setBounds(juce::Rectangle<int>(177, yButton, wButton, hButton));
    bassBoostEnableButton.setBounds(juce::Rectangle<int>(307, yButton, wButton, hButton));
    trebleAttenEnableButton.setBounds(juce::Rectangle<int>(437, yButton, wButton, hButton));
}

std::vector<juce::Component*> PlayBackEQAudioProcessorEditor::getComps()
{
    return
    {
        &gainSlider,
        &balanceSlider,
        &bassShelfSlider,
        &bassTurnOverSlider,
        &trebleSlider,
        &highPassFreqSlider,
        &lowPassFreqSlider,
        
        &responseCurveComponent,
        
        &hpSlopeSelect,
        &lpSlopeSelect,
        
        &bassShelfEnableButton,
        &bassBoostEnableButton,
        &trebleAttenEnableButton,

        &inputMeters,
        &outputMeters,
        
        &presetPanel
    };
}

void PlayBackEQAudioProcessorEditor::setupSliders ()
{
    auto userPreferences = audioProcessor.getUserPreferences();
    
    trebleSlider.labels.clear();
    switch( userPreferences.trebleUnits )
    {
        case UnitSelect::freq: {
            trebleSlider.labels.add({0.f, "1kHz"});
            trebleSlider.labels.add({1.f, "10kHz"});
            trebleSlider.title = "Treble\nTurnover\nFreq.";
            trebleSlider.setTextValueSuffix("Hz");
            break;
        }
        case UnitSelect::TC: {
            trebleSlider.labels.add({0.f, juce::CharPointer_UTF8("160µs")});
            trebleSlider.labels.add({1.f, juce::CharPointer_UTF8("8µs")});
            trebleSlider.title = "Treble\nTime\nConstant";
            trebleSlider.setTextValueSuffix(juce::CharPointer_UTF8("µs"));
            break;
        }
        case UnitSelect::atten: {
            trebleSlider.labels.add({0.f, "-20dB"});
            trebleSlider.labels.add({1.f, "-1dB"});
            trebleSlider.title = "Atten.\nat\n10kHz";
            trebleSlider.setTextValueSuffix("dB");
            break;
        }
    }
    trebleSlider.setColour(juce::Slider::textBoxTextColourId, colors.textColor);
    trebleSlider.setColour(juce::Slider::thumbColourId, colors.cloud.withBrightness(0.9f));
    
    bassTurnOverSlider.labels.clear();
    switch( userPreferences.bassBoostUnits )
    {
        case UnitSelect::freq: {
            bassTurnOverSlider.labels.add({0.f, "100"});
            bassTurnOverSlider.labels.add({1.f, "1kHz"});
            bassTurnOverSlider.title = "Bass\nTurnover\nFreq.";
            bassTurnOverSlider.setTextValueSuffix("Hz");
            break;
        }
        case UnitSelect::TC: {
            bassTurnOverSlider.labels.add({0.f, juce::CharPointer_UTF8("1.6ms")});
            bassTurnOverSlider.labels.add({1.f, juce::CharPointer_UTF8("160µs")});
            bassTurnOverSlider.title = "Bass\nBoost\nTime\nConst.";
            bassTurnOverSlider.setTextValueSuffix(juce::CharPointer_UTF8("µs"));
            break;
        }
        case UnitSelect::atten: {
            break;
        }
    }
    bassTurnOverSlider.setColour(juce::Slider::textBoxTextColourId, colors.textColor);
    bassTurnOverSlider.setColour(juce::Slider::thumbColourId, colors.cloud.withBrightness(0.9f));
    
    bassShelfSlider.labels.clear();
    switch( userPreferences.bassShelfUnits )
    {
        case UnitSelect::freq: {
            bassShelfSlider.labels.add({0.f, "10Hz"});
            bassShelfSlider.labels.add({1.f, "110Hz"});
            bassShelfSlider.title = "Bass\nShelving\n Freq.";
            bassShelfSlider.setTextValueSuffix("Hz");
            break;
        }
        case UnitSelect::TC: {
            bassShelfSlider.labels.add({0.f, "16ms"});
            bassShelfSlider.labels.add({1.f, "1.4ms"});
            bassShelfSlider.title = "Bass\nShelving\nTime\nConst.";
            bassShelfSlider.setTextValueSuffix(juce::CharPointer_UTF8("µs"));
            break;
        }
        case UnitSelect::atten: {
            break;
        }
    }
    bassShelfSlider.setColour(juce::Slider::textBoxTextColourId, colors.textColor);
    bassShelfSlider.setColour(juce::Slider::thumbColourId, colors.cloud.withBrightness(0.9f));
    
    highPassFreqSlider.labels.clear();
    switch( userPreferences.highPassUnits )
    {
        case UnitSelect::freq: {
            highPassFreqSlider.labels.add({0.f, "10Hz"});
            highPassFreqSlider.labels.add({1.f, "400Hz"});
            highPassFreqSlider.title = "High\nPass\nFreq.";
            highPassFreqSlider.setTextValueSuffix("Hz");
            break;
        }
        case UnitSelect::TC: {
            highPassFreqSlider.labels.add({0.f, "16ms"});
            highPassFreqSlider.labels.add({1.f, juce::CharPointer_UTF8("400µs")});
            highPassFreqSlider.title = "High\nPass\nTime\nConst.";
            highPassFreqSlider.setTextValueSuffix(juce::CharPointer_UTF8("µs"));
            break;
        }
        case UnitSelect::atten: {
            break;
        }
    }
    highPassFreqSlider.setColour(juce::Slider::textBoxTextColourId, colors.textColor);
    highPassFreqSlider.setColour(juce::Slider::thumbColourId, colors.cloud.withBrightness(0.9f));
};
