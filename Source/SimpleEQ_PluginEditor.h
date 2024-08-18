/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <juce_gui_basics/juce_gui_basics.h>
#include <sound_meter/sound_meter.h>
#include <sound_meter/sound_meter.h>
#include "PluginProcessor.h"
#include "Gui/PresetPanel.h"
#include "matkat/matkatEditor.h"

class PlayBackEQAudioProcessorEditor
: public juce::AudioProcessorEditor,
private juce::Timer
{
public:
    PlayBackEQAudioProcessorEditor (PlayBackEQAudioProcessor&);
    ~PlayBackEQAudioProcessorEditor() override;
    void paint (juce::Graphics&) override;
    void resized() override;
    void setupSliders ();
    
private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    PlayBackEQAudioProcessor& audioProcessor;
    
    Gui::colorPalette colors;
        
    RotarySliderWithLabels gainSlider,
        balanceSlider,
        bassShelfSlider,
        trebleSlider,
        bassTurnOverSlider,
        highPassFreqSlider,
        lowPassFreqSlider,
        hpSlopeSelect,
        lpSlopeSelect;
            
    ResponseCurveComponent responseCurveComponent;
            
    using APVTS = juce::AudioProcessorValueTreeState;
    using Attachment = APVTS::SliderAttachment;
    
    Attachment gainSliderAttachment,
        balanceSliderAttachment,
        bassShelfSliderAttachment,
        trebleSliderAttachment,
        bassTurnOverSliderAttachment,
        highPassFreqSliderAttachment,
        lowPassFreqSliderAttachment,
        hpSlopeSelectAttachment,
        lpSlopeSelectAttachment;
    
    LookAndFeel lnf;
    
    Gui::AnimatedToggleSlider bassShelfEnableButton, bassBoostEnableButton, trebleAttenEnableButton;
    
    using ButtonAttachment = APVTS::ButtonAttachment;
    
    ButtonAttachment bassShelfEnableButtonAttachment,
                        bassBoostEnableButtonAttachment,
                        trebleAttenEnableButtonAttachment;
                            
    std::vector<juce::Component*> getComps();
    
    sd::SoundMeter::MetersComponent inputMeters, outputMeters;

    void timerCallback() override;
    
    Gui::MaskComponent mask;
    
    Gui::SaveBoxComponent saveBoxComp;
    
    Gui::UnlockBoxComponent unlockBoxComp;
    
    Gui::MessageBoxComponent messageBoxComp, trialMessage;
    
    Gui::PreferencesBoxComponent preferencesBoxComp;

    Gui::PresetPanel presetPanel;
        
    juce::ImageButton saveButton;
    
    juce::ImageComponent logo, logoGrey;
        
    bool firstMeterFlag = true; // Meters segment and clip indicator colors will be updated at the first callback. This ensures that if the plugin loads with bypass on, the meters will be grey
    
    bool firstTrialFlag = true;
    
    ControlSettings prevSettings;
                
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PlayBackEQAudioProcessorEditor)
};
