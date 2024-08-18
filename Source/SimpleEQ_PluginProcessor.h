/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>
#include "Service/PresetManager.h"
#include "Service/GetSupportFilesLoc.h"
#include "matkat/matkatProcessor.h"
#include "OfflineUnlockStatus.h"

enum Slope
{
    Off,
    Slope_6,
    Slope_12,
    Slope_18,
    Slope_24
};

enum UnitSelect
{
    freq,
    TC,
    atten
};

enum OvrSampAmt
{
    none,
    minimum,
    low,
    medium,
    high,
    maximum
};

//==============================================================================
class MyOfflineUnlockStatus  : public juce::OfflineUnlockStatus
{
public:
    MyOfflineUnlockStatus() = default;

    juce::String getProductID() override
    {
        return ProjectInfo::projectName;
    }
    
    bool doesProductIDMatch (const juce::String& returnedIDFromServer) override
    {
        return getProductID() == returnedIDFromServer;
    }
    
    juce::RSAKey getPublicKey() override
    {
        return juce::RSAKey ("3,98708d2b169b77311fcaf1909f7ffe0fcabf6a5a8e323a9b12a79b9b65e9141d");
    }

    juce::String getWebsiteName() override
    {
        return "obscuresignals.com";
    }
};

//==============================================================================
/** Replaces left and right channels with sum of both / 2 */
template <typename SampleType>
class monoSum {
public:
    template <typename ProcessContext>
    void process (const ProcessContext& contextLeft, const ProcessContext& contextRight) {
        const auto& inputBlockLeft = contextLeft.getInputBlock();
        auto& outputBlockLeft = contextLeft.getOutputBlock();
        
        const auto& inputBlockRight = contextRight.getInputBlock();
        auto& outputBlockRight = contextRight.getOutputBlock();

        const auto numInputChannelsLeft  = inputBlockLeft.getNumChannels();
        const auto numInputChannelsRight  = inputBlockRight.getNumChannels();
        
        jassert (numInputChannelsLeft == 1);
        jassert (numInputChannelsRight == 1);

        const auto numOutputChannelsLeft = outputBlockLeft.getNumChannels();
        
        jassert (numOutputChannelsLeft == 1);
        
        const auto numOutputChannelsRight= outputBlockRight.getNumChannels();
        
        jassert (numOutputChannelsRight == 1);
        
        [[maybe_unused]] const auto numOutSamplesLeft = outputBlockLeft.getNumSamples();

        jassert (inputBlockLeft.getNumSamples() == numOutSamplesLeft);
        
        [[maybe_unused]] const auto numOutSamplesRight = outputBlockRight.getNumSamples();

        jassert (inputBlockRight.getNumSamples() == numOutSamplesRight);

        outputBlockLeft.replaceWithSumOf(inputBlockLeft, inputBlockRight);
        outputBlockLeft.multiplyBy(0.5f);
        
        outputBlockRight.copyFrom(outputBlockLeft);
    }
};

//==============================================================================
/**
When learning is enabled, takes in left and right frames of audio and estimates the delay between them in samples. Each delay estimate is collected anda final estimation is determined after the data is cleaned and averaged.
 
When learning is enabled, and a frame is delivered, the process is as follows:
1. Check if both channels of frame meet static RMS level threshold - ignore otherwise - We do not want to attemt to estimate the delay if we are delivered a signal with no or little information.
2. Collect frames until a larger block size is reached - estimations are more reliable for longer streches of audio, so we will gather a few frames before we try
3. Window the collected data
4. Attemt an estimation
    a. Perform an FFT on the windowed data
    b. In the frequency domain, we are going to multiply the left channel with the complex conjugate of the right channel (and, if whitening is enabled, divide the product by the magnitude of itself in order to set all the magnitudes to 1). We are only going to do this for frequencies between lowCutOffFreq and highCutOffFreq - everything else will be zeros. Note that whitening is disabled by default, so there is no division by the magnitude of product. The frequency-domain data will have 5 sectors:
        Sector 0 = Y[0]:Y[startIdxBin-1] = 0, these are positive frequencies below lowCutOffFreq
        Sector 1 = Y[startIdxBin]:Y[stopIdxBin-1] = ( fft(xl) • conj(fft(xr)) ) / |fft(xl) • conj(fft(xr))|, these are postive frequencies between lowCutOffFreq and highCutOffFreq
        Sector 2 = Y[stopIdxBin-1]:Y[bigBlockSize•2 - stopIdxBin] = 0, these are postive frequencies above highCutOffFreq and negative frequencies below -1•highCutOffFreq
        Sector 3 = Y[bigBlockSize•2 - stopIdxBin + 2]:Y[bigBlockSize•2 - startIdxBin] = ( fft(xl) • conj(fft(xr)) ) / |fft(xl) • conj(fft(xr))|, these are negative frequencies between -1•highCutOffFreq and -1•lowCutOffFreq
        Sector 4 = Y[bigBlockSize•2 - startIdxBin + 2] = 0, these are negative frequencies above -1 • lowCutOffFreq
    c. Perform an inversve fft on Y, yielding y
    d. Find the maximum value in y and the delay to which it corresponds
    e. If this delay is less than the maximum value allowed, store it - for this application we not expecting delays larger than 1ms - if we encounter an estimation larger than that, it is likely erroneous
5. While learning is taking place, we report out a delay value that is a weighted average of the delays accumulated since learning began, where the weights are the max value of the cross correlation for that delay rasied to weightExponent - Higher max values mean that the alignment at that delay was more complete - Raising this value to weightExponent icreases the influence of this weighting
6. When learning is stopped the collected data is cleaned and a final estimation is returned
    a.  Remove delay values for blocks where the RMS level of either channel is more than one standard deviation below the mean RMS level of that channel - thowing out  lower-level/unreliable frames
    b.  Remove delay values that are more than one standard deviation away from the mean delay value - removing outliers
    c.  Retrun weighted average of remaining delay estimates, where the the weights are the maximum value in the cross-correlation raised to weightExponent
*/

// TODO: This all needs to be looked at for speed!!!!!!!!
template<typename BlockType>
struct GCCPHAT
{
    // Constructor:
    GCCPHAT()
    {
        prepared.set(false);
    }

    void update(const juce::dsp::AudioBlock<BlockType> leftBlock, const juce::dsp::AudioBlock<BlockType> rightBlock)
    {
        jassert(prepared.get());

        if (!learningFlag){ // We only need to do update if we are in learning mode.
            return;
        }
        
        // Get input RMS level for each channel
        const auto RMSCurrL = getRMS_block(leftBlock);
        const auto RMSCurrR = getRMS_block(rightBlock);
        
        // Convert to DB to check against static level threshold
        const auto RMSdBLCurr =  juce::Decibels::gainToDecibels(RMSCurrL);
        const auto RMSdBRCurr =  juce::Decibels::gainToDecibels(RMSCurrR);
        
        // If static level threshold not met by either channel return
        if (RMSdBLCurr < staticLevelThresh || RMSdBRCurr < staticLevelThresh){ // Do not add blocks if either one is below the static level threshold
            return;
        }
        
        // Add RMS value of each block in big block to a rms value array
        RMS_Larray.add(RMSCurrL);
        RMS_Rarray.add(RMSCurrR);

        const int startIdxFrame = frameCount*blockSize; // The location at which the blocks need to be copied into xl and xr depends on frame count
        
        // Copy buffer data to local arrays.
        for( int idx = 0; idx < blockSize; idx++ ) {
            xl.setUnchecked(startIdxFrame+idx, leftBlock.getSample(0, idx));
            xr.setUnchecked(startIdxFrame+idx, rightBlock.getSample(0, idx));
        }
        
        // Collect frames until the numFrames have been accumulated - don't process until then
        frameCount++;
        if (frameCount == numFrames) {
            frameCount = 0;
        }
        else {
            return;
        }
        
        // The first half of xl and xr should now be filled with real input data. The second half is to accomodate the imaginary parts after the DFT

        
        // Get RMS value of big block = (RMS(RMS(each small block)))
        const auto RMS_L = getRMS_array(RMS_Larray);
        const auto RMS_R = getRMS_array(RMS_Rarray);
    
        RMS_Larray.clear();
        RMS_Rarray.clear();
        
        // Time domain windowing
        window->multiplyWithWindowingTable(xl.data(), bigBlockSize);
        window->multiplyWithWindowingTable(xr.data(), bigBlockSize);
                
        // References for clarity, capitals are freq domain
        auto &Xl = xl;
        auto &Xr = xr;
        
        // TODO: If there was any way to calculate only the frequencies from lowCutOffFreq to highCutOffFreq and still take advantage of the FFT algorithm/optimized FFT engines, that would be epic
        // Perform FFTs.
        // These arrays should now be completely filled with interleaved real and imaginary data.
        // For example, the first complex number in the left array, a + jb, is Xl[0] + j*Xl[1]
        fft->performRealOnlyForwardTransform (Xl.data(),false);
        fft->performRealOnlyForwardTransform (Xr.data(),false);

        
        /** conjProd = Xl•conj(Xr)
        (a + jb) (c + jd) = (ac - bd) + j(ad + bc).
        (a + jb) (c - jd) = (ac + bd) + j(bc - ad) */
        auto &conjProd = xl; // Product with complex conjugate, reference for clarity
        
        /** magConjProd = |conjProd|
        |(a + ib)| = (a^2 + b^2)^0.5
        a = conjProd[idx]
        b = conjProd[idx+1]
        Starting at magConjProd[0], every even index is |conjProd| - every odd index is unused */
        auto &magConjProd = xr; // Magnitude of product with complex conjugate, reference for clarity
        
        /** Y = conjProd/magConjProd */
        auto &Y = xl; // Product with complex conjugate divided by magnitude of product with complex conjugate, reference for clarity

        /** Sector 0 = Y[0]:Y[startIdxBin-1] = 0, these are positive frequencies below lowCutOffFreq */
        for( int idx = 0; idx < startIdxBin; idx++ ) {
            conjProd.setUnchecked(idx, 0.f);
        }
        
        /** Sector 1 = Y[startIdxBin]:Y[stopIdxBin-1] = ( fft(xl) • conj(fft(xr)) ) / |fft(xl) • conj(fft(xr))|, these are postive frequencies between lowCutOffFreq and highCutOffFreq */
        for( int idx = startIdxBin; idx < stopIdxBin; idx+=2 ) {
            // Capture values because we are going to alter these values in the array before we are done using them
            auto a = Xl[idx];
            auto b = Xl[idx+1];
            auto c = Xr[idx];
            auto d = Xr[idx+1];
            conjProd.setUnchecked(idx, a*c + b*d);
            conjProd.setUnchecked(idx+1, b*c - a*d);
            if (whiten){
                magConjProd.setUnchecked(idx, sqrt(powf(conjProd[idx],2) + powf(conjProd[idx+1],2)));
                Y.setUnchecked(idx, conjProd[idx]/magConjProd[idx]);
                Y.setUnchecked(idx+1, conjProd[idx+1]/magConjProd[idx]);

            }
            else {
                Y.setUnchecked(idx, conjProd[idx]);
                Y.setUnchecked(idx+1, conjProd[idx+1]);
            }
        }
        
        /** Sector 2 = Y[stopIdxBin-1]:Y[bigBlockSize•2 - stopIdxBin] = 0, these are postive frequencies above highCutOffFreq and negative frequencies below -1•highCutOffFreq */
        for( int idx = stopIdxBin; idx < Xl.size() - stopIdxBin + 2; idx++ ) {
            conjProd.setUnchecked(idx, 0.f);
        }
        
        /** Sector 3 = Y[bigBlockSize•2 - stopIdxBin + 2]:Y[bigBlockSize•2 - startIdxBin] = ( fft(xl) • conj(fft(xr)) ) / |fft(xl) • conj(fft(xr))|, these are negative frequencies between -1•highCutOffFreq and -1*lowCutOffFreq */
        for( int idx = Xl.size() - stopIdxBin + 2; idx <  Xl.size() - startIdxBin + 2; idx+=2 ) {
            // Capture values because we are going to alter these values in the array before we are done using them
            auto a = Xl[idx];
            auto b = Xl[idx+1];
            auto c = Xr[idx];
            auto d = Xr[idx+1];
            conjProd.setUnchecked(idx, a*c + b*d);
            conjProd.setUnchecked(idx+1, b*c - a*d);
            if (whiten){
                magConjProd.setUnchecked(idx, sqrt(powf(conjProd[idx],2) + powf(conjProd[idx+1],2)));
                Y.setUnchecked(idx, conjProd[idx]/magConjProd[idx]);
                Y.setUnchecked(idx+1, conjProd[idx+1]/magConjProd[idx]);
            }
            else {
                Y.setUnchecked(idx, conjProd[idx]);
                Y.setUnchecked(idx+1, conjProd[idx+1]);
            }

        }
        
        /** Sector 4 = Y[bigBlockSize•2 - startIdxBin + 2] = zeros, these are negative frequencies above -1•lowCutOffFreq */
        for( int idx = Xl.size() - startIdxBin + 2; idx < Xl.size(); idx++ ) {
            conjProd.setUnchecked(idx, 0.f);
        }
        
        auto &y = xl; // Time domain cross-correlation output, reference for clarity
        /** y = ifft(Y) */
        fft->performRealOnlyInverseTransform (y.data());
        
        auto maxIdx = findIdxOfMax(y, bigBlockSize);
        float maxVal = y[maxIdx];
        int currentDelay = idx2Delay(maxIdx);
        
        if (std::abs(currentDelay) <= delayLengthThreshSamples){ // Only set delay if delay time threshold is met
            // Add the big block RMS values to the level arrays
            leftLevels.add(RMS_L);
            rightLevels.add(RMS_R);
            
            // Add the big block RMS values to a stats accumulator - we will need the mean and standard deviation for post processing when learning is stopped
            levelStatsLeft.addValue(RMS_L);
            levelStatsRight.addValue(RMS_R);
            
            setDelay(currentDelay, maxVal);
        }
    }

    void prepare(int newBlockSize, int sampleRate, size_t upExp) // This should get called in prepareToPlay
    {
        prepared.set(false);
        
        float upSampleRate = sampleRate*pow(2,upExp);
        
        numFrames = std::ceil(2048/newBlockSize); // We want bigBlockSize to be >= 2048 TODO: should this be a time instead of a number of samples?
        
        blockSize = newBlockSize*pow(2,upExp); // BlockSize at upsample rate
        bigBlockSize = blockSize*numFrames; // Blocksize after accumulation over numFrames frames
        
        fftOrder = log2(bigBlockSize);
        fft = std::make_unique<juce::dsp::FFT>(fftOrder);
        window = std::make_unique<juce::dsp::WindowingFunction<float>>(bigBlockSize, juce::dsp::WindowingFunction<float>::hann, false);
        
        xl.clear();
        xr.clear();
        
        RMS_Larray.ensureStorageAllocated(numFrames);
        RMS_Rarray.ensureStorageAllocated(numFrames);
        maxNumDelays = ceil((maxLearnTime*upSampleRate)/bigBlockSize);
        delays.ensureStorageAllocated(maxNumDelays);
        rightLevels.ensureStorageAllocated(maxNumDelays);
        leftLevels.ensureStorageAllocated(maxNumDelays);
        weights.ensureStorageAllocated(maxNumDelays);
        
        // Set input array sizes
        // These need to be 2*blocksize*numFrames so the imaginary parts have somewhere to go
        xl.resize(bigBlockSize*2); // xl and xr is the only data being claimed here, but they're 2 * numFrames * the blocksize
        xr.resize(bigBlockSize*2);
                
        startIdxBin = std::ceil((lowCutOffFreq/sampleRate)*newBlockSize*numFrames)*2; // The bin corresponding to lowCutOffFreq
        stopIdxBin = std::ceil((highCutOffFreq/sampleRate)*newBlockSize*numFrames)*2; // The bin corresponding to highCutOffFreq
        
        delayLengthThreshSamples = ceil(delayLengthThreshTime*upSampleRate);
        
        frameCount = 0;
        
        prepared.set(true);
    }
    void startLearning(){
        RMS_Larray.clear();
        RMS_Rarray.clear();
        xl.clear();
        xr.clear();
        xl.resize(bigBlockSize*2);
        xr.resize(bigBlockSize*2);
        frameCount = 0;
        delays.clear();
        leftLevels.clear();
        rightLevels.clear();
        levelStatsLeft.reset();
        levelStatsRight.reset();
        weights.clear();
        weightedDelaySum = 0;
        weightsSum = 0;
        delayCount = 0;
        learningFlag = true;
    }
    
    void stopLearning(){
        learningFlag = false;
        
        // Remove delays that don't meet level threshold
        juce::Array<float> delaysLevelCleaned;
        juce::Array<float> weightsLevelCleaned;
        juce::Array<float> delaysStdDevCleaned;
        juce::Array<float> weightsStdDevCleaned;
        
        delaysLevelCleaned.ensureStorageAllocated(delays.size());
        weightsLevelCleaned.ensureStorageAllocated(delays.size());
        delaysStdDevCleaned.ensureStorageAllocated(delays.size());
        weightsStdDevCleaned.ensureStorageAllocated(delays.size());
        
        const float RMSLeftTotal = getRMS_array(leftLevels);
        const float RMSRightTotal = getRMS_array(rightLevels);
        const float stdDevLevelLeft = levelStatsLeft.getStandardDeviation();
        const float stdDevLevelRight = levelStatsRight.getStandardDeviation();
        const float levelThreshLeft = RMSLeftTotal-stdDevLevelLeft;
        const float levelThreshRight = RMSRightTotal-stdDevLevelRight;
        
        juce::StatisticsAccumulator<float> delaysLevelCleanedStats;

        // Remove delay values for frames where the RMS level of either channel is more than one standard deviation below the mean RMS level of that channel
        int numLevelReject = 0;
        for (int idx = 0; idx < delays.size(); idx++){
            if (leftLevels[idx] > levelThreshLeft && rightLevels[idx] > levelThreshRight){
                delaysLevelCleaned.add(delays[idx]);
                delaysLevelCleanedStats.addValue(delays[idx]);
                weightsLevelCleaned.add(weights[idx]);
            }
            else{
                numLevelReject++;
            }
        }
        
        int numStdDevReject = 0;
        const float stdDevDelay = delaysLevelCleanedStats.getStandardDeviation();
        const float delayMean = delaysLevelCleanedStats.getAverage();
        const float delayThreshLow = delayMean-stdDevDelay;
        const float delayThreshHigh = delayMean+stdDevDelay;
        
        // Remove delay values that are more than one standard deviation away from the mean delay value
        for (int idx = 0; idx < delaysLevelCleaned.size(); idx++){
            if (delaysLevelCleaned[idx] >= delayThreshLow && delaysLevelCleaned[idx] <= delayThreshHigh){
                delaysStdDevCleaned.add(delaysLevelCleaned[idx]);
                weightsStdDevCleaned.add(weightsLevelCleaned[idx]);
            }
            else{
                numStdDevReject++;
            }
        }
        
        // Weight the remaining delay values by the maximum value in the cross-correlation raised to weightExponent
        float weightedDelaySumCleaned = 0;
        float weightsSumCleaned = 0;
        for (int idx = 0; idx < delaysStdDevCleaned.size(); idx++){
            weightedDelaySumCleaned+= weightsStdDevCleaned[idx]*delaysStdDevCleaned[idx];
            weightsSumCleaned+= weightsStdDevCleaned[idx];
        }
        delay = weightedDelaySumCleaned/weightsSumCleaned;
        if (isnan(delay)){
            delay = 0.f;
        }
        newDelayValueAvailable = true;
    }
    
    float getDelay(bool resetNewDelayValueAvailable){
        if (resetNewDelayValueAvailable) {
            newDelayValueAvailable = false;
        }
        return delay;
    }
    
    int getMaxDelay(){
        return delayLengthThreshSamples;
    }
    
    bool getNewDelayValueAvailable(){
        return newDelayValueAvailable;
    }
        
    bool islearning(){
        return learningFlag;
    }
    
    int getFFTSize() const { return 1 << fftOrder; }
    bool isPrepared() const { return prepared.get(); }
    
private:
    void setDelay(int newDelay, float maxValue)
    {
        delays.add(newDelay);
        delayCount++;
        const auto weight = powf(maxValue,weightExponent);
        weights.add(weight);
        weightedDelaySum+= weight*newDelay;
        weightsSum+= weight;
        delay = weightedDelaySum/weightsSum;
        newDelayValueAvailable = true;
        if (delayCount == maxNumDelays){
            stopLearning();
        }
    }

    float getRMS_block (const juce::dsp::AudioBlock<BlockType> block)
    {
        /** Returns the rms value of a mono AudioBlock*/
        double sum = 0.0;

        for( int idx = 0; idx < blockSize; idx++ ) {
            auto sample = block.getSample(0, idx);
            sum += sample * sample;
        }
        return std::sqrt (sum / blockSize);
    }
    
    float getRMS_array (const juce::Array<float> array)
    {
        /** Returns the rms value of  a juce::array*/
        
        double sum = 0.0;
        for( int idx = 0; idx < array.size(); idx++ ) {
            sum += array[idx] * array[idx];
        }
        return std::sqrt (sum / array.size());
    }
    
    
    int findIdxOfMin (const juce::Array<float> data, int numValues)
    {
        /** Returns the index of the minimum value in a juce array. */
        int currentMinIdx = 0;
        auto currentMinVal = data[0];
        for (int idx = 1; idx < numValues; idx++)
        {
            if (data[idx] < currentMinVal)
            {
                currentMinIdx = idx;
                currentMinVal = data[idx];
            }
        }
        return currentMinIdx;
    }
    
    int findIdxOfMax (juce::Array<float> data, int numValues)
    {
        /** Returns the index of the maximum value in a juce array. */
        int currentMaxIdx = 0;
        auto currentMaxVal = data[0];
        for (int idx = 1; idx < numValues; idx++)
        {
            if (data[idx] > currentMaxVal)
            {
                currentMaxIdx = idx;
                currentMaxVal = data[idx];
            }
        }
        return currentMaxIdx;
    }
    
    int idx2Delay (int idx)
    {
        /** Converts an index of the cross correlation to a delay value in samples */
        int delay;
        
        if (idx <= bigBlockSize/2-1){
            delay = idx;
        }
        else if (idx == bigBlockSize/2){
            delay = 0;
        }
        else {
            delay = idx - bigBlockSize;
        }
        return delay;
    }
    
    int blockSize; // Block size at upsample rate
    int bigBlockSize; // Block size after numframes accumulations
    int fftOrder;
    std::unique_ptr<juce::dsp::FFT> fft;
    std::unique_ptr<juce::dsp::WindowingFunction<float>> window;
    juce::Array<float> xl, xr; // Left and right channel arrays
    juce::Atomic<bool> prepared = false;
    
    float delay { 0 }; // Delay value in samples to be consumed externally
    int learningFlag { 0 }; // The algorithm is currently learning
    
    int numFrames { 1 }; // How many frames should be accumulated before processing
    int frameCount { 0 }; // How many frames have been accumulated so far
    
    const float delayLengthThreshTime { 0.0009999 } ; // Maximum delay in seconds - larger values are rejected
    int delayLengthThreshSamples ; // Maximum delay in samples - larger values are rejected
    
    float lowCutOffFreq { 80 }; // Low cutoff for normalization in Hz
    float highCutOffFreq { 4000 }; // High cutoff for normalization in Hz
    int startIdxBin; // Index corrosponding to lowCutOffFreq Hz
    int stopIdxBin; // Index corrosponding to highCutOffFreq Hz
    
    int delayCount { 0 }; // How many delays have been set during learning
    
    juce::Array<float> delays; // Calculated delays before cleaning
    juce::Array<float> leftLevels; // RMS level for each big block
    juce::Array<float> rightLevels; // RMS level for each big block
    juce::StatisticsAccumulator<float> levelStatsLeft;
    juce::StatisticsAccumulator<float> levelStatsRight;
    int maxLearnTime { 300 }; // Maximum time to learn delay in seconds
    int maxNumDelays; // Maximum number of delays before learning is stopped automatically
    juce::Array<float> RMS_Larray; // Array of numFrames rms levels
    juce::Array<float> RMS_Rarray; // Array of numFrames rms levels
    const float staticLevelThresh { -60 }; // Level threshold at input (return if not met)
    juce::Array<float> weights; // Array of weights that correspond to delays variable
    float weightedDelaySum { 0 }; // Sum of weighted delay estimates while learning
    float weightsSum { 0 };  // Sum of weights while learning
    const float weightExponent { 0.5 }; // Weight exponent determines how max values are interpreted as weights - a higher value translates to a larger difference in the weight assigned to a small max value vs a large max value
    bool newDelayValueAvailable { false };
    bool whiten { false }; // Divide freq domain xcor by magnitude
};

//==============================================================================
struct UserPreferences
{
    UnitSelect trebleUnits { UnitSelect::freq}, bassBoostUnits { UnitSelect::freq}, bassShelfUnits { UnitSelect::freq}, highPassUnits { UnitSelect::freq};
    bool showRIAA { false }, overSample { true };
    OvrSampAmt ovrSampAmtSetting { OvrSampAmt::medium };
};

//==============================================================================
struct ControlSettings
{
    float gain { -30.f }, bass_shelf_freq{ 10.f }, trebleTurnOverFreq { 2.00f }, delayTime { 0.f };
    float bassTurnOver { 300.f }, highPassFreq { 300.f }, lowPassFreq { 20000.f }, balance { 0.f };
    Slope highPassSlope { Slope::Slope_12 }, lowPassSlope { Slope::Slope_12 };
    bool bassBoostEnable { true }, trebleCutEnable { true }, bassShelfEnable { true };
    bool RIAAEnable { false }, monoSum { false }, vertical { false }, delayEnable { false }, bypassEnable { true }, eqEnable { true};
};

//==============================================================================

const float w_10k = 2*juce::MathConstants<float>::pi*10000.f; // 10kHz converted to rad/sec

const float RIAAT0 = 0.000318f;
const float RIAAT1 = 0.000075f;
const float RIAAT2 = 0.003180f;
const float RIAAT3 = 0.0000016f; // for shelf @ ~100kHz

using Filter = juce::dsp::IIR::Filter<float>;

using Coefficients = juce::dsp::IIR::Coefficients<float>;

std::array<float, 6> getBassCoeffs(ControlSettings settings, double sampleRate);
std::array<float, 6> getTrebleAttenCoeffs(ControlSettings settings, double sampleRate);
std::array<float, 6> getHighPasscoeffs(ControlSettings settings, double sampleRate);
std::array<float, 6> getRIAACoeffs(ControlSettings settings, double sampleRate);

// 1st-order bilinear transform
std::array<float, 6> biLin1(float bs0, float bs1, float as0, float as1, double sampleRate);

// 2nd-order bilinear transform
std::array<float, 6> biLin2(float bs0, float bs1, float bs2, float as0, float as1, float as2, double sampleRate);

inline auto getLowPasscoeffs(const ControlSettings& settings, double sampleRate )
{
    int order;
    if (settings.lowPassSlope == 0){ // filter is bypassed, so order doesn't matter, but order = 0 will throw jassert
        order = 1;
    }
    else {
        order = settings.lowPassSlope;
    }
    return juce::dsp::FilterDesign<float>::designIIRLowpassHighOrderButterworthMethod(settings.lowPassFreq,
                                                                                      sampleRate,
                                                                                      order);
}

inline auto getHighPassCoeffsButter(const ControlSettings& settings, double sampleRate )
{
    int order;
    if (settings.highPassSlope == 0){ // filter is bypassed, so order doesn't matter, but order = 0 will throw jassert
        order = 1;
    }
    else {
        order = settings.highPassSlope;
    }
    return juce::dsp::FilterDesign<float>::designIIRHighpassHighOrderButterworthMethod(settings.highPassFreq,
                                                                                      sampleRate,
                                                                                      order);
}

//==============================================================================
class PlayBackEQAudioProcessor  : public juce::AudioProcessor
                            #if JucePlugin_Enable_ARA
                             , public juce::AudioProcessorARAExtension
                            #endif
{
public:
    //==============================================================================
    PlayBackEQAudioProcessor();
    ~PlayBackEQAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;
    
    static juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout(UserPreferences& uPrefs);
    
    UserPreferences userPreferences;
    
    juce::AudioProcessorValueTreeState apvts {*this, nullptr, "Parameters", createParameterLayout(userPreferences)};
        
    double getUpSampleRate();
    
    size_t getUpSampleExp(); // calculate upsample exponent, upsample rate is sampleRate*2^upsampleExponent
    
    using BlockType = juce::AudioBuffer<float>;
    StereoSummingSampleFifo<BlockType> spectrumAnalyzerSampleFifo{};
    
    juce::AudioChannelSet getOutputFormat() const noexcept { return m_outputFormat; }

    void setPeakLevelInput(int channelIndex, float peakLevel);
    void setPeakLevelOutput(int channelIndex, float peakLevel);

    float getPeakLevelInput(int channelIndex);
    float getPeakLevelOutput(int channelIndex);
        
    Service::PresetManager& getPresetManager() { return *presetManager; };
        
    const UserPreferences& getUserPreferences() { return userPreferences; };
    
    void setUserPreferences(UserPreferences newPrefs) { userPreferences = newPrefs;
                                                        auto userSettings = ap.getUserSettings();
                                                        userSettings->setValue("Treble Units", userPreferences.trebleUnits);
                                                        userSettings->setValue("Bass Boost Units", userPreferences.bassBoostUnits);
                                                        userSettings->setValue("Bass Shelf Units", userPreferences.bassShelfUnits);
                                                        userSettings->setValue("High-Pass Units", userPreferences.highPassUnits);
                                                        userSettings->setValue("Show RIAA in Response", userPreferences.showRIAA);
                                                        userSettings->setValue("Oversampling Amount", userPreferences.ovrSampAmtSetting);
                                                        ap.saveIfNeeded(); };
    
    ControlSettings getSettings(juce::AudioProcessorValueTreeState& apvts);
        
    void prepareOverSamplingProcessors();
    
    void prepareOverSamplers();
    
    void setCurrentSamplesPerBlock(int spb){    currentSamplesPerBlock = spb;  };
    
    // Simply return the property for currentSamplePerBlock, a way to get around the private nature of the property
    int getCurrentSamplesPerBlock(){    return currentSamplesPerBlock;  };
    
    GCCPHAT<float> gccPhat;
    
    void setOverSamplingProcessorsPrepared(const bool status) { overSamplingProcessorsPrepared = status; }
    bool getOverSamplingProcessorsPrepared() { return overSamplingProcessorsPrepared; }
    
    void setTriallMode( bool newTrialModeVal ) {
        trialMode = newTrialModeVal;
    }
    
    bool getTrialMode() {
        return trialMode;
    }
    
    void setUnlocked( bool isUnlocked ) {
        unlocked = isUnlocked;
    }
    
    bool getUnlocked() {
        return unlocked;
    }
    
    void setTrialExpired( bool newTrialExpired ) {
        trialExpired = newTrialExpired;
    }
    
    bool isTrialExpired() {
        return trialExpired;
    }
    
    void setEven( double newEven ) {
        evenMoreObscure = newEven;
    }
    
    auto getEven() {
        return evenMoreObscure;
    }

    void setEmail( const juce::String email){
        auto userSettings = ap.getUserSettings();
        userSettings->setValue("Email Address", email);
        ap.saveIfNeeded();
    }
    
    juce::String getEmail(){
        auto userSettings = ap.getUserSettings();
        return userSettings->getValue("Email Address", "");
    }
    
    void setLicenseString( const juce::String licenseString){
        auto userSettings = ap.getUserSettings();
        userSettings->setValue("License String", licenseString);
        ap.saveIfNeeded();
    }
    
    juce::String getLicenseString(){
        auto userSettings = ap.getUserSettings();
        return userSettings->getValue("License String", "");
    }

private:
    
    const float spectrumAnalyzerFifoBlockTime = 0.04f; // length of audio to collect before pushing it to the spectrum analyzer in seconds
        
    juce::dsp::ProcessorChain<Filter> leftBassChain,
                rightBassChain,
                leftTrebleChain, 
                rightTrebleChain,
                leftRIAAChain,
                rightRIAAChain;
    
    juce::dsp::ProcessorChain<Filter, Filter, Filter> leftHPChain, rightHPChain;
    
    juce::dsp::ProcessorChain<Filter, Filter> leftLPChain, rightLPChain;
    
    juce::dsp::Gain<float> leftGain, rightGain, invert;
    
    juce::dsp::Panner<float> panner;
    
    juce::dsp::DelayLine< float, juce::dsp::DelayLineInterpolationTypes::Thiran > leftDelay, rightDelay;

    monoSum<float> monoSummer;
        
    juce::dsp::Oversampling<float> overSamplerLeft {1, 1, juce::dsp::Oversampling< float >::FilterType::filterHalfBandFIREquiripple, true, false };
    juce::dsp::Oversampling<float> overSamplerRight {1, 1, juce::dsp::Oversampling< float >::FilterType::filterHalfBandFIREquiripple, true, false };
    
    void updateBassFilter(const ControlSettings& settings);
    void updateTrebleAttenFilter(const ControlSettings& settings);
    void updateRIAAFilter(const ControlSettings& settings);

    void updateHighPass(const ControlSettings& settings);
    void updateLowPass(const ControlSettings& settings);

    void updateAllFilters();
    
    juce::AudioChannelSet m_outputFormat;
    std::array<std::atomic<float>, 128> peakLevelsInput;
    std::array<std::atomic<float>, 128> peakLevelsOutput;
    
    std::unique_ptr<Service::PresetManager> presetManager;
    
    juce::ApplicationProperties ap;
    
    int currentSamplesPerBlock;
    
    bool overSamplingProcessorsPrepared {false};
    
    bool trialMode {true}, unlocked {false}, trialExpired {true};
    
    double evenMoreObscure;
    
    MyOfflineUnlockStatus unlockStatusProcessor;

    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PlayBackEQAudioProcessor)
};






