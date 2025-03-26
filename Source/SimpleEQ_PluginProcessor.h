
#pragma once

#include <juce_dsp/juce_dsp.h>

// template <typename T>
// struct Fifo
// {
//     void prepare (int numChannels, int numSamples)
//     {
//         static_assert (std::is_same_v<T, juce::AudioBuffer<float>>,
//             "prepare(numChannels, numSamples) should only be used when the Fifo is holding juce::AudioBuffer<float>");
//         for (auto& buffer : buffers)
//         {
//             buffer.setSize (numChannels,
//                 numSamples,
//                 false, //clear everything?
//                 true, //including the extra space?
//                 true); //avoid reallocating if you can?
//             buffer.clear();
//         }
//     }
//
//     void prepare (size_t numElements)
//     {
//         static_assert (std::is_same_v<T, std::vector<float>>,
//             "prepare(numElements) should only be used when the Fifo is holding std::vector<float>");
//         for (auto& buffer : buffers)
//         {
//             buffer.clear();
//             buffer.resize (numElements, 0);
//         }
//     }
//
//     bool push (const T& t)
//     {
//         auto write = fifo.write (1);
//         if (write.blockSize1 > 0)
//         {
//             buffers[write.startIndex1] = t;
//             return true;
//         }
//
//         return false;
//     }
//
//     bool pull (T& t)
//     {
//         auto read = fifo.read (1);
//         if (read.blockSize1 > 0)
//         {
//             t = buffers[read.startIndex1];
//             return true;
//         }
//
//         return false;
//     }
//
//     int getNumAvailableForReading() const
//     {
//         return fifo.getNumReady();
//     }
//
// private:
//     static constexpr int Capacity = 30;
//     std::array<T, Capacity> buffers;
//     juce::AbstractFifo fifo { Capacity };
// };

// enum Channel {
//     Left, //effectively 0
//     Right //effectively 1
// };
//
// template <typename BlockType>
// struct StereoSummingSampleFifo
// {
//     StereoSummingSampleFifo()
//     {
//         prepared.set (false);
//     }
//
//     void update (const BlockType& buffer)
//     {
//         jassert (prepared.get());
//
//         if (buffer.getNumChannels() == 1) // mono instance
//         {
//             auto* leftChannelPtr = buffer.getReadPointer (0);
//             for (int i = 0; i < buffer.getNumSamples(); ++i)
//             {
//                 const float samp = leftChannelPtr[i];
//                 pushNextSampleIntoFifo (samp);
//             }
//         }
//         else // stereo or ? instance
//         {
//             auto* leftChannelPtr = buffer.getReadPointer (0);
//             auto* rightChannelPtr = buffer.getReadPointer (1);
//             // sum L and R
//             for (int i = 0; i < buffer.getNumSamples(); ++i)
//             {
//                 const float samp = leftChannelPtr[i] + rightChannelPtr[i];
//                 pushNextSampleIntoFifo (samp);
//             }
//         }
//     }
//
//     void prepare (int bufferSize)
//     {
//         prepared.set (false);
//         size.set (bufferSize);
//
//         bufferToFill.setSize (1, //channel
//             bufferSize, //num samples
//             false, //keepExistingContent
//             true, //clear extra space
//             true); //avoid reallocating
//         audioBufferFifo.prepare (1, bufferSize);
//         fifoIndex = 0;
//         prepared.set (true);
//     }
//     //==============================================================================
//     int getNumCompleteBuffersAvailable() const { return audioBufferFifo.getNumAvailableForReading(); }
//     bool isPrepared() const { return prepared.get(); }
//     int getSize() const { return size.get(); }
//     //==============================================================================
//     bool getAudioBuffer (BlockType& buf) { return audioBufferFifo.pull (buf); }
//
// private:
//     int fifoIndex = 0;
//     Fifo<BlockType> audioBufferFifo;
//     BlockType bufferToFill;
//     juce::Atomic<bool> prepared = false;
//     juce::Atomic<int> size = 0;
//
//     void pushNextSampleIntoFifo (float sample)
//     {
//         if (fifoIndex == bufferToFill.getNumSamples())
//         {
//             auto ok = audioBufferFifo.push (bufferToFill);
//
//             juce::ignoreUnused (ok);
//
//             fifoIndex = 0;
//         }
//
//         bufferToFill.setSample (0, fifoIndex, sample);
//         ++fifoIndex;
//     }
// };

using coefficientsPointer = juce::dsp::IIR::Filter<float>::CoefficientsPtr;

inline void updateCoefficients (const coefficientsPointer& old, const coefficientsPointer& replacements)
{
    *old = *replacements;
};
