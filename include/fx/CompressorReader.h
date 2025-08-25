#pragma once

#include "fx/EffectReader.h"
#include <vector>

AUD_NAMESPACE_BEGIN

class CompressorReader : public EffectReader
{
public:
    CompressorReader(std::shared_ptr<IReader> reader, float threshold, float ratio, float attack, float release, float makeupGain, float kneeWidth, float lookaheadMs);
    // Override read method to apply compression
    virtual void read(int& length, bool& eos, sample_t* buffer) override;

private:
    float m_threshold;
    float m_ratio;
    float m_attack;
    float m_release;
    float m_makeupGain;
    float m_kneeWidth;
    int m_lookaheadSamples;
    std::vector<sample_t> m_delayBuffer;
    int m_delayBufferWritePos;


    float m_attackCoeff;
    float m_releaseCoeff;
    std::vector<float> m_rmsState;

    int m_channels;
    int m_windowSize;
    std::vector<float> m_envelope;

    CompressorReader(const CompressorReader&) = delete;
    CompressorReader& operator=(const CompressorReader&) = delete;
};

AUD_NAMESPACE_END
