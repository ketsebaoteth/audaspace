#include "fx/CompressorReader.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>


AUD_NAMESPACE_BEGIN

CompressorReader::CompressorReader(std::shared_ptr<IReader> reader, float threshold, float ratio, float attack, float release, float makeupGain, float kneeWidth, float lookaheadMs) :
    EffectReader(reader),
    m_threshold(threshold),
    m_ratio(ratio),
    m_attack(attack),
    m_release(release),
    m_makeupGain(makeupGain),
    m_kneeWidth(kneeWidth)
{
    Specs specs = m_reader->getSpecs();
    m_channels = specs.channels;
    m_windowSize = std::max(static_cast<int>(0.02f * specs.rate), 1); // 20ms window
    m_lookaheadSamples = 0;
    m_delayBufferWritePos = 0;

    // Convert attack/release from ms to seconds for coefficient calculation
    float attack_sec = std::max(m_attack * 0.001f, 1e-6f);
    float release_sec = std::max(m_release * 0.001f, 1e-6f);
    m_attackCoeff = std::exp(-1.0f / (attack_sec * specs.rate));
    m_releaseCoeff = std::exp(-1.0f / (release_sec * specs.rate));

    m_envelope.resize(m_channels, 0.0f);
    m_rmsState.resize(m_channels, 0.0f);

    // Lookahead buffer initialization
    if (lookaheadMs > 0.0f && m_channels > 0) {
        m_lookaheadSamples = static_cast<int>(lookaheadMs * 0.001f * specs.rate);
        if (m_lookaheadSamples > 0) {
            m_delayBuffer.resize(m_lookaheadSamples * m_channels, 0.0f);
        }
    }
}

void CompressorReader::read(int& length, bool& eos, sample_t* buffer)
{
    if (m_lookaheadSamples > 0) {
        // With lookahead, we read into a temporary "sidechain" buffer for analysis
        std::vector<sample_t> sidechainBuffer(length * m_channels);
        m_reader->read(length, eos, sidechainBuffer.data());

        if (length == 0) return;

        Specs specs = m_reader->getSpecs();
        const int channels = specs.channels;
        const int total_samples = length * channels;

        float threshold = m_threshold;
        float ratio = m_ratio;
        float makeup = pow(10.0f, m_makeupGain / 20.0f);
        float knee = m_kneeWidth;

        if (m_rmsState.size() != channels) m_rmsState.assign(channels, 0.0f);
        if (m_envelope.size() != channels) m_envelope.assign(channels, 0.0f);

        float rms_coeff = std::exp(-1.0f / (0.02f * specs.rate));
        float min_rms = 1e-12f;

        float attack_coeff = m_attackCoeff;
        float release_coeff = m_releaseCoeff;

        for (int i = 0; i < total_samples; ++i) {
            int ch = i % channels;
            float sample = sidechainBuffer[i]; // Use the unprocessed "future" signal for detection

            float detector = 0.0f;
            
            float sq = sample * sample;
            m_rmsState[ch] = rms_coeff * m_rmsState[ch] + (1.0f - rms_coeff) * sq;
            detector = std::sqrt(std::max(m_rmsState[ch], min_rms));

            float detector_db = 20.0f * log10(detector + min_rms);

            // Soft knee gain reduction calculation 
            float kneestart = threshold - knee / 2.0f;
            float kneeend = threshold + knee / 2.0f;
            float gainReduction_db = 0.0f;

            if (detector_db < kneestart)
                gainReduction_db = 0.0f;
            else if (detector_db > kneeend)
                gainReduction_db = (detector_db - threshold) * (1.0f - 1.0f / ratio);
            else {
                float x = detector_db - kneestart;
                float y = x * x / (knee * 2.0f); // Using quadratic interpolation inside the knee
                gainReduction_db = y * (1.0f / ratio - 1.0f);
            }

            // Envelope follower for smooth gain reduction
            float env = m_envelope[ch];
            if (gainReduction_db > env)
                env = attack_coeff * env + (1.0f - attack_coeff) * gainReduction_db;
            else
                env = release_coeff * env + (1.0f - release_coeff) * gainReduction_db;
            m_envelope[ch] = env;

            // Get delayed sample and apply gain (Your flicker-free logic)
            int delayIndex = (m_delayBufferWritePos + i - m_lookaheadSamples * channels + m_delayBuffer.size()) % m_delayBuffer.size();
            float delayedSample = m_delayBuffer[delayIndex];
            buffer[i] = delayedSample; // Output the delayed sample before processing

            m_delayBuffer[(m_delayBufferWritePos + i) % m_delayBuffer.size()] = sidechainBuffer[i]; // Store current sample in delay buffer

            float smoothedGainReduction = -env;
            float gain = pow(10.0f, smoothedGainReduction / 20.0f) * makeup;
            buffer[i] *= gain; // Apply gain to the delayed sample

            // Gentle soft knee limiter at ±1.0 (on the output)
            const float limit_knee = 0.1f;
            float out = buffer[i];
            if (out > 1.0f - limit_knee) {
                if (out < 1.0f + limit_knee)
                    buffer[i] = 1.0f - (out - 1.0f - limit_knee) * (out - 1.0f - limit_knee) / (4.0f * limit_knee);
                else
                    buffer[i] = 1.0f;
            } else if (out < -1.0f + limit_knee) {
                if (out > -1.0f - limit_knee)
                    buffer[i] = -1.0f + (out + 1.0f + limit_knee) * (out + 1.0f + limit_knee) / (4.0f * limit_knee);
                else
                    buffer[i] = -1.0f;
            }
        }
        m_delayBufferWritePos = (m_delayBufferWritePos + total_samples) % m_delayBuffer.size();

    } else {
        // no lookahead path
        m_reader->read(length, eos, buffer);
        if (length == 0) return;

        Specs specs = m_reader->getSpecs();
        const int channels = specs.channels;
        const int total_samples = length * channels;

        float threshold = m_threshold;
        float ratio = m_ratio;
        float makeup = pow(10.0f, m_makeupGain / 20.0f);
        float knee = m_kneeWidth;

        if (m_rmsState.size() != channels) m_rmsState.assign(channels, 0.0f);
        if (m_envelope.size() != channels) m_envelope.assign(channels, 0.0f);

        float rms_coeff = std::exp(-1.0f / (0.02f * specs.rate));
        float min_rms = 1e-12f;

        float attack_coeff = m_attackCoeff;
        float release_coeff = m_releaseCoeff;

        for (int i = 0; i < total_samples; ++i) {
            int ch = i % channels;
            float sample = buffer[i];

            float detector = 0.0f;

            float sq = sample * sample;
            m_rmsState[ch] = rms_coeff * m_rmsState[ch] + (1.0f - rms_coeff) * sq;
            detector = std::sqrt(std::max(m_rmsState[ch], min_rms));

            float detector_db = 20.0f * log10(detector + min_rms);

            // Soft knee gain reduction calculation
            float kneestart = threshold - knee / 2.0f;
            float kneeend = threshold + knee / 2.0f;
            float gainReduction_db = 0.0f;

            if (detector_db < kneestart)
                gainReduction_db = 0.0f;
            else if (detector_db > kneeend)
                gainReduction_db = (detector_db - threshold) * (1.0f - 1.0f / ratio);
            else {
                float x = detector_db - kneestart;
                float y = x * x / (knee * 2.0f);
                gainReduction_db = y * (1.0f / ratio - 1.0f);
            }

            //  Envelope follower for smooth gain reduction
            float env = m_envelope[ch];
            if (gainReduction_db > env)
                env = attack_coeff * env + (1.0f - attack_coeff) * gainReduction_db;
            else
                env = release_coeff * env + (1.0f - release_coeff) * gainReduction_db;
            m_envelope[ch] = env;

            // Apply gain reduction and makeup gain
            float smoothedGainReduction = -env;
            float gain = pow(10.0f, smoothedGainReduction / 20.0f) * makeup;
            float out = sample * gain;

            // Gentle soft knee limiter at ±1.0
            const float limit_knee = 0.1f;
            if (out > 1.0f - limit_knee) {
                if (out < 1.0f + limit_knee)
                    buffer[i] = 1.0f - (out - 1.0f - limit_knee) * (out - 1.0f - limit_knee) / (4.0f * limit_knee);
                else
                    buffer[i] = 1.0f;
            } else if (out < -1.0f + limit_knee) {
                if (out > -1.0f - limit_knee)
                    buffer[i] = -1.0f + (out + 1.0f + limit_knee) * (out + 1.0f + limit_knee) / (4.0f * limit_knee);
                else
                    buffer[i] = -1.0f;
            } else {
                buffer[i] = out;
            }
        }
    }
}

AUD_NAMESPACE_END