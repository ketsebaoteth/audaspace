#pragma once

#include "fx/Effect.h"

AUD_NAMESPACE_BEGIN

/*
 * This effect compresses the audio signal. based on the specified threshold and ratio and other parameters.
 */

class AUD_API Compressor : public Effect
{
private:
    float m_threshold;
    float m_ratio;
    float m_attack;
    float m_release;
    float m_makeupGain;
    float m_kneeWidth;
    float m_lookaheadMs;

    //delete copy constructor and assignment operator
    Compressor(const Compressor&) = delete;
    Compressor& operator=(const Compressor&) = delete;
public:
    /**
    * @brief Construct a new Compressor object
    * 
    * @param threshold Threshold level in dBFS (e.g., -18.0)
    * @param ratio Compression ratio (e.g., 4.0 means 4:1)
    * @param attack Attack time in milliseconds
    * @param release Release time in milliseconds
    * @param gain Output gain in dB
    * @param kneeWidth Knee width in dB
    * @param lookaheadMs Lookahead time in milliseconds
    */
    Compressor(std::shared_ptr<ISound> sound,float threshold, float ratio, float attack, float release, float makeupGain, float kneeWidth,float lookaheadMs);
    /**
     * @brief Get the threshold level
     * @return Threshold level in dBFS
     */
    float getThreshold() const;
    /**
     * @brief Get the compression ratio
     * @return Compression ratio
     */
    float getRatio() const;
    /**
     * @brief Get the attack time
     * @return Attack time in milliseconds
     */
    float getAttack() const;
    /**
     * @brief Get the release time
     * @return Release time in milliseconds
     */
    float getRelease() const;

    /**
     * @brief Get the output gain
     * @return Output gain in dB
     */
    float getGain() const;

    virtual std::shared_ptr<IReader> createReader();
};

AUD_NAMESPACE_END
