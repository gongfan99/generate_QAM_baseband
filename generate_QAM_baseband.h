/**
    generate_QAM_baseband.h
    Purpose: Generate the baseband IQ signal for QAM; download to the signal generator

    @author Fan Gong
    @version 1.0 07/03/18 
*/

#pragma once

#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include <string>

/**
    Generate random QAM symbols.

    @param Isymbol The container of the returned I symbols.
    @param Qsymbol The container of the returned Q symbols.
    @param totalSymbols The total number of symbols to generate.
    @param QAM The order of QAM e.g. 16 for 16QAM, 64 for 64QAM.
    @return none
*/
void generate_IQ_symbol_QAM_random(std::vector<double>& Isymbol, std::vector<double>& Qsymbol, unsigned int totalSymbols, unsigned int QAM){
  Isymbol.resize(totalSymbols);
  Qsymbol.resize(totalSymbols);
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd(); use gen(std::mt19937::default_seed) for consistent random number for test
  std::uniform_int_distribution<> dis(0, std::sqrt(QAM) - 1);
  int i;
  for (i = 0; i < Isymbol.size(); i++){
    Isymbol[i] = dis(gen) - std::sqrt(QAM) * 0.5 + 0.5;
    Qsymbol[i] = dis(gen) - std::sqrt(QAM) * 0.5 + 0.5;
  }
}

/**
    Generate baseband samples from the I/Q symbols.

    @param Isample The container of the returned I samples.
    @param Qsample The container of the returned Q samples.
    @param Isymbol The input I symbols.
    @param Qsymbol The input Q symbols.
    @param rollOff The roll-off factor of the raised cosine filter.
    @param symbolRate The symbol rate in Hz.
    @param samplesPerSymbol The number of samples per symbol.
    @param rrcTapNum The tap number of the raised cosine filter.
    @return none
*/
void generate_baseband_IQ_waveform(std::vector<double>& Isample, std::vector<double>& Qsample, std::vector<double>& Isymbol, std::vector<double>& Qsymbol, double rollOff, double symbolRate, unsigned int samplesPerSymbol, unsigned int rrcTapNum){
  constexpr double PI = 3.141592653589793;
  constexpr double sqrt2 = 1.41421356237;

  unsigned int totalSymbols = Isymbol.size();
  double sampleRate = symbolRate * samplesPerSymbol;
  double samplePeriod = 1 / sampleRate;

  // generate the impulse response of the raised cosine filter
  std::vector<double> rrcImpulseResp(rrcTapNum, 0);
  int i, j;
  double specialRollOffPoint = 1 / (4 * rollOff);
  double t;
  for (i = 0; i < rrcTapNum; i++){
    t = (i - rrcTapNum * 0.5) * samplePeriod;
    if (abs(t * symbolRate) < 1e-6){
      rrcImpulseResp[i] = symbolRate * (1 + rollOff * (4 / PI - 1));
    } else if (abs(t * symbolRate - specialRollOffPoint) < 1e-6 || abs(t * symbolRate + specialRollOffPoint) < 1e-6){
      rrcImpulseResp[i] = (symbolRate * rollOff / sqrt2) * ((1 + 2 / PI) * sin(PI / (4 * rollOff)) + (1 - 2 / PI) * cos(PI / (4 * rollOff)));
    } else {
      rrcImpulseResp[i] = (sin(PI * symbolRate * t * (1 - rollOff)) + 4 * rollOff * symbolRate * t * cos(PI * symbolRate * t * (1 + rollOff))) / (PI * t * (1 - 16 * rollOff * rollOff * symbolRate * symbolRate * t * t));
    }
  }

  // prepare container of the output I/Q samples
  unsigned int totalSamples = totalSymbols * samplesPerSymbol;
  Isample.resize(totalSamples + rrcTapNum, 0);
  Qsample.resize(totalSamples + rrcTapNum, 0);

  // apply the raised cosine filter
  double *pIsample, *pQsample;
  double Isymbol_i, Qsymbol_i;
  for (i = 0; i < totalSymbols; i++){
    pIsample = Isample.data() + samplesPerSymbol * i;
    pQsample = Qsample.data() + samplesPerSymbol * i;
    Isymbol_i = Isymbol[i];
    Qsymbol_i = Qsymbol[i];
    for (j = 0; j < rrcTapNum; j++){
      // use pIsample instead of Isample to allow vectorization optimization
      pIsample[j] += Isymbol_i * rrcImpulseResp[j];
      pQsample[j] += Qsymbol_i * rrcImpulseResp[j];
    }
  }
  for (j = 0; j < rrcTapNum; j++){
    Isample[j] += Isample[totalSamples + j];
    Qsample[j] += Qsample[totalSamples + j];
  }
  Isample.resize(totalSamples);
  Qsample.resize(totalSamples);
}

/**
    Generate baseband samples for two carriers from the I/Q symbols.

    @param Isample The container of the returned I samples.
    @param Qsample The container of the returned Q samples.
    @param Isymbol_1 The input I symbols for the first carrier.
    @param Qsymbol_1 The input Q symbols for the first carrier.
    @param Isymbol_2 The input I symbols for the second carrier.
    @param Qsymbol_2 The input Q symbols for the second carrier.
    @param rollOff The roll-off factor of the raised cosine filter.
    @param symbolRate The symbol rate in Hz.
    @param samplesPerSymbol The number of samples per symbol.
    @param rrcTapNum The tap number of the raised cosine filter.
    @param freqOffset_1 The frequency offset of the first carrier.
    @param freqOffset_2 The frequency offset of the second carrier.
    @param scale_1 The scale factor of the first carrier.
    @param scale_2 The scale factor of the second carrier.
    @return none
*/
void generate_baseband_IQ_waveform_2carriers(std::vector<double>& Isample, std::vector<double>& Qsample, std::vector<double>& Isymbol_1, std::vector<double>& Qsymbol_1, std::vector<double>& Isymbol_2, std::vector<double>& Qsymbol_2, double rollOff, double symbolRate, unsigned int samplesPerSymbol, unsigned int rrcTapNum, double freqOffset_1, double freqOffset_2, double scale_1, double scale_2){
  constexpr double PI = 3.141592653589793;

  unsigned int totalSymbols = Isymbol_1.size();
  unsigned int totalSamples = totalSymbols * samplesPerSymbol;
  double sampleRate = symbolRate * samplesPerSymbol;
  double samplePeriod = 1 / sampleRate;

  std::vector<double> Isample_1, Qsample_1;
  generate_baseband_IQ_waveform(Isample_1, Qsample_1, Isymbol_1, Qsymbol_1, rollOff, symbolRate, samplesPerSymbol, rrcTapNum);
  std::vector<double> Isample_2, Qsample_2;
  generate_baseband_IQ_waveform(Isample_2, Qsample_2, Isymbol_2, Qsymbol_2, rollOff, symbolRate, samplesPerSymbol, rrcTapNum);

  Isample.resize(totalSamples, 0);
  Qsample.resize(totalSamples, 0);
  double sinWT1, cosWT1, sinWT2, cosWT2;
  for (int i = 0; i < totalSamples; i++){
    sinWT1 = sin(2 * PI * freqOffset_1 * i * samplePeriod);
    cosWT1 = cos(2 * PI * freqOffset_1 * i * samplePeriod);
    sinWT2 = sin(2 * PI * freqOffset_2 * i * samplePeriod);
    cosWT2 = cos(2 * PI * freqOffset_2 * i * samplePeriod);
    Isample[i] = scale_1 * (Isample_1[i] * cosWT1 - Qsample_1[i] * sinWT1) + scale_2 * (Isample_2[i] * cosWT2 - Qsample_2[i] * sinWT2);
    Qsample[i] = scale_1 * (Isample_1[i] * sinWT1 + Qsample_1[i] * cosWT1) + scale_2 * (Isample_2[i] * sinWT2 + Qsample_2[i] * cosWT2);
  }
}

/**
    Download the I/Q waveform to the signal generator.

    @param Isample The I samples.
    @param Qsample The Q samples.
    @param sampleRate The sample rate in Hz
    @param rsrcName The SCPI name of the signal generator e.g. TCPIP0::192.168.0.29::inst0::INSTR or GPIB0::19::INSTR
    @param IdcOffset The I channel DC offset of the signal generator.
    @param QdcOffset The Q channel DC offset of the signal generator.
    @param IQphaseImbalance The I/Q phase imbalance of the signal generator.
    @param IQamplitudeImbalance The I/Q amplitude imbalance of the signal generator.
    @return none
*/
void download_IQ_sample_VSG(std::vector<double>& Isample, std::vector<double>& Qsample, double sampleRate, std::string& rsrcName, double IdcOffset, double QdcOffset, double IQphaseImbalance, double IQamplitudeImbalance){
  auto minmaxIsample = std::minmax_element(Isample.begin(), Isample.end());
  auto minmaxQsample = std::minmax_element(Qsample.begin(), Qsample.end());
  double temp1[4] = {abs(*minmaxIsample.first), *minmaxIsample.second, abs(*minmaxQsample.first), *minmaxQsample.second};
  double maxIQsample = *std::max_element(temp1, temp1 + 4);

  // generate the binary block data
  unsigned int totalSamples = Isample.size();
  std::string interleaveIQ(4 * totalSamples, NULL);
  long int16Value;
  double Inormalized, Qnormalized;
  for (int i = 0; i < totalSamples; i++){
    Inormalized = Isample[i] * 23000 / maxIQsample;
    Qnormalized = Qsample[i] * 23000 / maxIQsample;
    Inormalized -= IdcOffset;
    Qnormalized -= QdcOffset;
    Inormalized -= Qnormalized * IQphaseImbalance * 0.5;
    Qnormalized -= Inormalized * IQphaseImbalance * 0.5;
    Inormalized *= 1 + IQamplitudeImbalance * 0.5;
    Qnormalized *= 1 - IQamplitudeImbalance * 0.5;
    int16Value = std::lround(Inormalized);
    // intel uC is little-endian while the signal generator is big-endian
    interleaveIQ[4 * i] = (int16Value >> 8) & 0xFF;
    interleaveIQ[4 * i + 1] = int16Value & 0xFF;
    int16Value = std::lround(Qnormalized);
    interleaveIQ[4 * i + 2] = (int16Value >> 8) & 0xFF;
    interleaveIQ[4 * i + 3] = int16Value & 0xFF;
  }

  // download the binary block data to the signal generator
  ViStatus status;
  ViUInt32 retCount;
  ViSession defaultRM, instr;
  status = viOpenDefaultRM(&defaultRM);
  status = viOpen(defaultRM, (ViRsrc)rsrcName.data(), VI_NULL, VI_NULL, &instr);
  
  std::string cmd;
  cmd = ":SOURce:RADio:ARB:STATe OFF\n";
  status = viWrite(instr, (ViBuf)cmd.data(), (ViUInt32)cmd.size(), (ViPUInt32)&retCount);
  
  cmd = std::to_string(interleaveIQ.size());
  cmd = ":MEMory:DATA \"WFM1:FILE1\",#" + std::to_string(cmd.size()) + cmd + interleaveIQ + "\n";
  status = viWrite(instr, (ViBuf)cmd.data(), (ViUInt32)cmd.size(), (ViPUInt32)&retCount);
  
  cmd = ":SOURce:RADio:ARB:WAVeform \"WFM1:FILE1\"\n";
  status = viWrite(instr, (ViBuf)cmd.data(), (ViUInt32)cmd.size(), (ViPUInt32)&retCount);
  
  cmd = ":SOURce:RADio:ARB:SCLock:RATE " + std::to_string(sampleRate) + "\n";
  status = viWrite(instr, (ViBuf)cmd.data(), (ViUInt32)cmd.size(), (ViPUInt32)&retCount);
  
  cmd = ":SOURce:RADio:ARB:STATe ON\n";
  status = viWrite(instr, (ViBuf)cmd.data(), (ViUInt32)cmd.size(), (ViPUInt32)&retCount);
  
  cmd = ":SYSTem:COMMunicate:GTLocal\n";
  status = viWrite(instr, (ViBuf)cmd.data(), (ViUInt32)cmd.size(), (ViPUInt32)&retCount);
  
  viClose(instr);
  viClose(defaultRM);
}