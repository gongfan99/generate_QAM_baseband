#pragma once

#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include <string>

void generate_IQ_symbol_QAM_random(std::vector<double>& Isymbol, std::vector<double>& Qsymbol, unsigned int totalSymbols, unsigned int QAM){
  Isymbol.resize(totalSymbols);
  Qsymbol.resize(totalSymbols);
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> dis(0, std::sqrt(QAM) - 1);
  int i;
  for (i = 0; i < Isymbol.size(); i++){
    Isymbol[i] = dis(gen) - std::sqrt(QAM) * 0.5 + 0.5;
    Qsymbol[i] = dis(gen) - std::sqrt(QAM) * 0.5 + 0.5;
  }
}

void generate_baseband_IQ_waveform(std::vector<double>& Isample, std::vector<double>& Qsample, std::vector<double>& Isymbol, std::vector<double>& Qsymbol, double rollOff, double symbolRate, unsigned int samplesPerSymbol, unsigned int rrcTapNum){
  constexpr double PI = 3.141592653589793;
  constexpr double sqrt2 = 1.41421356237;

  unsigned int totalSymbols = Isymbol.size();
  double sampleRate = symbolRate * samplesPerSymbol;
  double samplePeriod = 1 / sampleRate;

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

  unsigned int totalSamples = totalSymbols * samplesPerSymbol;
  Isample.resize(totalSamples + rrcTapNum, 0);
  Qsample.resize(totalSamples + rrcTapNum, 0);

  for (i = 0; i < totalSymbols; i++){
    for (j = 0; j < rrcTapNum; j++){
      Isample[samplesPerSymbol * i + j] += Isymbol[i] * rrcImpulseResp[j];
      Qsample[samplesPerSymbol * i + j] += Qsymbol[i] * rrcImpulseResp[j];
    }
  }
  for (j = 0; j < rrcTapNum; j++){
    Isample[j] += Isample[totalSamples + i];
    Qsample[j] += Qsample[totalSamples + i];
  }
  Isample.resize(totalSamples);
  Qsample.resize(totalSamples);
}

void generate_baseband_IQ_waveform_2carriers(std::vector<double>& Isample, std::vector<double>& Qsample, std::vector<double>& Isymbol_1, std::vector<double>& Qsymbol_1, std::vector<double>& Isymbol_2, std::vector<double>& Qsymbol_2, double rollOff, double symbolRate, unsigned int samplesPerSymbol, unsigned int rrcTapNum, double freqOffset_1, double freqOffset_2){
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
    Isample[i] = Isample_1[i] * cosWT1 - Qsample_1[i] * sinWT1 + Isample_2[i] * cosWT2 - Qsample_2[i] * sinWT2;
    Qsample[i] = Isample_1[i] * sinWT1 + Qsample_1[i] * cosWT1 + Isample_2[i] * sinWT2 + Qsample_2[i] * cosWT2;
  }
}

void download_IQ_sample_VSG(std::vector<double>& Isample, std::vector<double>& Qsample, double sampleRate, std::string& rsrcName){
  auto minmaxIsample = std::minmax_element(Isample.begin(), Isample.end());
  auto minmaxQsample = std::minmax_element(Qsample.begin(), Qsample.end());
  double temp1[4] = {abs(*minmaxIsample.first), *minmaxIsample.second, abs(*minmaxQsample.first), *minmaxQsample.second};
  double maxIQsample = *std::max_element(temp1, temp1 + 4);

  unsigned int totalSamples = Isample.size();
  std::string interleaveIQ(4 * totalSamples, NULL);
  long int16Value;
  for (int i = 0; i < totalSamples; i++){
    int16Value = std::lround(Isample[i] * 23000 / maxIQsample);
    interleaveIQ[4 * i] = (int16Value >> 8) & 0xFF;
    interleaveIQ[4 * i + 1] = int16Value & 0xFF;
    int16Value = std::lround(Qsample[i] * 23000 / maxIQsample);
    interleaveIQ[4 * i + 2] = (int16Value >> 8) & 0xFF;
    interleaveIQ[4 * i + 3] = int16Value & 0xFF;
  }

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