#include <iostream>
#include <string>
#include <cmath>
#include <random>

#include "visa.h"
#include "generate_QAM_baseband.h"

int main(int argc, char* argv[]) {
  unsigned int QAM = 16; // 64QAM
  double rollOff = 0.18;
  double symbolRate = 32000;
  unsigned int totalSymbols = 3200;
  double freqOffset_1 = -0.3e6;
  double freqOffset_2 = 0.3e6;
  double scale_1 = 1;
  double scale_2 = 1;
  bool singleCarrier = true;
  std::string rsrcName = "TCPIP0::192.168.0.30::inst0::INSTR"; // or "GPIB::1::INSTR"

  unsigned int samplesPerSymbol = 4;
  unsigned int rrcTapNum = 128 * samplesPerSymbol;

  if (argc < 3){
    std::cout << "Usage: test totalSymbols samplesPerSymbol" << std::endl;
  } else {
    totalSymbols = atoi(argv[1]);
    samplesPerSymbol = atoi(argv[2]);
  }

  std::vector<double> Isample, Qsample;
  if (singleCarrier){
    std::vector<double> Isymbol, Qsymbol;
    generate_IQ_symbol_QAM_random(Isymbol, Qsymbol, totalSymbols, QAM);

    generate_baseband_IQ_waveform(Isample, Qsample, Isymbol, Qsymbol, rollOff, symbolRate, samplesPerSymbol, rrcTapNum);
  } else {
    samplesPerSymbol = 4 * (abs(freqOffset_1) > abs(freqOffset_2)? abs(freqOffset_1) : abs(freqOffset_2)) / symbolRate;
    rrcTapNum = 128 * samplesPerSymbol;

    std::vector<double> Isymbol_1, Qsymbol_1;
    generate_IQ_symbol_QAM_random(Isymbol_1, Qsymbol_1, totalSymbols, QAM);
    std::vector<double> Isymbol_2, Qsymbol_2;
    generate_IQ_symbol_QAM_random(Isymbol_2, Qsymbol_2, totalSymbols, QAM);

    generate_baseband_IQ_waveform_2carriers(Isample, Qsample, Isymbol_1, Qsymbol_1, Isymbol_2, Qsymbol_2, rollOff, symbolRate, samplesPerSymbol, rrcTapNum, freqOffset_1, freqOffset_2, scale_1, scale_2);
  }
  
  double IdcOffset = 0, QdcOffset = 0, IQphaseImbalance = 0, IQamplitudeImbalance = 0;
  double sampleRate = symbolRate * samplesPerSymbol;
  download_IQ_sample_VSG(Isample, Qsample, sampleRate, rsrcName, IdcOffset, QdcOffset, IQphaseImbalance, IQamplitudeImbalance);

  std::cout << "generate_QAM_baseband" << std::endl;
}