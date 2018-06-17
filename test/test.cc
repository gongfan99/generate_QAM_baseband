#include <iostream>
#include <string>
#include <cmath>
#include <random>

#include "visa.h"
#include "generate_QAM_baseband.h"

int main(int argc, char* argv[]) {
  unsigned int QAM = 64;
  double rollOff = 0.13;
  double symbolRate = 168000;
  unsigned int totalSymbols = 16;
  unsigned int samplesPerSymbol = 2;
  unsigned int rrcTapNum = 32;
  double freqGap = 1e6;
  bool singleCarrier = true;
  std::string rsrcName = "GPIB::1::INSTR";

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
    std::vector<double> Isymbol_1, Qsymbol_1;
    generate_IQ_symbol_QAM_random(Isymbol_1, Qsymbol_1, totalSymbols, QAM);
    std::vector<double> Isymbol_2, Qsymbol_2;
    generate_IQ_symbol_QAM_random(Isymbol_2, Qsymbol_2, totalSymbols, QAM);

    generate_baseband_IQ_waveform_2carriers(Isample, Qsample, Isymbol_1, Qsymbol_1, Isymbol_2, Qsymbol_2, rollOff, symbolRate, samplesPerSymbol, rrcTapNum, freqGap);
  }
  
  double sampleRate = symbolRate * samplesPerSymbol;
  download_IQ_sample_VSG(Isample, Qsample, sampleRate, rsrcName);

  std::cout << "generate_QAM_baseband" << std::endl;
}