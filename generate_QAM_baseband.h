#pragma once

#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include <string>

std::string generate_QAM_baseband(){

const double PI = 3.141592653589793;

unsigned int QAM = 8;

double rollOff = 0.13;
double symbolRate = 168000;
unsigned int totalSymbols = 1680;
unsigned short samplesPerSymbol = 5;
constexpr unsigned short rrcTapNum = 32;

double sampleRate = symbolRate * samplesPerSymbol;
double samplePeriod = 1 / sampleRate;

double rrcImpulseResp[rrcTapNum];

int i, j;
double specialRollOffPoint = 1 / (4 * rollOff);
double t;
for (i = 0; i < rrcTapNum; i++){
  t = (i - rrcTapNum / 2) * samplePeriod;
  if (abs(t * symbolRate) < 1e-6){
    rrcImpulseResp[i] = symbolRate * (1 + rollOff * (4 / PI - 1));
  } else if (abs(t * symbolRate - specialRollOffPoint) < 1e-6 || abs(t * symbolRate + specialRollOffPoint) < 1e-6){
    rrcImpulseResp[i] = (symbolRate * rollOff / 1.41421356237) * ((1 + 2 / PI) * sin(PI / (4 * rollOff)) + (1 - 2 / PI) * cos(PI / (4 * rollOff)));
  } else {
    rrcImpulseResp[i] = (sin(PI * symbolRate * t * (1 - rollOff)) + 4 * rollOff * symbolRate * t * cos(PI * symbolRate * t * (1 + rollOff))) / (PI * t * (1 - 16 * rollOff * rollOff * symbolRate * symbolRate * t * t));
  }
}

unsigned int totalSamples = totalSymbols * samplesPerSymbol;
std::vector<double> Isample(totalSamples, 0);
std::vector<double> Qsample(totalSamples, 0);

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_int_distribution<> dis(-QAM, QAM);
int symbolIValue, symbolQValue;
for (i = 0; i < totalSymbols; i++){
  symbolIValue = dis(gen);
  symbolQValue = dis(gen);
  for (j = 0; j < rrcTapNum; j++){
    Isample[(samplesPerSymbol * i + j) % totalSamples] += symbolIValue * rrcImpulseResp[j];
    Qsample[(samplesPerSymbol * i + j) % totalSamples] += symbolQValue * rrcImpulseResp[j];
  }
}

auto minmaxIsample = std::minmax_element(Isample.begin(), Isample.end());
auto minmaxQsample = std::minmax_element(Qsample.begin(), Qsample.end());
double temp1[4] = {abs(*minmaxIsample.first), *minmaxIsample.second, abs(*minmaxQsample.first), *minmaxQsample.second};
double maxIQsample = *std::max_element(temp1, temp1 + 4);
std::string interleaveIQ(4 * totalSamples, NULL);
short int int16Value;
for (i = 0; i < totalSamples; i++){
  int16Value = Isample[i] * 23000 / maxIQsample;
  interleaveIQ[4 * i] = (int16Value >> 8) & 0xFF;
  interleaveIQ[4 * i + 1] = int16Value & 0xFF;
  int16Value = Qsample[i] * 23000 / maxIQsample;
  interleaveIQ[4 * i + 2] = (int16Value >> 8) & 0xFF;
  interleaveIQ[4 * i + 3] = int16Value & 0xFF;
}
for (double temp2 : rrcImpulseResp){
std::cout << temp2 << std::endl;
}

return interleaveIQ;
}