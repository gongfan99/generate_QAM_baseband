#pragma once

#include <cmath>
#include <random>
#include <vector>
#include <algorithm>

void generate_QAM_baseband(){}

const double PI = 3.141592653589793;

unsigned int QAM = 8;

double rollOff = 0.13;
double symbolRate = 168000;
unsigned int totalSymbols = 1680;
unsigned short samplesPerSymbol = 5;
unsigned short rrcTapNum = 128;

double sampleRate = symbolRate * samplesPerSymbol;
double samplePeriod = 1 / sampleRate;

double rrcImpulseResp[rrcTapNum];

int i, j;
double specialRollOffPoint = 1 / (4 * rollOff * symbolRate);
double t;
for (i = 0; i < rrcTapNum; i++){
  t = (i - rrcTapNum / 2) * samplePeriod;
  if (abs(t) < 1e-6){
    rrcImpulseResp[i] = symbolRate * (1 + rollOff * (4 / PI - 1));
  } else if (abs(t - specialRollOffPoint) < 1e-6 || abs(t + specialRollOffPoint) < 1e-6){
    rrcImpulseResp[i] = (symbolRate * rollOff / 1.41421) * ((1 + 2 / PI) * sin(PI / (4 * rollOff)) + (1 - 2 / PI) * cos(PI / (4 * rollOff)));
  } else {
    rrcImpulseResp[i] = (sin(PI * symbolRate * t * (1 - rollOff)) + 4 * roolOff * symbolRate * t * cos(PI * symbolRate * t * (1 + rollOff))) / (PI * t * (1 - 16 * roolOff * roolOff * symbolRate * symbolRate * t * t);
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
std::vector<char> interweaveIQ(4 * totalSamples);
for (i = 0; i < totalSamples; i++){
}