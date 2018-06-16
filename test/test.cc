#include <iostream>
#include <string>

#include "visa.h"
#include "generate_QAM_baseband.h"

int main() {
  std::string interleaveIQ(generate_QAM_baseband());
  std::string cmd = std::to_string(interleaveIQ.size());
  cmd = ":MEM:DATA \"WFM1:FILE1\", #%d%d" + std::to_string(cmd.size()) + cmd;
  ViStatus status;
  ViUInt32 retCount;
  ViSession vi;
  status = viWrite(vi, (ViBuf)cmd.data(), (ViUInt32)cmd.size(), (ViPUInt32)&retCount);
  status = viWrite(vi, (ViBuf)interleaveIQ.data(), (ViUInt32)interleaveIQ.size(), (ViPUInt32)&retCount);
  status = viWrite(vi, (ViBuf)std::string("\n").data(), 1, (ViPUInt32)&retCount);
  std::cout << "generate_QAM_baseband" << std::endl;
}