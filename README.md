# generate_QAM_baseband
Generate the baseband IQ signal for QAM with C++. It has been verified on Windows 10.

# build example
```c
// Windows 10, VS 2015 installed
git clone https://github.com/gongfan99/generate_QAM_baseband.git
cd test
build
test
```

# usage
More usage cases can be found in test/test.cc

It can generate single carrier or two carriers QAM baseband signal. 

# Reference
https://en.wikipedia.org/wiki/Root-raised-cosine_filter
https://en.wikipedia.org/wiki/Quadrature_amplitude_modulation