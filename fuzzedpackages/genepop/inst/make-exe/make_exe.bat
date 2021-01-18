:: not *.cpp because of RcppExports.cpp
g++ -std=c++11 -O3 -o genepop.exe bootstrap.cpp   conversions.cpp CT_tests.cpp F_est.cpp genepop.cpp GenepopS.cpp HW_tests.cpp multimig.cpp proba.cpp settings.cpp tools.cpp
