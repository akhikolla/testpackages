library(Rmixmod)
data(heterodatatest)
data(heterodatatrain)
learn<-mixmodLearn(heterodatatrain[-1], knownLabels=heterodatatrain$V1)
# predict<-mixmodPredict(heterodatatest[-1], classificationRule=learn["bestResult"])
# missclassified<-sum(as.integer(predict@partition)-as.integer(heterodatatest$V1))

# TODO
# ERROR on CRAN on clang-UBSAN, gcc-UBSAN: Tests of memory access errors using Undefined Behavior Sanitizer
# https://www.stats.ox.ac.uk/pub/bdr/memtests/README.txt
#
# > predict<-mixmodPredict(heterodatatest[-1], classificationRule=learn["bestResult"])
# Kernel/Parameter/GaussianEDDAParameter.cpp:297:95: runtime error: downcast of address 0x60e00000c0c0 which does not point to an object of type 'GaussianGeneralParameter'
# 0x60e00000c0c0: note: object is of type 'XEM::GaussianDiagParameter'
# 02 00 80 36  20 0b e0 b0 c7 7f 00 00  02 00 00 00 00 00 00 00  02 00 00 00 00 00 00 00  30 e8 06 00
# ^~~~~~~~~~~~~~~~~~~~~~~
#   vptr for 'XEM::GaussianDiagParameter'
# #0 0x7fc7b0333d29 in XEM::GaussianEDDAParameter::initUSER(XEM::Parameter*) Kernel/Parameter/GaussianEDDAParameter.cpp:297
# #1 0x7fc7b03d1abf in XEM::GaussianDiagParameter::initUSER(XEM::Parameter*) Kernel/Parameter/GaussianDiagParameter.cpp:172
# #2 0x7fc7b04b6e16 in XEM::Model::initUSER(XEM::Parameter*) Kernel/Model/Model.cpp:998
# #3 0x7fc7b030a21d in XEM::PredictStrategy::run(XEM::Model*) DiscriminantAnalysis/Predict/PredictStrategy.cpp:59
# #4 0x7fc7b0308fad in XEM::PredictMain::run(XEM::IoMode, int, int) DiscriminantAnalysis/Predict/PredictMain.cpp:141
# #5 0x7fc7b0190024 in predictMain /data/gannet/ripley/R/packages/tests-gcc-SAN/Rmixmod/src/predictMain.cpp:295
# #6 0x57b493 in R_doDotCall /data/gannet/ripley/R/svn/R-devel/src/main/dotcode.c:598
# #7 0x5839ac in do_dotcall /data/gannet/ripley/R/svn/R-devel/src/main/dotcode.c:1281
# #8 0x61e951 in bcEval /data/gannet/ripley/R/svn/R-devel/src/main/eval.c:7105
# #9 0x66b6ff in Rf_eval /data/gannet/ripley/R/svn/R-devel/src/main/eval.c:723
# #10 0x670b85 in R_execClosure /data/gannet/ripley/R/svn/R-devel/src/main/eval.c:1888
# #11 0x6732a4 in Rf_applyClosure /data/gannet/ripley/R/svn/R-devel/src/main/eval.c:1814
# #12 0x66bb48 in Rf_eval /data/gannet/ripley/R/svn/R-devel/src/main/eval.c:846
# #13 0x679581 in do_set /data/gannet/ripley/R/svn/R-devel/src/main/eval.c:2960
# #14 0x66c02c in Rf_eval /data/gannet/ripley/R/svn/R-devel/src/main/eval.c:798
# #15 0x6ea51d in Rf_ReplIteration /data/gannet/ripley/R/svn/R-devel/src/main/main.c:264
# #16 0x6ea51d in Rf_ReplIteration /data/gannet/ripley/R/svn/R-devel/src/main/main.c:200
# #17 0x6eac18 in R_ReplConsole /data/gannet/ripley/R/svn/R-devel/src/main/main.c:314
# #18 0x6ead64 in run_Rmainloop /data/gannet/ripley/R/svn/R-devel/src/main/main.c:1113
# #19 0x6eadb2 in Rf_mainloop /data/gannet/ripley/R/svn/R-devel/src/main/main.c:1120
# #20 0x419368 in main /data/gannet/ripley/R/svn/R-devel/src/main/Rmain.c:29
# #21 0x7fc7c2bddf42 in __libc_start_main (/lib64/libc.so.6+0x23f42)
# #22 0x41baad in _start (/data/gannet/ripley/R/gcc-SAN/bin/exec/R+0x41baad)

