#!/bin/bash
g++ run.C -o run.exe \
      `/net/hisrv0001/home/abaty/jetGrooming/CMSSW_9_0_0/src/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` -lfastjetcontribfragile $(root-config --cflags --libs)
