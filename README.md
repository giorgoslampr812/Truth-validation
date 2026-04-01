Brief description 

Rivet routine that works on run 2 EVNT samples with Rivet version 3.1.6

VBS_fid.cc - Rivet routine to make llqq ntuples. 

analyse.C - Code for processing of llqq ntuples - gives reasonable results for ratio of resolved/merged regions (page 651 of https://cds.cern.ch/record/2808769). 

VBS_lv.cc - Rivet routine to make lvqq ntuples. 

analyse_1lep.C - Code for processing of lvqq ntuples - gives reasonable results for ratio of resolved/merged regions. 

pyt_script.py - Python script example to run the different rivet routines. 

Running example 

setupATLAS -c centos7 

cd /data/glamprinoudis/RIVET

asetup AthGeneration,22.6.24,here

source setupRivet

rivet-build RivetVBS_fid.so VBS_fid.cc $(root-config --cflags --libs) 

athena pyt_script.py


