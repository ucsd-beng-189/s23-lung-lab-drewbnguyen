%filename: lung.m (main program)
clear all
clf
global Pstar cstar n maxcount M Q camax RT cI;
setup_lung
cvsolve
outchecklung
