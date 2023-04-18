function Md=Mdiff(cv,r)
%filename: Mdiff.m
global M Q;
%Md=M-sum(i=1:n)(Q(i)*(carterial(cv,r(i))-cv)):
Md=M-Q'*(carterial(cv,r)-cv)
