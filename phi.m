function c=phi(ca,cv,r)
%filename:phi.m
global RT cI;
c=ca-cv+r.*(H(ca)/RT-cI);
