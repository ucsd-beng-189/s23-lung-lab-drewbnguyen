%filename: cvsolve.m
%This script solves the equation 
%Mdiff(cv,r)=0 
%by the method of bisection.
if(Mdiff(0,r)>0)
  %stop program and issue error message:
  error('M is too large')
end
%initialize:
cv_lower=0;            %lower limit of interval
cv_upper=camax;        %upper limit of interval
cv=(cv_lower+cv_upper)/2; %midpoint of interval
%begin main loop:
for count=1:maxcount
  %check sign of Mdiff at midpoint of interval
  pos=(Mdiff(cv,r)>0);
  neg=~pos;
  %if(pos) cv_upper=cv; else cv_lower=cv:
  cv_upper=cv_upper*neg+cv*pos 
  cv_lower=cv_lower*pos+cv*neg 
  %set cv=midpoint of current interval:
  cv=(cv_lower+cv_upper)/2;
end
