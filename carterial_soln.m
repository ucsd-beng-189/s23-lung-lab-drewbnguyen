function ca=carterial(cv,r)
%filename: carterial.m
global camax n maxcount;
%initialize:
ca_lower=zeros(n,1);         %lower end of interval
ca_upper= ones(n,1)*camax;   %upper end of interval
ca=(ca_lower+ca_upper)/2;    %midpoint of interval
%start main loop:
for count=1:maxcount
  %check sign of function at midpoint of interval:
  pos=(phi(ca,cv,r)>0); 
  neg=~pos;
  %if (phi>0) ca_upper=ca; else ca_lower=ca:             
  ca_upper=ca_upper.*neg+ca.*pos;    
  ca_lower=ca_lower.*pos+ca.*neg;
  %ca=midpoint of new interval:    
  ca=(ca_lower+ca_upper)/2;
end
