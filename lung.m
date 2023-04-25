%filename: lung.m (main program)
clear all
clf
global Pstar cstar n maxcount M Q camax RT cI;

setup_lung
cvsolve
outchecklung

%% Task 3
clear all
clf
global Pstar cstar n maxcount M Q camax RT cI;

inppOa = zeros(10,1);
mAlvppOa = zeros(10,1);
mArtppOa = zeros(10,1);
vppOa = zeros(10,1);
betaA = linspace(0,1,10);

for i = 1:10
    beta = i*0.1;
    setup_lung
    cvsolve
    outchecklung
    inppOa(i,1) = PI;
    mAlvppOa(i,1) = PAbar; % mean alveolar oxygen partial pressure
    mArtppOa(i,1) = Pabar; % mean arterial oxygen partial pressure
    vppOa(i,1) = Pv;
end

figure
plot(betaA,inppOa)
hold on
plot(betaA,mAlvppOa)
plot(betaA,mArtppOa)
plot(betaA,vppOa)
hold off
title('beta v.s. Partial Pressure of Oxygen')
legend('Inspired','Alveolar','Arterial','Venous')
xlabel('beta')
ylabel('Partial Pressure of Oxygen (mmHg)')

%% Task 4
clear all
clf
global Pstar cstar n maxcount M Q camax RT cI;

betaA = linspace(0,1,10);
M_A = zeros(10,1);
cref=0.2/(22.4*(310/273));

count = 1;
for i = 1:10
    beta = i*0.1;
    M=0.25*cref*5.6;
    while i>0
        setup_lung % setup lung with new M
        try %run cvsolve, if error, then break
            cvsolve
        catch %if the error is thrown out, then break the while loop
            break
        end
        M= M+0.0005;
    end
    M_A(i,1) = M;
    count = count + 1;
end
figure
plot(betaA,M_A)
title('Maximum sustainable rate of oxygen consumption for different Betas')
xlabel('beta')
ylabel('M, rate of oxygen consumption (moles/minute)')

%% Task 5
clear all
clf
global Pstar cstar n maxcount M Q camax RT cI;

% beta set to 0.5
cI = 0.0079; % initial value if cI = cref

cI_A = zeros(100,1); %inspired air oxygen concentration
pV_A = zeros(100,1); %oxygen partial pressure in venous blood
pa_A = zeros(100,1); %mean arterial oxygen partial pressure
pA_A = zeros(100,1); %mean alveolar oxygen partial pressure

for i=1:100
    setup_lung
    cvsolve
    outchecklung
    cI_A(i,1) = cI; % update arrays with values from outchecklung.m
    pV_A(i,1) = Pv;
    pa_A(i,1) = Pabar;
    pA_A(i,1) = PAbar;
    cI = cI + 0.0001;
end
figure
plot(cI_A,pV_A,'.')
hold on
plot(cI_A,pa_A,'.')
plot(cI_A,pA_A,'.')
hold off
legend('Venous Blood','Mean Arterial','Mean Alveolar')
title('cI v.s. Oxygen Partial Pressure')
xlabel('cI (moles per liter)')
ylabel('Oxygen Partial Pressure (mmHg)')

%% Task 6
clear all
clf
global Pstar cstar n maxcount M Q camax RT cI;

%cref=0.2/(22.4*(310/273));
%RT=760*22.4*(310/273);
%cI=cref;
%PI = cI*cref;

PI=152; %initialize PI with 152 mmHg, no change in elevation

pI_A = zeros(1000,1); % inspired air partial pressure of oxygen
pV_A = zeros(1000,1); %oxygen partial pressure in venous blood
pa_A = zeros(1000,1); %mean arterial oxygen partial pressure
pA_A = zeros(1000,1); %mean alveolar oxygen partial pressure

% Define cvzero: normal resting value of the venous oxygen concentration
% Both cardiac output and total alveolar ventilation increase during exercise. 
% We can make a simplified model of this by 
% assuming that both quantities are inversely proportional 
% to the venous oxygen concentration
% cvzero is lung evaluated at beta = 0

cvzero = 0.0058;
cv_A = zeros(1000,1);

for i=1:1000
    setup_lung
    try
        cvsolve
    catch
        break
    end
    outchecklung
    if (cvzero-cv) > 0.001 % deviation away from resting oxygen consumption
        break
    end
    cv_A(i,1) = cv; % update arrays with values from outchecklung.m
    pI_A(i,1) = PI;
    pV_A(i,1) = Pv;
    pa_A(i,1) = Pabar;
    pA_A(i,1) = PAbar;
    PI = PI-1;
end
% remove all nonzero elements from arrays, 
% allows for plotting if an error cuts the loop short
pI_A = nonzeros(pI_A); % remove all nonzero elements from arrays
pV_A = nonzeros(pV_A);
pa_A = nonzeros(pa_A);
pA_A = nonzeros(pA_A);

figure
plot(pI_A,pV_A,'.')
hold on
plot(pI_A,pa_A,'.')
plot(pI_A,pA_A,'.')
hold off
legend('Venous Blood','Mean Arterial','Mean Alveolar')
title('pI (Altitude) v.s. Oxygen Partial Pressure in Blood')
xlabel('pI (mmHg)')
ylabel('Partial Pressure in Blood (mmHg)')

%% Task 7
clear all
clf
global Pstar cstar n maxcount M Q camax RT cI;

%cref=0.2/(22.4*(310/273));
%RT=760*22.4*(310/273);
%cI=cref;
%PI = cI*cref;

PI=152; %initialize PI with 152 mmHg, no change in elevation

pI_A = zeros(1000,1); % inspired air partial pressure of oxygen
pV_A = zeros(1000,1); %oxygen partial pressure in venous blood
pa_A = zeros(1000,1); %mean arterial oxygen partial pressure
pA_A = zeros(1000,1); %mean alveolar oxygen partial pressure

% Define cvzero: normal resting value of the venous oxygen concentration
% Both cardiac output and total alveolar ventilation increase during exercise. 
% We can make a simplified model of this by 
% assuming that both quantities are inversely proportional 
% to the venous oxygen concentration
% cvzero is lung evaluated at beta = 0

cvzero = 0.0099; %cvzero value adjusted for cstar=1.5*cref
cv_A = zeros(1000,1);

for i=1:1000
    setup_lung
    try
        cvsolve
    catch
        break %break in case M is too large and throws an error
    end
    outchecklung
    if (cvzero-cv) > 0.001 % deviation away from resting oxygen consumption
        break
    end
    cv_A(i,1) = cv; % update arrays with values from outchecklung.m
    pI_A(i,1) = PI;
    pV_A(i,1) = Pv;
    pa_A(i,1) = Pabar;
    pA_A(i,1) = PAbar;
    PI = PI-1;
end
% remove all nonzero elements from arrays, 
% allows for plotting if an error cuts the loop short
pI_A = nonzeros(pI_A); % remove all nonzero elements from arrays
pV_A = nonzeros(pV_A);
pa_A = nonzeros(pa_A);
pA_A = nonzeros(pA_A);

figure
plot(pI_A,pV_A,'.')
hold on
plot(pI_A,pa_A,'.')
plot(pI_A,pA_A,'.')
hold off
legend('Venous Blood','Mean Arterial','Mean Alveolar')
title('pI (Altitude) v.s. Oxygen Partial Pressure in Blood')
xlabel('pI (mmHg)')
ylabel('Partial Pressure in Blood (mmHg)')

%% Task 8
clear all
clf
global Pstar cstar n maxcount M Q camax RT cI;

%cref=0.2/(22.4*(310/273));
%RT=760*22.4*(310/273);
%cI=cref;
%PI = cI*cref;

PI=152; %initialize PI with 152 mmHg, no change in elevation

pI_A = zeros(1000,1); % inspired air partial pressure of oxygen
pV_A = zeros(1000,1); %oxygen partial pressure in venous blood
pa_A = zeros(1000,1); %mean arterial oxygen partial pressure
pA_A = zeros(1000,1); %mean alveolar oxygen partial pressure

% Define cvzero: normal resting value of the venous oxygen concentration
% Both cardiac output and total alveolar ventilation increase during exercise. 
% We can make a simplified model of this by 
% assuming that both quantities are inversely proportional 
% to the venous oxygen concentration
% cvzero is lung evaluated at beta = 0

cvzero = 0.0058;
cv_A = zeros(1000,1);

for i=1:1000
    setup_lung
    try
        cvsolve
    catch
        break %break in case M is too large and throws an error
    end
    outchecklung
    if (cvzero-cv) > 0.001 % deviation away from resting oxygen consumption
        break
    end
    cv_A(i,1) = cv; % update arrays with values from outchecklung.m
    pI_A(i,1) = PI;
    pV_A(i,1) = Pv;
    pa_A(i,1) = Pabar;
    pA_A(i,1) = PAbar;
    PI = PI-1;
end
% remove all nonzero elements from arrays, 
% allows for plotting if an error cuts the loop short
pI_A = nonzeros(pI_A); % remove all nonzero elements from arrays
pV_A = nonzeros(pV_A);
pa_A = nonzeros(pa_A);
pA_A = nonzeros(pA_A);

figure
plot(pI_A,pV_A,'.')
hold on
plot(pI_A,pa_A,'.')
plot(pI_A,pA_A,'.')
hold off
legend('Venous Blood','Mean Arterial','Mean Alveolar')
title('pI (Altitude) v.s. Oxygen Partial Pressure in Blood')
xlabel('pI (mmHg)')
ylabel('Partial Pressure in Blood (mmHg)')

%% Task 9

clear all
clf
global Pstar cstar n maxcount M Q camax RT cI;

cstar=0.2/(22.4*(310/273)); %initialize cstar value

pI_A = zeros(200,1); % inspired air partial pressure of oxygen
pV_A = zeros(200,1); % oxygen partial pressure in venous blood
pa_A = zeros(200,1); % mean arterial oxygen partial pressure
pA_A = zeros(200,1); % mean alveolar oxygen partial pressure

cabar_A = zeros(200,1); % mean arterial oxygen concentration
cAbar_A = zeros(200,1); % mean alveolar oxygen concentration

cstar_A = zeros(200,1); % cstar collection array

cvzero = 0.0058; % initialize cvzero to use during deviation from resting
cv_A = zeros(200,1); % venous blood oxygen concentration

for i=1:200
    setup_lung
    try
        cvsolve
    catch
        break %break in case M is too large and throws an error
    end
    outchecklung
    if (cvzero-cv) > 0.001 % deviation away from resting oxygen consumption
        break
    end
    cv_A(i,1) = cv; % update arrays with values from outchecklung.m
    pI_A(i,1) = PI;
    pV_A(i,1) = Pv;
    pa_A(i,1) = Pabar;
    pA_A(i,1) = PAbar;

    cabar_A(i,1) = cabar;
    cAbar_A(i,1) = cAbar;
    cstar_A(i,1) = cstar;
    cstar = cstar - 0.000005;
end

% remove all nonzero elements from arrays, 
% allows for plotting if an error cuts the loop short
cstar_A = nonzeros(cstar_A);

% Pressures
pI_A = nonzeros(pI_A); % remove all nonzero elements from arrays
pV_A = nonzeros(pV_A);
pa_A = nonzeros(pa_A);
pA_A = nonzeros(pA_A);

figure
plot(cstar_A,pV_A,'.')
hold on
plot(cstar_A,pa_A,'.')
plot(cstar_A,pA_A,'.')
plot(cstar_A,pI_A,'.')
hold off
legend('Venous Blood','Mean Arterial','Mean Alveolar','Inspired')
title('c-star v.s. Oxygen Partial Pressure in Blood')
xlabel('c-star (moles/liter)')
ylabel('Partial Pressure in Blood (mmHg)')

% Blood concentration
cabar_A = nonzeros(cabar_A); % remove all nonzero elements from arrays
cAbar_A = nonzeros(cAbar_A);
cv_A = nonzeros(cv_A);

figure
plot(cstar_A,cabar_A,'.')
hold on
plot(cstar_A,cAbar_A,'.')
plot(cstar_A,cv_A,'.')
hold off
legend('Mean Arterial','Mean Alveolar','Venous Blood')
title('c-star v.s. Oxygen Concentration in Blood')
xlabel('c-star (moles/liter)')
ylabel('Concentration in Blood (moles/liter)')

%% Task 10
clear all
clf
global Pstar cstar n maxcount M Q camax RT cI;

%cref=0.2/(22.4*(310/273));
%RT=760*22.4*(310/273);
%cI=cref;
%PI = cI*cref;

% beta defined at beta = 0.5

PI=152; %initialize PI with 152 mmHg, no change in elevation
cstar = 0.0072; %initialize cstar with cref BUT CHANGE FOR TASK 10

pI_A = zeros(300,1); % inspired air partial pressure of oxygen
pV_A = zeros(300,1); %oxygen partial pressure in venous blood
pa_A = zeros(300,1); %mean arterial oxygen partial pressure
pA_A = zeros(300,1); %mean alveolar oxygen partial pressure

% Define cvzero: normal resting value of the venous oxygen concentration
% Both cardiac output and total alveolar ventilation increase during exercise. 
% We can make a simplified model of this by 
% assuming that both quantities are inversely proportional 
% to the venous oxygen concentration
% cvzero is lung evaluated at beta = 0

cvzero = 0.0058;
cv_A = zeros(300,1);

for i=1:300
    setup_lung
    try
        cvsolve
    catch
        break
    end
    outchecklung
    if (cvzero-cv) > 0.001 % deviation away from resting oxygen consumption
        break
    end
    cv_A(i,1) = cv; % update arrays with values from outchecklung.m
    pI_A(i,1) = PI;
    pV_A(i,1) = Pv;
    pa_A(i,1) = Pabar;
    pA_A(i,1) = PAbar;
    PI = PI-1;
end
% remove all nonzero elements from arrays, 
% allows for plotting if an error cuts the loop short
pI_A = nonzeros(pI_A); 
pV_A = nonzeros(pV_A);
pa_A = nonzeros(pa_A);
pA_A = nonzeros(pA_A);

figure
plot(pI_A,pV_A,'.')
hold on
plot(pI_A,pa_A,'.')
plot(pI_A,pA_A,'.')
hold off
legend('Venous Blood','Mean Arterial','Mean Alveolar')
title('pI (Altitude) v.s. Oxygen Partial Pressure in Blood')
xlabel('pI (mmHg)')
ylabel('Partial Pressure in Blood (mmHg)')
