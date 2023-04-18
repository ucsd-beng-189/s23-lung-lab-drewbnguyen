Lung Lab

Model requires: carterial.m, cvsolve.m, H.m, lung.m, Mdiff.m, outchecklung.m, phi.m, setup_lung.m

Instructions: To complete this assignment you must do 2 things. 
1)Answer the questions below in a Word document. Generate plots/measurements as necessary to explain your findings

2)Commit your code changes for each step and push to Github classroom. 
a. Clone the repository from the Github classroom.
b.For each task make sure to commit both changes to the provided code as well as any scripts you generate. Also add figures associated with your findings to the repo as .jpg or .png files.


Task 1: Describe the following files
•Describe in words what setup_lung.m plots on lines 60-67
•Describe in words what cvsolve.m does
•Compare and contrast carterial.m and cvsolve.m


Task 2: The programs introduced in the foregoing section make it possible to study the consequences of changing the extent to which airflow and bloodflow are coordinated throughout the lung. What does beta = 0 represent in the model? Similarly, what does beta=1 represent? Hint: look at lines 60-67 in setup_lung.m


Task 3: Modify the code so you can run the experiments with various beta values. Plot the inspired partial pressure of oxygen mean alveolar partial pressure of oxygen, the mean arterial partial pressure of oxygen, and the venous partial pressure of oxygen, all on the same graph as functions of beta. Explain the results in terms of the theory developed in this chapter.


Task 4: We can also ask, what is the maximum sustainable rate of oxygen consumption for a given set of parameters, and how does it depend on the extent to which ventilation is matched to perfusion? To study this, use trial and error to find the maximum sustainable rate of oxygen consumption for each of several values of beta, and plot this rate as function of beta. (Recall that the program will stop with an error message if the assumed rate of oxygen consumption is too high.) Keep all parameters other than beta and M constant throughout this study.


Task 5: Now consider the effects of altitude, with and without adaptation. The primary change at altitude is a reduction in the inspired partial pressure of oxygen, and also, therefore, in the corresponding concentration of oxygen in the inspired air. You can use the inspired concentration or partial pressure of oxygen as your principal independent parameter for this study. First vary cI with all other parameters constant to see the effects of acute exposure to altitude. Plot the mean alveolar, mean arterial, and venous partial pressures and concentrations of oxygen as functions of the cI.


Task 6: Convert your plots such that your plot the results as a function of altitude (instead of cI). At what altitude does it become impossible to sustain the normal resting rate of oxygen consumption (without breathing harder or increasing the cardiac output)? These results are for a person adapted to sea level.


Task 7: Now let’s evaluate someone who has been living at high enough altitude for so long that his/her red blood cells make up 60% of the blood volume instead of the usual 40%. To do so, set cstar=l.5*cref instead of cstar=cref. At what altitude does it become impossible even for this altitude-adapted person to sustain the normal resting rate of oxygen consumption (without breathing harder or increasing cardiac output)? Compare the results obtained for this altitude-adapted person to those obtained above for a person who lacks such adaptation to altitude.


Task 8: Say I (someone used to living at sea level), fly to Peru and unfortunately experience get a pulmonary embolism on the flight. This cause regions of my lungs to be hypoperfused. What is the effect on ventilation-perfusion matching? Revist your results in Task 6 for this case.


Task 9: Another condition involving a change in the hemoglobin concentration in blood is anemia. To simulate anemia in our model, simply reduce cstar while keeping all other parameters constant. Plot the inspired, mean alveolar, mean arterial, and venous partial pressures and concentrations of oxygen as functions of cstar as that parameter is reduced below cref while all other parameters are held at their normal resting (sea-level) values.

What is the minimum value of cstar at which the normal resting rate of oxygen consumption can be maintained (without breathing harder or increasing cardiac output). How do these results depend on the ventilation-perfusion mismatch parameter beta? (Repeat the study for several different beta to find out.)


Task 10: Revisit your findings in Task 6 for a patient who is anemic.
