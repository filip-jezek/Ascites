# Ascites
Ascites modeling, following the paper from Levitt and Levitt (Levitt, David G., and Michael D. Levitt. 2012. “Quantitative Modeling of the Physiology of Ascites in Portal Hypertension.” BMC Gastroenterology 12 (March): 26.), coupled with hemodynamics. Please see the model for the documentation.

One needs a recent Physiolibrary (https://github.com/filip-jezek/Physiolibrary/tree/v2.4_maintenance) to run the model.

## Model description
We have employed 0D modelling of the splanchnic hemodynamics. We have used the steady state Levitt’s and Levitt’s model (Levitt and Levitt 2012) to calculate the steady state ascites pressure. This model consists of resistive components in series (intestinal arteries, intestinal veins, liver, hepatic vein) and a model of peritoneal compartment. When the hepatic vein pressure is lower than the ascites pressure, the hepatic vein collapses, increasing its pressure drop.
The fixed pressure drops assumed for the liver (6 mmHg), intestine venules (3 mmHg) and intestinal arteries were replaced with resistances so that the pressure drop at the said components stays the same for nominal inflow of 1 L/min. Thus, we can observe the change of the pressure drop at different inflows.

### Natural Shunts
In patients with portal hypertension there are often observed shunts from portal vein into the vena cava or into the hepatic vein. We hypothesize that the elasticity of the venous wall, probably together with the long-term remodeling, causes otherwise tiny vessels to grow in diameter and by increasing its lumen to be able to diverge large amount of portal vein flow around the liver. By being parallel to the liver, the HVPG as well as portal vein pressure drops. 
The shunt is characterized by its linear resistance R. 

dP = R*Q

The resistance R of the shunt is given by Hagen-Poiseulle’s law:

R = 8*ni*L/(pi*r^4);

Where ni is dynamic viscosity of the blood (4e-3), L length of the shunt and r is the actual diameter. The actual diameter is calculated from the actual volume:

V = L*r

The actual volume is then given by the midpoint pressure P_inner = (p_in + P_out)/2 and assumed vessel compliance C

V = P_inner*C + V0

Where V0 is the zero pressure volume, given by zero pressure diameter

V0 = L*r0

We have estimated pressure drop across the TIPS and calculated its resistance.


### The liver venous tree
To evaluate the impact of portal vein embolism, we have created a simple portal vein tree, branching into 4 major liver segments: right Posterior and Anterior and left lateral and medial segments. The resistances are assumed to be proportional to the volume fractions of particular segments from (Leelaudomlipi et al. 2002).

### Experiment setup
For the first experiment, we fed the splanchnic circulation with a constant inflow rate of 1 L/min without any shunt. We then compared the HVPG for cases with and without the shunt applied in parallel to the liver and additionally with the TIPS applied.

For the embolism scenario, we have tested obstruction in the right anterior segment and compared healthy with ascitic with and without natural shunts.



