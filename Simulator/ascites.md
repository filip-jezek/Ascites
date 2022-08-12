Ascites


<bdl-fmi id="idfmi" mode="oneshot" src="Lymphatics_0Hemodynamics_0Experiments_0HVPGShuntsComparison.js" fminame="Lymphatics_0Hemodynamics_0Experiments_0HVPGShuntsComparison" tolerance="0.000001" starttime="0" fstepsize="1" stoptime="30" guid="{9dba4a4b-729a-40d8-8c81-0bd98c12ce55}" valuereferences="905969715,905969733,905969810,16777271" valuelabels="ascites_NoShunts.HVPG,ascites_Shunts.HVPG,ascites_ShuntStiff.HVPG,Shunt_Compliance" inputs="id1,16777243,7.5006e-09,1,t" inputlabels="ascites_Shunts.shunt.Comp" eventlisten="change"></bdl-fmi>

<bdl-range id="id1" title="Remodeling sensitivity" min="0.1" max="5" default="1" step="0.1"></bdl-range>

<bdl-chartjs-time width="400" height="400" fromid="idfmi" labels="ascites no shunts,ascites shunts,ascites shuntstiff" initialdata="0, 1, 2" refindex="0" refvalues="3" ylabel="HVPG (mmHg)" xlabel="Liver resistance" convertors="1,133.32;1,133.32;1,133.32"></bdl-chartjs-time>
