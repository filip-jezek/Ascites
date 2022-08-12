# Ascites simulator


<bdl-fmi id="idfmi" mode="oneshot" src="Lymphatics_0Hemodynamics_0Experiments_0HVPGShuntsComparison.js" fminame="Lymphatics_0Hemodynamics_0Experiments_0HVPGShuntsComparison" tolerance="0.000001" starttime="0" fstepsize="1" stoptime="30" guid="{9dba4a4b-729a-40d8-8c81-0bd98c12ce55}" valuereferences="905969715,905969733,905969714,905969765,905969776,905969768,905969767,100663313" valuelabels="ascites_NoShunts.HVPG,ascites_Shunts.HVPG,ascites_NoShunts.PPV,ascites_Shunts.PPV,ascites_Shunts.shunt.d,ascites_Shunts.Q_shunt,ascites_Shunts.Q_liver,ascites_NoShunts.Q_liver" inputs="id1,16777243,7.5006e-09,1,t;id2,16777270,1.6666666666667e-05,1,t" inputlabels="ascites_Shunts.shunt.Comp;Inflow" eventlisten="change"></bdl-fmi>

<bdl-range id="id1" title="Remodeling sensitivity" min="0.1" max="5" default="1" step="0.1"></bdl-range>

<bdl-range id="id2" title="Inflow" min="0.7" max="1.3" default="1" step="0.02"></bdl-range>

<div class="w3-row">
<div class="w3-half">
HVPG

<bdl-chartjs-time width="400" height="400" fromid="idfmi" labels="ascites no shunts,ascites shunts,ascites shuntstiff" initialdata="0, 1" refindex="0" refvalues="2" ylabel="HVPG (mmHg)" xlabel="Liver resistance (mmHg.min/L)" showLine="false" convertors="1,133.32;1,133.32" min="0" max="35"></bdl-chartjs-time>

</div>
<div class="w3-half">
PPV

<bdl-chartjs-time width="400" height="400" fromid="idfmi" labels="ascites no shunts,ascites shunts,ascites shuntstiff" initialdata="0, 1, 2" refindex="2" refvalues="2" ylabel="PPV (mmHg)" xlabel="Liver resistance (mmHg.min/L)" convertors="1,133.32;1,133.32;1,133.32" min=0 max=70></bdl-chartjs-time>

</div>
</div>

<div class="w3-row">
<div class="w3-half">
Flows

<bdl-chartjs-time width="400" height="400" fromid="idfmi" labels="Shunts, Liver, Total inflow" initialdata="0, 1" refindex="5" refvalues="3" ylabel="Flow (L/min)" xlabel="Liver resistance (mmHg.min/L)" showLine="false" convertors="6e4,1;6e4,1;6e4,1" min="0" max="1.3"></bdl-chartjs-time>

</div>
<div class="w3-half">
Shunt diameter

<bdl-chartjs-time width="400" height="400" fromid="idfmi"  initialdata="0," refindex="4" refvalues="1" ylabel="Diameter(mm)" xlabel="Liver resistance (mmHg.min/L)" convertors="1000,1" min="0" max="6"></bdl-chartjs-time>

</div>
</div>
