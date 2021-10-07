within ;
package Lymphatics
  model Ascites
    parameter Physiolibrary.Types.OsmoticPermeability L_Y=1.250102626409427e-10
        *(7.86/60);
        parameter Physiolibrary.Types.Pressure P_min=266.64477483
                                                     "nominal abodminal cavity pressure";
    Physiolibrary.Hydraulic.Sources.UnlimitedPump SplanchnicInflow(SolutionFlow=8.3333333333333e-06)
      annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));
    Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=266.64477483)
      annotation (Placement(transformation(extent={{100,-100},{80,-80}})));
    Physiolibrary.Hydraulic.Components.PumpPressureHead HV(p_head=266.64477483)
      "Hepatic vein pressure gradient"
      annotation (Placement(transformation(extent={{60,-100},{40,-80}})));
    Physiolibrary.Hydraulic.Components.PumpPressureHead HPVG(p_head=533.28954966)
      "Hepatic vein pressure gradient"
      annotation (Placement(transformation(extent={{20,-100},{0,-80}})));
    Physiolibrary.Hydraulic.Sensors.PressureMeasure P_HV
      "Pressure in hepatic vein"
      annotation (Placement(transformation(extent={{34,-78},{54,-58}})));
    Physiolibrary.Hydraulic.Sensors.PressureMeasure P_PV "portal vein pressure"
      annotation (Placement(transformation(extent={{-8,-78},{12,-58}})));
    Physiolibrary.Hydraulic.Components.PumpPressureHead P_IC(p_head=
          399.967162245)
      "Intestinal capillaries pressure drop"
      annotation (Placement(transformation(extent={{-40,-100},{-60,-80}})));
    Physiolibrary.Osmotic.Components.Membrane M_i(cond=1.250102626409427e-10*(
          6.25/60),
                useHydraulicPressureInputs=true)
      "Instestinal capillary lymph flow"
      annotation (Placement(transformation(extent={{-64,40},{-44,20}})));
    Physiolibrary.Osmotic.Components.OsmoticCell abodminalCompartment(
      NumberOfMembraneTypes=1,
      useImpermeableSolutesInput=false,
      ImpermeableSolutes={0.0001},
      volume(start=0.0002))
      annotation (Placement(transformation(extent={{48,20},{68,40}})));
    Physiolibrary.Osmotic.Sources.UnlimitedSolution Osm_P(Osm(displayUnit="mmol/l")=
           1) "Plasma osmolarity"
      annotation (Placement(transformation(extent={{-100,20},{-80,40}})));
    Physiolibrary.Hydraulic.Sensors.PressureMeasure P_C
      "Intestinal capillary pressure"
      annotation (Placement(transformation(extent={{-70,-78},{-50,-58}})));
    AbdominalCompliance abdominalCompliance(
      V0=0.0001,
      P0=P_min,
      D=6.0004926067653e-06)
      annotation (Placement(transformation(extent={{58,0},{38,20}})));
    Physiolibrary.Osmotic.Components.Membrane M_L(
      cond=1.250102626409427e-10*(10.3/60),
                useHydraulicPressureInputs=true,
      pBreak=0)                                  "Liver capillary lymph flow"
      annotation (Placement(transformation(extent={{-50,0},{-30,-20}})));
    Modelica.Blocks.Math.Add3 add3_1(
                                 k1=0.5, k2=0.5,
      k3=-1)
      annotation (Placement(transformation(extent={{-6,-28},{-26,-48}})));
    Physiolibrary.Osmotic.Sources.SolventOutflux solventOutflux(
        useSolutionFlowInput=true)
      annotation (Placement(transformation(extent={{108,20},{128,40}})));
    Modelica.Blocks.Sources.RealExpression realExpression(y=max(0, L_Y*(
          abdominalCompliance.y - CVP.p + 2)))
      annotation (Placement(transformation(extent={{154,40},{134,60}})));
    Physiolibrary.Chemical.Sources.UnlimitedSolutionStorage
      unlimitedSolutionStorage(Conc=1)
      annotation (Placement(transformation(extent={{-100,60},{-80,80}})));
    Physiolibrary.Chemical.Components.Substance substance(useNormalizedVolume=
          false, solute_start=2e-05)
      annotation (Placement(transformation(extent={{30,60},{50,80}})));
    Physiolibrary.Osmotic.Sensors.FlowMeasure flowMeasure
      annotation (Placement(transformation(extent={{82,40},{102,20}})));
    Physiolibrary.Chemical.Components.Clearance clearance(useSolutionFlowInput=
          true)
      annotation (Placement(transformation(extent={{82,80},{102,60}})));
    Physiolibrary.Osmotic.Sensors.FlowMeasure flowMeasure1
      annotation (Placement(transformation(extent={{2,40},{22,20}})));
    Physiolibrary.Chemical.Components.Diffusion diffusion(Conductance(
          displayUnit="l/day") = 2.3148148148148e-08)
      annotation (Placement(transformation(extent={{0,60},{20,80}})));
    Modelica.Blocks.Sources.RealExpression Pi_P(y=M_i.opi)
      annotation (Placement(transformation(extent={{80,-20},{100,0}})));
    Modelica.Blocks.Sources.RealExpression Pi_A(y=M_i.opo)
      annotation (Placement(transformation(extent={{80,-40},{100,-20}})));
    Physiolibrary.Types.Constants.PressureConst pressure(k=P_min)
      annotation (Placement(transformation(extent={{14,-34},{6,-26}})));
    Modelica.Blocks.Sources.BooleanExpression ascites(y=abdominalCompliance.y
           > P_HV.pressure)
      annotation (Placement(transformation(extent={{80,0},{100,20}})));
  equation
    connect(SplanchnicInflow.q_out, P_IC.q_out) annotation (Line(
        points={{-80,-90},{-60,-90}},
        color={0,0,0},
        thickness=1));
    connect(P_IC.q_in, HPVG.q_out) annotation (Line(
        points={{-40,-90},{0,-90}},
        color={0,0,0},
        thickness=1));
    connect(HPVG.q_in, HV.q_out) annotation (Line(
        points={{20,-90},{40,-90}},
        color={0,0,0},
        thickness=1));
    connect(HV.q_in, CVP.y) annotation (Line(
        points={{60,-90},{80,-90}},
        color={0,0,0},
        thickness=1));
    connect(P_HV.q_in, HV.q_out) annotation (Line(
        points={{40,-74},{30,-74},{30,-90},{40,-90}},
        color={0,0,0},
        thickness=1));
    connect(P_PV.q_in, HPVG.q_out) annotation (Line(
        points={{-2,-74},{-12,-74},{-12,-90},{0,-90}},
        color={0,0,0},
        thickness=1));
    connect(Osm_P.port, M_i.q_in) annotation (Line(
        points={{-80,30},{-64,30}},
        color={127,127,0},
        thickness=1));
    connect(P_C.q_in, P_IC.q_out) annotation (Line(
        points={{-64,-74},{-64,-90},{-60,-90}},
        color={0,0,0},
        thickness=1));
    connect(P_C.pressure, M_i.hydraulicPressureIn)
      annotation (Line(points={{-54,-72},{-54,22},{-62,22}}, color={0,0,127}));
    connect(abodminalCompartment.volume, abdominalCompliance.u) annotation (Line(
          points={{64,20},{72,20},{72,10},{60,10}}, color={0,0,127}));
    connect(M_i.hydraulicPressureOut, abdominalCompliance.y)
      annotation (Line(points={{-46,22},{-46,10},{37,10}}, color={0,0,127}));
    connect(add3_1.u2, P_HV.pressure) annotation (Line(points={{-4,-38},{60,-38},
            {60,-72},{50,-72}}, color={0,0,127}));
    connect(add3_1.u1, P_PV.pressure) annotation (Line(points={{-4,-46},{18,-46},
            {18,-72},{8,-72}}, color={0,0,127}));
    connect(M_L.q_in, M_i.q_in) annotation (Line(
        points={{-50,-10},{-72,-10},{-72,30},{-64,30}},
        color={127,127,0},
        thickness=1));
    connect(M_i.q_out, M_L.q_out) annotation (Line(
        points={{-44,30},{-4,30},{-4,-10},{-30,-10}},
        color={127,127,0},
        thickness=1));
    connect(add3_1.y, M_L.hydraulicPressureIn) annotation (Line(points={{-27,
            -38},{-48,-38},{-48,-18}}, color={0,0,127}));
    connect(M_L.hydraulicPressureOut, abdominalCompliance.y) annotation (Line(
          points={{-32,-18},{0,-18},{0,10},{37,10}},
          color={0,0,127}));
    connect(solventOutflux.solutionFlow, realExpression.y) annotation (Line(
          points={{118,37},{122,37},{122,50},{133,50}}, color={0,0,127}));
    connect(abodminalCompartment.volume, substance.solutionVolume) annotation (
        Line(points={{64,20},{72,20},{72,74},{36,74}}, color={0,0,127}));
    connect(substance.solute, abodminalCompartment.impermeableSolutes[1])
      annotation (Line(points={{46,60},{46,36},{50,36}}, color={0,0,127}));
    connect(flowMeasure.q_out, solventOutflux.q_in) annotation (Line(
        points={{102,30},{112,30}},
        color={127,127,0},
        thickness=1));
    connect(clearance.solutionFlow, flowMeasure.volumeFlowRate)
      annotation (Line(points={{92,63},{92,38}}, color={0,0,127}));
    connect(substance.q_out, clearance.q_in) annotation (Line(
        points={{40,70},{82,70}},
        color={107,45,134},
        thickness=1));
    connect(flowMeasure.q_in, abodminalCompartment.q_in[1]) annotation (Line(
        points={{82,30},{58,30}},
        color={127,127,0},
        thickness=1));
    connect(flowMeasure1.q_out, abodminalCompartment.q_in[1]) annotation (Line(
        points={{22,30},{58,30}},
        color={127,127,0},
        thickness=1));
    connect(flowMeasure1.q_in, M_L.q_out) annotation (Line(
        points={{2,30},{-4,30},{-4,-10},{-30,-10}},
        color={127,127,0},
        thickness=1));
    connect(diffusion.q_in, unlimitedSolutionStorage.q_out) annotation (Line(
        points={{0,70},{-80,70}},
        color={107,45,134},
        thickness=1));
    connect(diffusion.q_out, substance.q_out) annotation (Line(
        points={{20,70},{40,70}},
        color={107,45,134},
        thickness=1));
    connect(add3_1.u3, pressure.y)
      annotation (Line(points={{-4,-30},{5,-30}}, color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=36000,
        __Dymola_NumberOfIntervals=5000,
        Tolerance=1e-05,
        __Dymola_Algorithm="Dassl"));
  end Ascites;

  model AbdominalCompliance
    extends Modelica.Blocks.Interfaces.SISO(u(unit="m3"), y(unit="Pa"));

    parameter Physiolibrary.Types.Volume V0 "Nominal volume";
    parameter Physiolibrary.Types.Pressure P0 "Nominal pressure";
    parameter Physiolibrary.Types.HydraulicCompliance D = 800 "Compliance";
  equation
    u = (y - P0)*D + V0;


  end AbdominalCompliance;

  model LevittCase1SS "Model of dynamic ascites build-up by Levitt and Levitt (2012), PMID 22453061. Converted from Maple code."
  //   Real D,Ly,Pbreak,Va,R,F,Pra,Pp,Pc,Pa,Pi,Pl,Phv,Aa1,Aa2,Pa1,Pa2,Ai,Jl,Jla,Jy,m,Pmin,Jya,Ap,Pdel,Pgrad,nplot,maxLl,mPa,mLl,i,Ll,x1,x2,Papoint,Lypoint,Setupeqs,delLl,Lt,maxgrad,MaxPa,MaxAa,Pamin,mingrad,delgrad,maxN,Vpoint,Apoint,gradfract,m2Pa,Pgrad0,Vpoint2,delplot;
  //   Real D,Ly,Pbreak,Va,R,F,Pra,Pp,Pc,Pa,Pi,Pl,Phv,Aa1,Aa2,Pa1,Pa2,Jl,Jla,Jy,m,Pmin,Jya,Ap,Pdel,Pgrad,nplot,maxLl,mPa,mLl,i,Ll,x1,x2,Papoint,Lypoint,Setupeqs,delLl,Lt,maxgrad,MaxPa,MaxAa,Pamin,mingrad,delgrad,maxN,Vpoint,gradfract,m2Pa,Pgrad0,Vpoint2,delplot;
  //  Real Ji, Aa, Vmin;

  parameter Real Pra = 2 annotation(Evaluate = false);//right atrial pressure
  parameter Real Pamin =   2;//minimum ascites pressure when Jlymph = 0
  parameter Real Ap = 25; //Blood colloid osmotic pressure
  parameter Real m =  0.8;
  parameter Real Pmin = 2.0;//must be less than Pra
  parameter Real Pdel = 2.0 annotation(Evaluate=false);
  parameter Real Pbreak = 8;
  parameter Real Ly =  0.131; //ml/min/mm Hg
  parameter Real Ll =  0.172; //
  parameter Real Lt =  0.104;
  parameter Real D =  0.8; //New  - value of Henriksen and Lieberman
  parameter Real Vmin =  0.1;  //minimum volume when P = Pamin

  Real Phv;
  Real Pp;
  Real Pc;
  Real Pl;
  Real Jl;
  Real Ji;
  Real Jy;

  Real Pa, Aa;
  Real Pgrad = time annotation (Dialog(tab="General", group="Inputs"));

  Real Av = Vmin+D*(Pa-Pamin);
  equation

    //First condition, Pa&lt;Pra+Pdel;
    if Pa < Pra + Pdel then
      Phv = Pra+Pdel;
      Jy = if (Pmin+Pa-Pra<=0) then 0 else Ly*(Pmin+Pa-Pra);//0 Simple case Pa&gt;Pra

    else
      Phv = Pa;// Simple case - assume that there is ascites with high pressure - makes algabra simpler</Text-field>
      Jy = Ly*(Pmin+Pa-Pra);//0 Simple case Pa&gt;Pra</Text-field>

    end if;
  //Jy= max(0,Ly*(Pa - Pra + Pmin)) "Lymph flow, eq. 20";
  Pp = Phv+Pgrad;
  Pc = Pp+3;//intestinal capillary pressure
  Pl = (Pp+Phv)/2.0;//Liver sinusoid pressure
  Jl = Ll*(Pl-Pa-Pbreak);//Simple case, assume Pl-Pa &gt; Pbreak
  Ji = Lt*(Pc-Pa-Ap+Aa); //flux from intestine to ascites space

  Aa*Jy=m*Ap*Jl;//protein balance, //Pa =
  Jy=Jl+Ji;//fluid balance, //Aa =
  //else
      //Second condition, Pa&gt;Pra+Pdel;


  // Pp = Phv+Pgrad;
  // Pc = Pp+3;//intestinal capillary pressure
  // Pl = (Pp+Phv)/2.0;//Liver sinusoid pressure
  // Jl = Ll*(Pl-Pa-Pbreak);//Simple case, assume Pl-Pa &gt; Pbreak
  // Ji = Lt*(Pc-Pa-Ap+Aa); //flux from intestine to ascites space</Text-field>
  //Jy = Ly*(Pmin+Pa-Pra);//0 Simple case Pa&gt;Pra</Text-field>
  // eq1 = Aa*Jy=m*Ap*Jl;//protein balancefile:///C:/home/UMICH/ascites/Lymphatics.mo
  // eq2 = Jy=Jl+Ji;//fluid balance</Text-field>
  // solutions2 = solve({eq1,eq2},{Pa,Aa});



    annotation (experiment(
        StartTime=6,
        StopTime=25,
        __Dymola_NumberOfIntervals=5000,
        Tolerance=1e-05,
        __Dymola_Algorithm="Dassl"));
  end LevittCase1SS;

  model LevitCase2Dynamics "Model of dynamic ascites build-up by Levitt and Levitt (2012), PMID 22453061. Converted from Maple code."
  parameter Boolean paracentesis = false;

  //Independent variables
  parameter Real dayconv =  24.0 "hours/day   - thus paramteters in units of 1/day";
  parameter Real volconv = 0.001 "liter/ml";


  //volconv = 1.0 "liter/ml";
  parameter Real Ap = 25 "Plasma colloid osmoitc pressure";
  parameter Real Po =  2.0 "minimum pressure in the pressure Pa vs Volume relation";
  parameter Real Pvg = 2.0 "Hep[atic vein to CVP gradient";
  parameter Real Pmin = Po "minimum pressure in the lymph flow equation (in paper - assume that = Po";
  parameter Real Vmin =  0.1 "No lymph flow until ascitic volume ";
  parameter Real Pbreak =  8.0 "pressure required to break liver lymphatics";
  parameter Real D =  0.8 " Experimental Valuesml/mm Hg   Volume = Vmin+D*(Pa - Po),  in steady state  Vol = D*P";


  //Adjustable paramters:  These are the same values used in Ascites_state_grad9 and in paper
  parameter Real m = 0.8 "liver tissue protein = 80% of plasma";
  parameter Real Lt = 6.25*dayconv*volconv "ml/hour/mm Hg for total conductance from paper";
  parameter Real Lc = 2.0*Lt;
  parameter Real Li = Lc "distribute intestinal resistance equally between capilary and mesothelium";
  parameter Real Ll =   10.3*dayconv*volconv;
  parameter Real Ly0 =  7.86*dayconv*volconv "ml/hour mm Hg  steady state value";


  //parameters for intestinal tissue
  parameter Real Vimin = 100*volconv " = 100 ml = int. volume where Pi = Pa (i.e. 0 presssure gradient";
  parameter Real Vimin2 = 50*volconv " have constant negative pressure for Vi <Vimin2";
  parameter Real Di =  (0.75/100)/volconv "mm/ml Hg  P = Di*Vi  NOTE: units are the inverse of D - this corresponds to 0.75 mm Hg pressure increase for doubling of volume from 100 to 200 ml";
  parameter Real Lyi =  18*dayconv*volconv "int. tissue lymph flow - about twice the conductance of the peritoneal space";
  //Real Vi0 = 110*volconv "initial intestinal volume";
  parameter Real Perm =  2.0*volconv "capillary protein permeability";

  parameter Modelica.Units.SI.Time t_e =  30 "event time";
  parameter Real Pgrad1 = 16, Pgrad2 = 12.8;
  parameter Real Pra1 = 5,Pra2 = 2;

  Real Pgrad = if time < t_e then Pgrad1 else Pgrad2;
  Real Pra = if time < t_e then Pra1 else Pra2;
  Real Ly = if paracentesis and time > t_e and time < t_e + 0.5 then Ly0*50 else Ly0 "At Pgrad = 20.0 and Pra = 5.0;";



  //Initial condtions:
  parameter Real V0 = Vmin "initial ascites volume";
  // Real Pp = Phv+4 "Cannot use Phv+Pgrad1 because this make Ai <0 because of low Pa";
  parameter Real Vi0 =  Vimin "Initial intestine volume follow from Pi = Pa";
  //Real Ai = Ap+Pi-Pc "equilibrium across capillary";
  Real Amti0 = Ai*Vi0 "initial intestine amount - Ai*VI0";
  Real Amt0 = Aa*V0 "initial ascite amount";

  // ODEs
  Real Va(start = V0);
  Real Amt(start = Amt0);
  Real Vi(start = Vi0);
  Real Amti(start = Amti0);


  // Simple parameter relations (algebraic, not requiring diff. eq.
  Real Pp(start = Phv+4) = Phv+Pgrad "portal vein pressure";
  Real Pl = (Pp + Phv)/2.0 "liver tissue pressure = average sinusoidal pressure";

  Real Pa( start = Pmin) = if Va<=Vmin then Po else Po+(Va-Vmin)/D "linear,  see older version for quadratic";
  // Real Phv(start = Pra1+3) = if Pra+Pvg >Pa then Pra+Pvg else Pa "P hepatic vein = r. atrium +3 is this is greater in Pa, otherwise Pa";
  Real Phv(start = Pra1+3) = max(Pra+Pvg, Pa) "P hepatic vein = r. atrium +3 is this is greater in Pa, otherwise Pa";
  Real Pc = Pp+3 "Intestinal capillary pressure";
  Real Pi(start = Pa) = if Vi<=Vimin2 then Pa+Di*(Vimin2-Vimin) else Pa+Di*(Vi-Vimin) "constant negative pressure for Vi<Vimin2, varying negative for Vi<Vimin, positive for Vi&>Vimin ";

  Real Aa(start = Ai) = if Amt<=0 then 0.001 else Amt/Va "ascites protein conc.";
  Real Ai( start = Ap+Pi-Pc) =  if Amti<=0 then 0.001 else Amti/Vi "int. tissue protein conc.";

  Real Ji =   Li*(Aa - Ai+Pi-Pa) "volume flow across intestinal mesotheium into ascites";
  Real Jl =  if (Pl-Pa)<Pbreak then 0 else Ll*(Pl-Pa-Pbreak) "volume flow from liver";
  Real Jy = if Pa<=Po then 0 else Ly*(Pmin+Pa-Pra) "0 if Va<0, Pmin pressure if Pa<Pra,else Pmin+Pa-Pra";

  Real Jla = m*Ap*Jl "albumin flow from liver";
  Real Jya = Aa*Jy " peritoneal lymph protrein removal rate (Aa = ascites protein)";

  //New relations for intestinal tissue compartment
  Real Jc =  Lc*(Pc - Pi+ Ai -Ap) "capillary water flow";
  Real Jyi = if (Pi-Pa)<0 then 0 else Lyi*(Pi-Pa) "For some reason, this integrates rapidly gives reasonable results";
  Real JcA = Perm*(Ap - Ai) "New -add a capillary protein permeabilit to balance int. lymph flow protein";


  equation
  //Differential equations:
  der(Va) = Ji+Jl-Jy;
  der(Amt)=Jla-Jya;
  der(Vi)=Jc-Jyi-Ji;
  der(Amti)= JcA  -Jyi*Ai;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=60,
        __Dymola_NumberOfIntervals=5000,
        Tolerance=1e-05,
        __Dymola_Algorithm="Dassl"));
  end LevitCase2Dynamics;

  package Auxiliary
    model OncoticPressureToConc
      parameter Modelica.Units.SI.Temperature T = 310;
      parameter Real Wp(unit = "g/mol") = 65e3 "Weight of protein, kg/mol";
      Physiolibrary.Types.Pressure Ponc = 0 annotation (Dialog(group="Time varying output signal"));
      Physiolibrary.Types.Concentration c;
      Real conc_gdl(unit = "g/dl") = c*Wp*1e-3/10;

    equation
      Ponc = c*Modelica.Constants.R*T;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end OncoticPressureToConc;

    model POnc2Conc "Calculation acc to Nitta 1981, PMID 7324049, for 37C and pH of 7.4"
      parameter Boolean inputIsPressure=false
                                             "use inputValue as oncotic pressure or as alb conc [g/dl] otherwise"
        annotation (checkbox=true);
        Real inputValue=time "Pressure [Pa] if inputIsPresure, alb conc [g/dl] otherwise" annotation (Dialog(tab="General", group="Inputs"));
            Real alpha = AGf*beta;
            Real beta = 1-alpha;
            parameter Real AGf=4/3   "Albumin to globulin fraction";
      Real c;
      Real c_alb = c*alpha;
      Real c_glb = c*beta;
      Physiolibrary.Types.Pressure Pi "Oncotic pressure";
            parameter Physiolibrary.Types.Pressure mmHg=133.322387415;
    equation
      if inputIsPressure then
        inputValue = Pi;
      else
        inputValue = c_alb;
      end if;

      // This surprisingly does not fit exactly their calculated values. Maybe they also applied conversion for pH and or temp?
      Pi/mmHg = alpha*(2.8*c + 0.18*c^2 + 0.012*c^3) + beta*(0.9^c + 0.12*c^2 + 0.004
        *c^3);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StartTime=1,
          StopTime=10,
          __Dymola_Algorithm="Dassl"));
    end POnc2Conc;
  end Auxiliary;

  model LD_concs
    extends LevitCase2Dynamics(
      Ap=22.0,
      Pgrad2=5.0,
      Pvg=0.0);
    constant Physiolibrary.Types.Pressure mmHg=133.322387415;
    Auxiliary.POnc2Conc pOnc2Conc(
      inputIsPressure=false,
      inputValue=4,
      AGf=4/3)
      annotation (Placement(transformation(extent={{-40,-20},{-20,0}})));
    Auxiliary.POnc2Conc protein_asictes(
      inputIsPressure=true,
      inputValue=Aa,
      AGf=4/3) annotation (Placement(transformation(extent={{0,-20},{20,0}})));
    Auxiliary.POnc2Conc protein_intes(
      inputIsPressure=true,
      inputValue=Ai,
      AGf=4/3) annotation (Placement(transformation(extent={{40,-20},{60,0}})));
  end LD_concs;

  model LSS_concs
    extends LevittCase1SS;
    Auxiliary.POnc2Conc protein_plasma(inputIsPressure=true, inputValue=Ap)
      annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
    Auxiliary.POnc2Conc protein_asc(inputIsPressure=true, inputValue=Aa)
      annotation (Placement(transformation(extent={{-20,20},{0,40}})));
  end LSS_concs;

  package Test
    model Healthy
      extends LD_concs(Vmin=0.01);
    end Healthy;

    model PortalHT
      extends Healthy( Pgrad2=10);
    end PortalHT;

    model Nephrotic
      extends Healthy(Ap=3.6, pOnc2Conc(inputValue=1));
    end Nephrotic;

    model Cirrhosis
      extends Healthy(Pgrad2=20, pOnc2Conc(inputValue=2));
    end Cirrhosis;
  end Test;

  package SSCharacteristics
    model TestPra
      LevittCase1SS levittCase1SS
        annotation (Placement(transformation(extent={{-80,-2},{-60,18}})));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end TestPra;
  end SSCharacteristics;
  annotation (uses(Physiolibrary(version="2.4.1"), Modelica(version="4.0.0")));
end Lymphatics;
