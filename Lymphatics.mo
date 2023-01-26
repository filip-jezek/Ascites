within ;
package Lymphatics

  package AscitesLevitt "Ascites model by Levitt and Levitt"
    package OncoticPressuresConversion
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
    end OncoticPressuresConversion;

    partial model AscitesIO
      Physiolibrary.Types.RealIO.PressureInput Pp "portal vein pressure input"
         annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
      Physiolibrary.Types.RealIO.PressureInput Pra "Right atrial pressure input"
        annotation (Placement(transformation(extent={{-20,-20},{20,20}},
            rotation=180,
            origin={100,0})));
      Physiolibrary.Types.RealIO.PressureOutput Pa "Ascites pressure output"
        annotation (Placement(transformation(extent={{-20,-100},{20,-60}})));
      Physiolibrary.Types.RealIO.PressureInput Phv "Hepatic vein pressure input"
        annotation (Placement(transformation(extent={{-120,-100},{-80,-60}})));
      Physiolibrary.Types.RealIO.PressureInput Pc "Intestinal capillary pressure"
        annotation (Placement(transformation(extent={{-120,60},{-80,100}})));
      annotation (Icon(graphics={Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              lineThickness=0.5)}));
    end AscitesIO;

    model LevittCase1Ss
      "Model of dynamic ascites build-up by Levitt and Levitt (2012), PMID 22453061. Converted from Maple code."
    //   Real D,Ly,Pbreak,Va,R,F,Pra,Pp,Pc,Pa,Pi,Pl,Phv,Aa1,Aa2,Pa1,Pa2,Ai,Jl,Jla,Jy,m,Pmin,Jya,Ap,Pdel,Pgrad,nplot,maxLl,mPa,mLl,i,Ll,x1,x2,Papoint,Lypoint,Setupeqs,delLl,Lt,maxgrad,MaxPa,MaxAa,Pamin,mingrad,delgrad,maxN,Vpoint,Apoint,gradfract,m2Pa,Pgrad0,Vpoint2,delplot;
    //   Real D,Ly,Pbreak,Va,R,F,Pra,Pp,Pc,Pa,Pi,Pl,Phv,Aa1,Aa2,Pa1,Pa2,Jl,Jla,Jy,m,Pmin,Jya,Ap,Pdel,Pgrad,nplot,maxLl,mPa,mLl,i,Ll,x1,x2,Papoint,Lypoint,Setupeqs,delLl,Lt,maxgrad,MaxPa,MaxAa,Pamin,mingrad,delgrad,maxN,Vpoint,gradfract,m2Pa,Pgrad0,Vpoint2,delplot;
    //  Real Ji, Aa, Vmin;

    parameter Real Pra=5   annotation(Evaluate = false);//right atrial pressure
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
    end LevittCase1Ss;

    model LevittCase1SsSI
      "Model of dynamic ascites build-up by Levitt and Levitt (2012), PMID 22453061. Converted from Maple code. With SI units."
    //   Real D,Ly,Pbreak,Va,R,F,Pra,Pp,Pc,Pa,Pi,Pl,Phv,Aa1,Aa2,Pa1,Pa2,Ai,Jl,Jla,Jy,m,Pmin,Jya,Ap,Pdel,Pgrad,nplot,maxLl,mPa,mLl,i,Ll,x1,x2,Papoint,Lypoint,Setupeqs,delLl,Lt,maxgrad,MaxPa,MaxAa,Pamin,mingrad,delgrad,maxN,Vpoint,Apoint,gradfract,m2Pa,Pgrad0,Vpoint2,delplot;
    //   Real D,Ly,Pbreak,Va,R,F,Pra,Pp,Pc,Pa,Pi,Pl,Phv,Aa1,Aa2,Pa1,Pa2,Jl,Jla,Jy,m,Pmin,Jya,Ap,Pdel,Pgrad,nplot,maxLl,mPa,mLl,i,Ll,x1,x2,Papoint,Lypoint,Setupeqs,delLl,Lt,maxgrad,MaxPa,MaxAa,Pamin,mingrad,delgrad,maxN,Vpoint,gradfract,m2Pa,Pgrad0,Vpoint2,delplot;
    //  Real Ji, Aa, Vmin;
    constant Real mmHg = 133.322;

    parameter Physiolibrary.Types.Pressure  Pra(displayUnit="mmHg")=266.64477483 "right atrial pressure"
                                                   annotation(Evaluate = false);
    parameter Physiolibrary.Types.Pressure Pamin=266.64477483 "minimum ascites pressure when Jlymph = 0";
    parameter Physiolibrary.Types.Pressure Ap(displayUnit="mmHg")=3333.059685375 "Blood colloid osmotic pressure";
    parameter Real m =  0.8;
    parameter Physiolibrary.Types.Pressure Pmin(displayUnit="mmHg")=266.64477483 "Minmal ascites pressure. Must be less than Pra";
    parameter Physiolibrary.Types.Pressure Pdel(displayUnit="mmHg")=266.64477483 "Hepatic vein pressure drop. Assumed same as Pmin"
                                                      annotation(Evaluate=false);
    parameter Physiolibrary.Types.Pressure Pbreak(displayUnit="mmHg")=1066.57909932;
    parameter Physiolibrary.Types.HydraulicConductance Ly(displayUnit="ml/(mmHg.min)")=
         1.6376344405963e-11;                                       //ml/min/mm Hg
    parameter Physiolibrary.Types.HydraulicConductance Ll(displayUnit="ml/(mmHg.min)")=
         2.1501765174242e-11;                                       //
    parameter Physiolibrary.Types.HydraulicConductance Lt(displayUnit="ml/(mmHg.min)")=
         1.3001067314658e-11;
    parameter Physiolibrary.Types.HydraulicCompliance D(displayUnit="l/mmHg")=6.0004926067653e-06;
                                                                //New  - value of Henriksen and Lieberman. 0.8 L/mmHg
    parameter Physiolibrary.Types.Volume Vmin(displayUnit="l")=0.0001;
                                                       //minimum volume when P = Pamin

    parameter Physiolibrary.Types.Pressure PcGrad=399.967162245
                                                      "Gradient from portal to intestinal capillary pressure";
    Physiolibrary.Types.Pressure Phv "Hepatic vein pressure";
    Physiolibrary.Types.Pressure Pp "Portal vein pressure";
    Physiolibrary.Types.Pressure Pc "Intestinal capillary pressure";
    Physiolibrary.Types.Pressure Pl "Liver sinus pressure";
    Physiolibrary.Types.VolumeFlowRate Jl;
    Physiolibrary.Types.VolumeFlowRate Ji;
    Physiolibrary.Types.VolumeFlowRate Jy;

    Physiolibrary.Types.Pressure Pa;
    Physiolibrary.Types.Pressure Aa;
    Physiolibrary.Types.Pressure Pgrad = time*mmHg annotation (Dialog(tab="General", group="Inputs"));

    Physiolibrary.Types.Volume Av = Vmin+D*(Pa-Pamin);
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
    Pc = Pp+PcGrad;//intestinal capillary pressure
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
          StartTime=-5,
          StopTime=25,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Dassl"));
    end LevittCase1SsSI;

    model LevittCase1SsSiIo
      "Model of static ascites build-up by Levitt and Levitt (2012), PMID 22453061. Converted from Maple code into SI units with IO ports."
      extends AscitesIO;
    //   Real D,Ly,Pbreak,Va,R,F,Pra,Pp,Pc,Pa,Pi,Pl,Phv,Aa1,Aa2,Pa1,Pa2,Ai,Jl,Jla,Jy,m,Pmin,Jya,Ap,Pdel,Pgrad,nplot,maxLl,mPa,mLl,i,Ll,x1,x2,Papoint,Lypoint,Setupeqs,delLl,Lt,maxgrad,MaxPa,MaxAa,Pamin,mingrad,delgrad,maxN,Vpoint,Apoint,gradfract,m2Pa,Pgrad0,Vpoint2,delplot;
    //   Real D,Ly,Pbreak,Va,R,F,Pra,Pp,Pc,Pa,Pi,Pl,Phv,Aa1,Aa2,Pa1,Pa2,Jl,Jla,Jy,m,Pmin,Jya,Ap,Pdel,Pgrad,nplot,maxLl,mPa,mLl,i,Ll,x1,x2,Papoint,Lypoint,Setupeqs,delLl,Lt,maxgrad,MaxPa,MaxAa,Pamin,mingrad,delgrad,maxN,Vpoint,gradfract,m2Pa,Pgrad0,Vpoint2,delplot;
    //  Real Ji, Aa, Vmin;
    constant Real mmHg = 133.322;

    // parameter Physiolibrary.Types.Pressure  Pra(displayUnit="mmHg")=266.64477483 "right atrial pressure"
    //                                                annotation(Evaluate = false);
    parameter Physiolibrary.Types.Pressure Pamin=266.64477483 "minimum ascites pressure when Jlymph = 0";
    parameter Physiolibrary.Types.Pressure Ap(displayUnit="mmHg")=3333.059685375 "Blood colloid osmotic pressure";
    parameter Real m =  0.8;
    parameter Physiolibrary.Types.Pressure Pmin(displayUnit="mmHg")=266.64477483 "Minmal ascites pressure. Must be less than Pra";
    // parameter Physiolibrary.Types.Pressure Pdel(displayUnit="mmHg")=266.64477483 "Hepatic vein pressure drop. Assumed same as Pmin"
    //                                                   annotation(Evaluate=false);
    parameter Physiolibrary.Types.Pressure Pbreak(displayUnit="mmHg")=1066.57909932;
    parameter Physiolibrary.Types.HydraulicConductance Ly(displayUnit="ml/(mmHg.min)")=
         1.6376344405963e-11;                                       //ml/min/mm Hg
    parameter Physiolibrary.Types.HydraulicConductance Ll(displayUnit="ml/(mmHg.min)")=
         2.1501765174242e-11;                                       //
    parameter Physiolibrary.Types.HydraulicConductance Lt(displayUnit="ml/(mmHg.min)")=
         1.3001067314658e-11;
    parameter Physiolibrary.Types.HydraulicCompliance D(displayUnit="l/mmHg")=6.0004926067653e-06;
                                                                //New  - value of Henriksen and Lieberman. 0.8 L/mmHg
    parameter Physiolibrary.Types.Volume Vmin(displayUnit="l")=0.0001;
                                                       //minimum volume when P = Pamin

    // parameter Physiolibrary.Types.Pressure PcGrad=399.967162245
    //                                                   "Gradient from portal to intestinal capillary pressure";
    // Physiolibrary.Types.Pressure Phv "Hepatic vein pressure";
    // Physiolibrary.Types.Pressure Pp "Portal vein pressure";
    // Physiolibrary.Types.Pressure Pc "Intestinal capillary pressure";
    Physiolibrary.Types.Pressure Pl "Liver sinus pressure";
    Physiolibrary.Types.VolumeFlowRate Jl;
    Physiolibrary.Types.VolumeFlowRate Ji;
    Physiolibrary.Types.VolumeFlowRate Jy;

    // Physiolibrary.Types.Pressure Pa;
    Physiolibrary.Types.Pressure Aa;

    Physiolibrary.Types.Volume Av = Vmin+ max(0, D*(Pa-Pamin));

    equation

        Jy = if (Pmin+Pa-Pra<=0) then 0 else Ly*(Pmin+Pa-Pra);//0 Simple case Pa&gt;Pra
    Pl = (Pp+Phv)/2.0;//Liver sinusoid pressure
    Jl = Ll*(Pl-Pa-Pbreak);//Simple case, assume Pl-Pa &gt; Pbreak
    Ji = Lt*(Pc-Pa-Ap+Aa); //flux from intestine to ascites space

    Aa*Jy=m*Ap*Jl;//protein balance, //Pa =
    Jy=Jl+Ji;//fluid balance, //Aa =

      annotation (experiment(
          StartTime=6,
          StopTime=25,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Dassl"), Icon(graphics={Text(
              extent={{-100,-100},{100,100}},
              textColor={28,108,200},
              textString="Levitt
Steady-state
ascites")}));
    end LevittCase1SsSiIo;

    model LevitCase2Dynamics
      "Model of dynamic ascites build-up by Levitt and Levitt (2012), PMID 22453061. Converted from Maple code."
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
    parameter Real Pc0 = Pra1+(3+4+3);
    parameter Real Amti0 = (Ap+Pmin-Pc0)*Vi0 "initial intestine amount - Ai*VI0";
    parameter Real Amt0 = (Ap+Pmin-Pc0)*V0 "initial ascite amount";

    // ODEs
    Real Va(start = V0);
    Real Amt(start = Amt0);
    Real Vi(start = Vi0);
    Real Amti(start = Amti0);

    // Simple parameter relations (algebraic, not requiring diff. eq.
    Real Pp = Phv+Pgrad "portal vein pressure";
    Real Pl = (Pp + Phv)/2.0 "liver tissue pressure = average sinusoidal pressure";

    Real Pa = Po+max((Va-Vmin)/D,0) "linear,  see older version for quadratic";
    // Real Phv(start = Pra1+3) = if Pra+Pvg >Pa then Pra+Pvg else Pa "P hepatic vein = r. atrium +3 is this is greater in Pa, otherwise Pa";
    Real Phv = max(Pra+Pvg, Pa) "P hepatic vein = r. atrium +3 is this is greater in Pa, otherwise Pa";
    Real Pc = Pp+3 "Intestinal capillary pressure";
    Real Pi = if Vi<=Vimin2 then Pa+Di*(Vimin2-Vimin) else Pa+Di*(Vi-Vimin) "constant negative pressure for Vi<Vimin2, varying negative for Vi<Vimin, positive for Vi&>Vimin ";

    Real Aa = if Amt<=0 then 0.001 else Amt/Va "ascites protein conc.";
    Real Ai =  if Amti<=0 then 0.001 else Amti/Vi "int. tissue protein conc.";

    Real Ji =   Li*(Aa - Ai+Pi-Pa) "volume flow across intestinal mesotheium into ascites";
    Real Jl =  max(Ll*(Pl-Pa-Pbreak), 0) "volume flow from liver";
    Real Jy = if Pa<=Po then 0 else Ly*(Pmin+Pa-Pra) "0 if Va<0, Pmin pressure if Pa<Pra,else Pmin+Pa-Pra";

    Real Jla = m*Ap*Jl "albumin flow from liver";
    Real Jya = Aa*Jy " peritoneal lymph protrein removal rate (Aa = ascites protein)";

    //New relations for intestinal tissue compartment
    Real Jc =  Lc*(Pc - Pi+ Ai -Ap) "capillary water flow";
    Real Jyi = max(Lyi*(Pi-Pa),0) "For some reason, this integrates rapidly gives reasonable results";
    Real JcA = Perm*(Ap - Ai) "New -add a capillary protein permeabilit to balance int. lymph flow protein";

    equation
    //Differential equations:
    der(Va) = Ji+Jl-Jy;
    //Va = 10;
    der(Amt)=Jla-Jya;
    der(Vi)=Jc-Jyi-Ji;

    //Vi = 0.1;
    der(Amti)= JcA  -Jyi*Ai;

    if Va < 1e-6 then
      terminate("Ascites at zero");
    end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=60,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Dassl"));
    end LevitCase2Dynamics;

    model LevitCase2DynamicsSI
      "Model of dynamic ascites build-up by Levitt and Levitt (2012), PMID 22453061. Converted from Maple code."
      import Physiolibrary.Types.*;
      parameter Boolean paracentesis = false;
      parameter Real speedup = 60*60*24;
    //Independent variables
    parameter Real dayconv =  24.0 "hours/day   - thus paramteters in units of 1/day";
    parameter Real volconv = 0.001 "liter/ml";

    //volconv = 1.0 "liter/ml";
    parameter Pressure Ap=3333.059685375
                               "Plasma colloid osmoitc pressure";
    parameter Pressure Po=266.64477483
                                 "minimum pressure in the pressure Pa vs Volume relation";
    parameter Pressure Pvg=266.64477483
                                 "Hep[atic vein to CVP gradient";
    parameter Pressure Pmin = Po "minimum pressure in the lymph flow equation (in paper - assume that = Po";
    parameter Volume Vmin(displayUnit="l")=0.0001
                                 "No lymph flow until ascitic volume ";
    parameter Pressure Pbreak=1066.57909932
                                     "pressure required to break liver lymphatics";
    parameter Physiolibrary.Types.HydraulicCompliance D(displayUnit="l/mmHg")=6.0004926067653e-06
                                                               " Experimental Valuesml/mm Hg   Volume = Vmin+D*(Pa - Po),  in steady state  Vol = D*P";

    //Adjustable paramters:  These are the same values used in Ascites_state_grad9 and in paper
    parameter Real m = 0.8 "liver tissue protein = 80% of plasma";
    parameter Physiolibrary.Types.HydraulicConductance Lt(displayUnit="ml/(mmHg.h)")=
         1.3021902358432e-11                                                     "ml/hour/mm Hg for total conductance from paper";
    parameter HydraulicConductance Lc = 2.0*Lt;
    parameter HydraulicConductance Li = Lc "distribute intestinal resistance equally between capilary and mesothelium";
    parameter HydraulicConductance Ll(displayUnit="ml/(mmHg.h)")=2.1460095086695e-11;
    parameter HydraulicConductance Ly0(displayUnit="ml/(mmHg.h)")=1.6376344405963e-11
                                                               "ml/hour mm Hg  steady state value";

    //parameters for intestinal tissue
    parameter Volume Vimin=0.0001        " = 100 ml = int. volume where Pi = Pa (i.e. 0 presssure gradient";
    parameter Volume Vimin2=5e-05        " have constant negative pressure for Vi <Vimin2";
    parameter Physiolibrary.Types.HydraulicElastance Di=133322387.415*((0.75/100))
                                                                              "mm Hg/ml P = Di*Vi  NOTE: units are the inverse of D - this corresponds to 0.75 mm Hg pressure increase for doubling of volume from 100 to 200 ml";
    parameter HydraulicConductance Lyi(displayUnit="ml/(mmHg.h)")=3.7503078792283e-11
                                                             "int. tissue lymph flow - about twice the conductance of the peritoneal space";
    //Real Vi0 = 110*volconv "initial intestinal volume";
    parameter Physiolibrary.Types.DiffusionPermeability Perm(displayUnit="ml/day")=2.3148148148148e-11
                                                                            "capillary protein permeability ml/d";

    parameter Modelica.Units.SI.Time t_e(displayUnit="s")=30
                                               "event time";
    parameter Pressure Pgrad1=2133.15819864,
                                    Pgrad2=1706.526558912;
    parameter Pressure Pra1=666.611937075,
                                Pra2=266.64477483;

    Pressure Pgrad = if time < t_e then Pgrad1 else Pgrad2;
    Pressure Pra = if time < t_e then Pra1 else Pra2;
    HydraulicConductance Ly = if paracentesis and time > t_e and time < t_e + 0.5 then Ly0*50 else Ly0 "At Pgrad = 20.0 and Pra = 5.0;";

    //Initial condtions:
    parameter Volume V0 = Vmin "initial ascites volume";
    // Real Pp = Phv+4 "Cannot use Phv+Pgrad1 because this make Ai <0 because of low Pa";
    parameter Volume Vi0 =  Vimin "Initial intestine volume follow from Pi = Pa";
    //Real Ai = Ap+Pi-Pc "equilibrium across capillary";
    parameter Pressure Pc0 = Pra1+133.32*(3+4+3);
    parameter Real Amti0 = (Ap+Pmin-Pc0)*Vi0 "initial intestine amount - Ai*VI0";
    parameter Real Amt0 = (Ap+Pmin-Pc0)*V0 "initial ascite amount";

    // ODEs
    Volume Va(start = V0);
    Real Amt(start = Amt0, unit = "Pa.m3", displayUnit = "mmHg.l");
    Real Amt_dbg = Amt/133.322*1000;
    Volume Vi(start = Vi0);
    Real Amti(start = Amti0, unit = "Pa.m3", displayUnit = "mmHg.l");
    Real Amti_dbg = Amti/133.322*1000;

    // Simple parameter relations (algebraic, not requiring diff. eq.
    Pressure Pp = Phv+Pgrad "portal vein pressure";
    Pressure Pl = (Pp + Phv)/2.0 "liver tissue pressure = average sinusoidal pressure";

    Pressure Pa = Po+max((Va-Vmin)/D, 0) "linear,  see older version for quadratic";
    // Real Phv(start = Pra1+3) = if Pra+Pvg >Pa then Pra+Pvg else Pa "P hepatic vein = r. atrium +3 is this is greater in Pa, otherwise Pa";
    Pressure Phv = max(Pra+Pvg, Pa) "P hepatic vein = r. atrium +3 is this is greater in Pa, otherwise Pa";
    Pressure Pc = Pp+3*133.32 "Intestinal capillary pressure";
    Pressure Pi = if Vi<=Vimin2 then Pa+Di*(Vimin2-Vimin) else Pa+Di*(Vi-Vimin) "constant negative pressure for Vi<Vimin2, varying negative for Vi<Vimin, positive for Vi&>Vimin ";

    Pressure Aa = if Amt<=0 then 0.001 else Amt/Va "ascites protein conc.";
    Pressure Ai =  if Amti<=0 then 0.001 else Amti/Vi "int. tissue protein conc.";

    VolumeFlowRate Ji =   Li*(Aa - Ai+Pi-Pa) "volume flow across intestinal mesotheium into ascites";
    VolumeFlowRate Jl =  if (Pl-Pa)<Pbreak then 0 else Ll*(Pl-Pa-Pbreak) "volume flow from liver";
    VolumeFlowRate Jy = if Pa<=Po then 0 else Ly*(Pmin+Pa-Pra) "0 if Va<0, Pmin pressure if Pa<Pra,else Pmin+Pa-Pra";

    Real Jla = m*Ap*Jl "albumin flow from liver";
    Real Jya = Aa*Jy " peritoneal lymph protrein removal rate (Aa = ascites protein) - mmHg*m3/s";

    //New relations for intestinal tissue compartment
    VolumeFlowRate Jc =  Lc*(Pc - Pi+ Ai -Ap) "capillary water flow";
    VolumeFlowRate Jyi = if (Pi-Pa)<0 then 0 else Lyi*(Pi-Pa) "For some reason, this integrates rapidly gives reasonable results";
    Real JcA = Perm*(Ap - Ai) "New -add a capillary protein permeabilit to balance int. lymph flow protein";

    equation
    //Differential equations:
    der(Va)/speedup = Ji+Jl-Jy;
    //Va = 1e-3*10;
    der(Amt)/speedup=Jla-Jya;
    der(Vi)/speedup=Jc-Jyi-Ji;
    //Vi = 0.1*1e-3;
    der(Amti)/speedup= JcA  -Jyi*Ai;

    if Va < 1e-6 then
      terminate("Ascites at zero");
    end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=60,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-09,
          __Dymola_Algorithm="Dassl"));
    end LevitCase2DynamicsSI;

    model LevitCase2DynamicsSiIo
      "Model of dynamic ascites build-up by Levitt and Levitt (2012), PMID 22453061. Converted from Maple code, units transformed to SI and cut into variable inputs"
      import Physiolibrary.Types.*;
      extends AscitesIO;

      parameter Boolean paracentesis = false;
      parameter Real speedup = 60*60*24;
    //Independent variables
    parameter Real dayconv =  24.0 "hours/day   - thus paramteters in units of 1/day";
    parameter Real volconv = 0.001 "liter/ml";

    //volconv = 1.0 "liter/ml";
    parameter Pressure Ap=3333.059685375
                               "Plasma colloid osmoitc pressure";
    parameter Pressure Po=266.64477483
                                 "minimum pressure in the pressure Pa vs Volume relation";
    parameter Pressure Pvg=266.64477483
                                 "Hep[atic vein to CVP gradient";
    parameter Pressure Pmin = Po "minimum pressure in the lymph flow equation (in paper - assume that = Po";
    parameter Volume Vmin(displayUnit="l")=0.0001
                                 "No lymph flow until ascitic volume ";
    parameter Pressure Pbreak=1066.57909932
                                     "pressure required to break liver lymphatics";
    parameter Physiolibrary.Types.HydraulicCompliance D(displayUnit="l/mmHg")=
        6.0004926067653e-06                                    " Experimental Valuesml/mm Hg   Volume = Vmin+D*(Pa - Po),  in steady state  Vol = D*P";

    //Adjustable paramters:  These are the same values used in Ascites_state_grad9 and in paper
    parameter Real m = 0.8 "liver tissue protein = 80% of plasma";
    parameter Physiolibrary.Types.HydraulicConductance Lt(displayUnit=
            "ml/(mmHg.h)")=1.3021902358432e-11                                   "ml/hour/mm Hg for total conductance from paper";
    parameter HydraulicConductance Lc = 2.0*Lt;
    parameter HydraulicConductance Li = Lc "distribute intestinal resistance equally between capilary and mesothelium";
    parameter HydraulicConductance Ll(displayUnit="ml/(mmHg.h)")=
        2.1460095086695e-11;
    parameter HydraulicConductance Ly0(displayUnit="ml/(mmHg.h)")=
        1.6376344405963e-11                                    "ml/hour mm Hg  steady state value";

    //parameters for intestinal tissue
    parameter Volume Vimin=0.0001        " = 100 ml = int. volume where Pi = Pa (i.e. 0 presssure gradient";
    parameter Volume Vimin2=5e-05        " have constant negative pressure for Vi <Vimin2";
    parameter Physiolibrary.Types.HydraulicElastance Di=133322387.415*((0.75/
          100))                                                               "mm Hg/ml P = Di*Vi  NOTE: units are the inverse of D - this corresponds to 0.75 mm Hg pressure increase for doubling of volume from 100 to 200 ml";
    parameter HydraulicConductance Lyi(displayUnit="ml/(mmHg.h)")=
        3.7503078792283e-11                                  "int. tissue lymph flow - about twice the conductance of the peritoneal space";
    //Real Vi0 = 110*volconv "initial intestinal volume";
    parameter Physiolibrary.Types.DiffusionPermeability Perm(displayUnit=
            "ml/day")=2.3148148148148e-11                                   "capillary protein permeability ml/d";

    parameter Modelica.Units.SI.Time t_e(displayUnit="s")=30
                                               "event time";


    HydraulicConductance Ly = if paracentesis and time > t_e and time < t_e + 0.5 then Ly0*50 else Ly0 "At Pgrad = 20.0 and Pra = 5.0;";

    //Initial condtions:
    parameter Volume V0 = Vmin "initial ascites volume";
    // Real Pp = Phv+4 "Cannot use Phv+Pgrad1 because this make Ai <0 because of low Pa";
    parameter Volume Vi0 =  Vimin "Initial intestine volume follow from Pi = Pa";
    //Real Ai = Ap+Pi-Pc "equilibrium across capillary";

    parameter Pressure Pc0=3999.67162245;
    parameter Real Amti0 = (Ap+Pmin-Pc0)*Vi0 "initial intestine amount - Ai*VI0";
    parameter Real Amt0 = (Ap+Pmin-Pc0)*V0 "initial ascite amount";

    // ODEs
    Volume Va(start = V0);
    Real Amt(start = Amt0, unit = "Pa.m3", displayUnit = "mmHg.l");
    Real Amt_dbg = Amt/133.322*1000;
    Volume Vi(start = Vi0);
    Real Amti(start = Amti0, unit = "Pa.m3", displayUnit = "mmHg.l");
    Real Amti_dbg = Amti/133.322*1000;

    // Simple parameter relations (algebraic, not requiring diff. eq.
    Pressure Pl = (Pp + Phv)/2.0 "liver tissue pressure = average sinusoidal pressure";

    Pressure Pi = if Vi<=Vimin2 then Pa+Di*(Vimin2-Vimin) else Pa+Di*(Vi-Vimin) "constant negative pressure for Vi<Vimin2, varying negative for Vi<Vimin, positive for Vi&>Vimin ";

    Pressure Aa = if Amt<=0 then 0.001 else Amt/Va "ascites protein conc.";
    Pressure Ai =  if Amti<=0 then 0.001 else Amti/Vi "int. tissue protein conc.";

    VolumeFlowRate Ji =   Li*(Aa - Ai+Pi-Pa) "volume flow across intestinal mesotheium into ascites";
    VolumeFlowRate Jl =  if (Pl-Pa)<Pbreak then 0 else Ll*(Pl-Pa-Pbreak) "volume flow from liver";
    VolumeFlowRate Jy = if Pa<=Po then 0 else Ly*(Pmin+Pa-Pra) "0 if Va<0, Pmin pressure if Pa<Pra,else Pmin+Pa-Pra";

    Real Jla = m*Ap*Jl "albumin flow from liver";
    Real Jya = Aa*Jy " peritoneal lymph protrein removal rate (Aa = ascites protein) - mmHg*m3/s";

    //New relations for intestinal tissue compartment
    VolumeFlowRate Jc =  Lc*(Pc - Pi+ Ai -Ap) "capillary water flow";
    VolumeFlowRate Jyi = if (Pi-Pa)<0 then 0 else Lyi*(Pi-Pa) "For some reason, this integrates rapidly gives reasonable results";
    Real JcA = Perm*(Ap - Ai) "New -add a capillary protein permeabilit to balance int. lymph flow protein";

    equation
    Pa = Po+max((Va-Vmin)/D, 0) "linear,  see older version for quadratic";


    //Differential equations:
    der(Va)/speedup = Ji+Jl-Jy;
    //Va = 1e-3*10;
    der(Amt)/speedup=Jla-Jya;
    der(Vi)/speedup=Jc-Jyi-Ji;
    //Vi = 0.1*1e-3;
    der(Amti)/speedup= JcA  -Jyi*Ai;

    if Va < 1e-6 then
      terminate("Ascites at zero");
    end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Text(
              extent={{-100,-100},{100,100}},
              textColor={28,108,200},
              textString="Levitt
transient
ascites")}),                                                         Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=60,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-09,
          __Dymola_Algorithm="Dassl"));
    end LevitCase2DynamicsSiIo;

    model LD_concs
      import Lymphatics;
      import Lymphatics;
      import Lymphatics;
      extends AscitesLevitt.LevitCase2Dynamics(
        Ap=22.0,
        Pgrad2=5.0,
        Pvg=0.0);
      constant Physiolibrary.Types.Pressure mmHg=133.322387415;
      Lymphatics.AscitesLevitt.OncoticPressuresConversion.POnc2Conc pOnc2Conc(
        inputIsPressure=false,
        inputValue=4,
        AGf=4/3)
        annotation (Placement(transformation(extent={{-40,-20},{-20,0}})));
      Lymphatics.AscitesLevitt.OncoticPressuresConversion.POnc2Conc
        protein_asictes(
        inputIsPressure=true,
        inputValue=Aa,
        AGf=4/3) annotation (Placement(transformation(extent={{0,-20},{20,0}})));
      Lymphatics.AscitesLevitt.OncoticPressuresConversion.POnc2Conc
        protein_intes(
        inputIsPressure=true,
        inputValue=Ai,
        AGf=4/3)
        annotation (Placement(transformation(extent={{40,-20},{60,0}})));
    end LD_concs;

    model LSS_concs
      import Lymphatics;
      import Lymphatics;
      extends Lymphatics.AscitesLevitt.LevittCase1Ss;
      Lymphatics.AscitesLevitt.OncoticPressuresConversion.POnc2Conc
        protein_plasma(inputIsPressure=true, inputValue=Ap)
        annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
      Lymphatics.AscitesLevitt.OncoticPressuresConversion.POnc2Conc protein_asc(
          inputIsPressure=true, inputValue=Aa)
        annotation (Placement(transformation(extent={{-20,20},{0,40}})));
    end LSS_concs;

    package Test
      extends Modelica.Icons.ExamplesPackage;
      model Healthy
        extends Modelica.Icons.Example;
        extends AscitesLevitt.LD_concs(
                         Vmin=0.01);
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

  end AscitesLevitt;

  package Hemodynamics
    package WholeBodyCirculation
      model Kofranek2014
        extends
          Physiolibrary.Hydraulic.Examples.Kofranek2014.NonPulsatileCirculation(
            RT(k(displayUnit="(mmHg.min)/l") = Modelica.Constants.inf));
        Physiolibrary.Hydraulic.Components.Conductor
                 nonMuscle(Conductance(displayUnit="m3/(Pa.s)") = 2.6e-9)
          annotation (Placement(transformation(extent={{-4,-112},{16,-92}})));
        Physiolibrary.Hydraulic.Components.Conductor
                 muscle(Conductance(displayUnit="l/(mmHg.min)")=
            1.3001067314658e-09)
          annotation (Placement(transformation(extent={{-4,-94},{16,-74}})));
        Components.Ascites_Resistance ascites_Resistance(
          liverConductance(y=1/(max(time/200, 5)*ascites_Resistance.mmHg/
                ascites_Resistance.Qnom)),
          TIPSS(enable=true, Resistance=Modelica.Constants.inf),
          Liver(enable=true, useConductanceInput=true))
          annotation (Placement(transformation(extent={{16,-154},{-4,-134}})));
        ADAN_main.Components.Subsystems.Systemic.Organs.Renal.Renal_P_Int_i
          renal_P_Int_i(halving=1)
          annotation (Placement(transformation(extent={{16,-130},{-4,-126}})));
        ADAN_main.Components.Interfaces.LeveledPressureFlowConverter
          leveledPressureFlowConverter
          annotation (Placement(transformation(extent={{-8,-122},{0,-114}})));
        ADAN_main.Components.Interfaces.LeveledPressureFlowConverter
          leveledPressureFlowConverter1
          annotation (Placement(transformation(extent={{22,-122},{14,-114}})));
      equation
        connect(muscle.q_in, SystemicVeins.q_in) annotation (Line(
            points={{-4,-84},{-10,-84},{-10,-60},{-36,-60}},
            color={0,0,0},
            thickness=1));
        connect(muscle.q_out, TotalSystemicResistance.q_in) annotation (Line(
            points={{16,-84},{22,-84},{22,-60},{16,-60}},
            color={0,0,0},
            thickness=1));
        connect(nonMuscle.q_in, SystemicVeins.q_in) annotation (Line(
            points={{-4,-102},{-10,-102},{-10,-60},{-36,-60}},
            color={0,0,0},
            thickness=1));
        connect(nonMuscle.q_out, TotalSystemicResistance.q_in) annotation (Line(
            points={{16,-102},{22,-102},{22,-60},{16,-60}},
            color={0,0,0},
            thickness=1));
        connect(leveledPressureFlowConverter1.port_a, nonMuscle.q_out)
          annotation (Line(
            points={{22,-118},{22,-102},{16,-102}},
            color={0,0,0},
            thickness=1));
        connect(leveledPressureFlowConverter.port_a, SystemicVeins.q_in)
          annotation (Line(
            points={{-8,-118},{-10,-118},{-10,-60},{-36,-60}},
            color={0,0,0},
            thickness=1));
        connect(renal_P_Int_i.port_a, leveledPressureFlowConverter1.leveledPort_b)
          annotation (Line(
            points={{16,-128},{16,-124},{14,-124},{14,-118}},
            color={162,29,33},
            thickness=0.5));
        connect(renal_P_Int_i.port_b, leveledPressureFlowConverter.leveledPort_b)
          annotation (Line(
            points={{-4,-128},{-4,-124},{0,-124},{0,-118}},
            color={162,29,33},
            thickness=0.5));
        connect(ascites_Resistance.q_in, leveledPressureFlowConverter1.port_a)
          annotation (Line(
            points={{16,-144},{22,-144},{22,-118}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Resistance.q_out, leveledPressureFlowConverter.port_a)
          annotation (Line(
            points={{-4,-144},{-8,-144},{-8,-118}},
            color={0,0,0},
            thickness=1));
      end Kofranek2014;

      model CVS_GCG
        extends Physiolibrary.Hydraulic.Examples.CardiovascularSystem_GCG(
            nonMuscle(Conductance(displayUnit="m3/(Pa.s)") = 2.6e-9));
        Components.Ascites_Resistance ascites_Resistance(
          liverConductance(y=1/(max(time/100, 5)*ascites_Resistance.mmHg/
                ascites_Resistance.Qnom)),
          TIPSS(enable=true, Resistance=Modelica.Constants.inf),
          Liver(enable=true, useConductanceInput=true))
          annotation (Placement(transformation(extent={{-4,-86},{-24,-66}})));
      equation
        connect(ascites_Resistance.q_in, arteries.q_in) annotation (Line(
            points={{-4,-76},{10,-76},{10,-36},{24,-36}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Resistance.q_out, nonMuscle.q_in) annotation (Line(
            points={{-24,-76},{-34,-76},{-34,-36},{-24,-36}},
            color={0,0,0},
            thickness=1));
        annotation (experiment(StopTime=2500, __Dymola_Algorithm="Dassl"));
      end CVS_GCG;

      model CVS_GCG_TIPS
        extends WholeBodyCirculation.CVS_GCG(
                        ascites_Resistance(TIPSS(Resistance=7999343.2449*(15/1.5))));
        annotation (experiment(StopTime=3000, __Dymola_Algorithm="Dassl"));
      end CVS_GCG_TIPS;
    end WholeBodyCirculation;

    package Components
      model partialAscites
        "Levitts ascites with nominal hemodynamics break down"
        replaceable
        AscitesLevitt.LevittCase1SsSiIo levittCase1SsSiIo constrainedby
          AscitesLevitt.AscitesIO
          annotation (Placement(transformation(extent={{42,14},{62,34}})));
        Physiolibrary.Hydraulic.Components.PumpPressureHead
                          HV(useExternalCollapsingPressure=true, p_head0=-266.64477483)
          "Hepatic vein (free)"
          annotation (Placement(transformation(extent={{52,-10},{72,10}})));
        Physiolibrary.Hydraulic.Sensors.PressureMeasure PRa
          "Right atrial pressure"
          annotation (Placement(transformation(extent={{86,18},{66,38}})));
        replaceable
        Physiolibrary.Hydraulic.Components.PumpPressureHead Liver(p_head0=-799.93432449)
          constrainedby Physiolibrary.Hydraulic.Interfaces.OnePort
          annotation (Placement(transformation(extent={{0,-10},{20,10}})));
        replaceable
        Physiolibrary.Hydraulic.Components.PumpPressureHead IntestineVenule(p_head0=-399.967162245)
          constrainedby Physiolibrary.Hydraulic.Interfaces.OnePort
          "Intestinal venule pressure drop"
          annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
        Physiolibrary.Hydraulic.Sensors.PressureMeasure Ppv
          annotation (Placement(transformation(extent={{-20,18},{0,38}})));
        replaceable
        Physiolibrary.Hydraulic.Components.PumpPressureHead IntestinesArt(p_head0=-2666.4477483)
          constrainedby Physiolibrary.Hydraulic.Interfaces.OnePort
          "Intestinal Arteriole resistance"
          annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
        Physiolibrary.Hydraulic.Sensors.PressureMeasure Pc
          annotation (Placement(transformation(extent={{-48,26},{-28,46}})));
        Physiolibrary.Hydraulic.Sensors.PressureMeasure Phv
          annotation (Placement(transformation(extent={{18,10},{38,30}})));
      equation
        connect(PRa.pressure, levittCase1SsSiIo.Pra) annotation (Line(points={{70,24},
                {62,24}},                 color={0,0,127}));
        connect(Liver.q_out, HV.q_in) annotation (Line(
            points={{20,0},{52,0}},
            color={0,0,0},
            thickness=1));
        connect(IntestineVenule.q_out, Liver.q_in) annotation (Line(
            points={{-20,0},{0,0}},
            color={0,0,0},
            thickness=1));
        connect(Ppv.q_in, IntestineVenule.q_out) annotation (Line(
            points={{-14,22},{-14,2.22045e-16},{-20,2.22045e-16}},
            color={0,0,0},
            thickness=1));
        connect(levittCase1SsSiIo.Pp, Ppv.pressure)
          annotation (Line(points={{42,24},{-4,24}},  color={0,0,127}));
        connect(IntestinesArt.q_out, IntestineVenule.q_in) annotation (Line(
            points={{-60,0},{-40,0}},
            color={0,0,0},
            thickness=1));
        connect(IntestinesArt.q_out, Pc.q_in) annotation (Line(
            points={{-60,2.22045e-16},{-42,2.22045e-16},{-42,30}},
            color={0,0,0},
            thickness=1));
        connect(Liver.q_out,Phv. q_in) annotation (Line(
            points={{20,2.22045e-16},{24,2.22045e-16},{24,14}},
            color={0,0,0},
            thickness=1));
        connect(Phv.pressure, levittCase1SsSiIo.Phv)
          annotation (Line(points={{34,16},{42,16}}, color={0,0,127}));
        connect(Pc.pressure, levittCase1SsSiIo.Pc)
          annotation (Line(points={{-32,32},{42,32}}, color={0,0,127}));
        connect(HV.EP, levittCase1SsSiIo.Pa) annotation (Line(points={{52,8},{
                52,16}},                 color={0,0,127}));
        connect(PRa.q_in, HV.q_out) annotation (Line(
            points={{80,22},{80,2.22045e-16},{72,2.22045e-16}},
            color={0,0,0},
            thickness=1));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          experiment(
            StartTime=6,
            StopTime=25,
            __Dymola_Algorithm="Dassl"));
      end partialAscites;

      model Ascites_Resistance "Levitt's model with resistances instead of fied pressure differences"
        extends partialAscites(
          redeclare replaceable Physiolibrary.Hydraulic.Components.Resistor Liver(
              useConductanceInput=true, Resistance=6*mmHg/Qnom),
          redeclare replaceable Physiolibrary.Hydraulic.Components.Resistor IntestinesArt(
              Resistance=84*mmHg/Qnom),
          redeclare replaceable Physiolibrary.Hydraulic.Components.Resistor IntestineVenule(
              Resistance=3*mmHg/Qnom));
        Modelica.Blocks.Sources.RealExpression liverConductance(y=1/((time + 4)
              *mmHg/Qnom))
          annotation (Placement(transformation(extent={{-40,4},{-20,24}})));

        parameter Physiolibrary.Types.VolumeFlowRate Qnom=1.666666666666667e-08*(5000*
            0.2) "Nominal flow through the splanchnic circulation";
        parameter Physiolibrary.Types.Pressure mmHg=133.322;
        replaceable
        Physiolibrary.Hydraulic.Components.Resistor TIPSS(useConductanceInput=false,
            Resistance(displayUnit="(mmHg.min)/l") = 7999343.2449*(15/1.5))
          if useTIPPS
          constrainedby Physiolibrary.Hydraulic.Interfaces.OnePort
          annotation (Placement(transformation(extent={{0,-44},{20,-24}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a
                             q_in "Volume inflow" annotation (Placement(
              transformation(extent={{-114,-14},{-86,14}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_b
                             q_out "Volume outflow"
                               annotation (Placement(
              transformation(extent={{86,-14},{114,14}})));
        parameter Boolean useTIPPS=false;
        Modelica.Blocks.Sources.RealExpression tipssConductance(y=0)
          annotation (Placement(transformation(extent={{-40,-34},{-20,-14}})));

        output Physiolibrary.Types.Pressure PPV = Liver.q_in.pressure "Portal vein pressure";
        output Physiolibrary.Types.Pressure HVPG = Liver.dp "Hepatic venous pressure gradient";
        output Physiolibrary.Types.Volume V_asc = levittCase1SsSiIo.Av "Ascites volume";
        output Physiolibrary.Types.Pressure P_abd = levittCase1SsSiIo.Pa "Abdominal (ascites) pressure";
        output Physiolibrary.Types.VolumeFlowRate Q_liver = Liver.q_in.q;

      equation
        connect(Liver.cond, liverConductance.y)
          annotation (Line(points={{10,6},{10,14},{-19,14}}, color={0,0,127}));
        connect(TIPSS.q_in, IntestineVenule.q_out) annotation (Line(
            points={{0,-34},{-14,-34},{-14,0},{-20,0}},
            color={0,0,0},
            thickness=1));
        connect(TIPSS.q_out, HV.q_in) annotation (Line(
            points={{20,-34},{24,-34},{24,0},{52,0}},
            color={0,0,0},
            thickness=1));
        connect(HV.q_out, q_out) annotation (Line(
            points={{72,0},{100,0}},
            color={0,0,0},
            thickness=1));
        connect(q_in, IntestinesArt.q_in) annotation (Line(
            points={{-100,0},{-80,0}},
            color={0,0,0},
            thickness=1));
        connect(tipssConductance.y, TIPSS.cond) annotation (Line(points={{-19,
                -24},{10,-24},{10,-28}}, color={0,0,127}));
        annotation (Documentation(info="<html>
<p>The TIPS resistance taken from TIPS flow and PVP from Su et al (2012, PMID 22099870).</p>
</html>"));
      end Ascites_Resistance;

      model Ascites_Resistance_Shunts
        extends Ascites_Resistance(levittCase1SsSiIo(D=9.0007389101479e-06));
        ResistancePressureDep splenorenalShunt(
          R_nom(displayUnit="(mmHg.min)/ml") = 7999343244.9,
          Comp(displayUnit="ml/mmHg") = 7.5006157584566e-09,
          P_nom=1199.901486735,
          side=Lymphatics.Hemodynamics.Components.SideEnum.Central,
          useExternalCollapsingPressure=false)
          annotation (Placement(transformation(extent={{0,-72},{20,-52}})));
        output Physiolibrary.Types.VolumeFlowRate Q_shunt=splenorenalShunt.q_in.q;

        Integer phase(start=1);
        parameter Physiolibrary.Types.Fraction collateralFlowFr=0.1    "Fraction of collateral flow to portal vein flow, marking the existence of varices";
        parameter Physiolibrary.Types.Volume V_asc_min(displayUnit="l") = 0.005
          "Minimal observable ascites, marking the existence of ascites";
        parameter Physiolibrary.Types.Pressure HVPG_bleed=2266.480586055
                                                          "Minimal HVPG which results in variceal bleeding";
        replaceable Physiolibrary.Hydraulic.Components.Resistor EsophagealVeins(
            useConductanceInput=true, Resistance=84*mmHg/Qnom)
          "Intestinal Arteriole resistance"
          annotation (Placement(transformation(extent={{20,-100},{40,-80}})));
        replaceable Physiolibrary.Hydraulic.Components.Resistor LeftGastricVein(
            useConductanceInput=true, Resistance=84*mmHg/Qnom)
          "Intestinal Arteriole resistance"
          annotation (Placement(transformation(extent={{-42,-90},{-22,-70}})));
        replaceable Physiolibrary.Hydraulic.Components.Resistor AzygousVein(
            Resistance(displayUnit="(mmHg.min)/l") = 7999343.2449)
          "Intestinal Arteriole resistance"
          annotation (Placement(transformation(extent={{52,-100},{72,-80}})));
        HagenPoiseulleConductance LGVhagenPoiseulleConductance(d_nominal=0.004,
            mu=6e-3)
          annotation (Placement(transformation(extent={{-54,-66},{-34,-46}})));
        HagenPoiseulleConductance LGVhagenPoiseulleConductance1(d_nominal=0.002,
            mu=6e-3)
          annotation (Placement(transformation(extent={{64,-80},{44,-60}})));
        replaceable Physiolibrary.Hydraulic.Components.Resistor gastricArt(
            Resistance=84*mmHg/Qnom) "Intestinal Arteriole resistance"
          annotation (Placement(transformation(extent={{-80,-90},{-60,-70}})));
      equation

        when not Liver.dp > HVPG_bleed and splenorenalShunt.q_in.q >
            collateralFlowFr*IntestinesArt.q_in.q then
          phase = 2;
        elsewhen levittCase1SsSiIo.Av > V_asc_min then
          phase = 3;
        elsewhen Liver.dp > HVPG_bleed and splenorenalShunt.q_in.q >
            collateralFlowFr*IntestinesArt.q_in.q then
          phase = 4;
        end when;


        connect(splenorenalShunt.q_in, IntestineVenule.q_out) annotation (Line(
            points={{0,-62},{-14,-62},{-14,0},{-20,0}},
            color={0,0,0},
            thickness=1));
        connect(splenorenalShunt.q_out, HV.q_in) annotation (Line(
            points={{20,-62},{24,-62},{24,0},{28,0},{28,2.22045e-16},{52,
                2.22045e-16}},
            color={0,0,0},
            thickness=1));
        connect(splenorenalShunt.P_ext, levittCase1SsSiIo.Pa) annotation (Line(
              points={{0,-53},{52,-53},{52,16}}, color={0,0,127}));
        connect(EsophagealVeins.q_out, AzygousVein.q_in) annotation (Line(
            points={{40,-90},{52,-90}},
            color={0,0,0},
            thickness=1));
        connect(AzygousVein.q_out, q_out) annotation (Line(
            points={{72,-90},{78,-90},{78,0},{100,0}},
            color={0,0,0},
            thickness=1));
        connect(LGVhagenPoiseulleConductance.hydraulicconductance,
          LeftGastricVein.cond) annotation (Line(points={{-33.8,-56},{-28,-56},
                {-28,-74},{-32,-74}}, color={0,0,127}));
        connect(LGVhagenPoiseulleConductance1.hydraulicconductance,
          EsophagealVeins.cond) annotation (Line(points={{43.8,-70},{30,-70},{
                30,-84}}, color={0,0,127}));
        connect(gastricArt.q_out, LeftGastricVein.q_in) annotation (Line(
            points={{-60,-80},{-42,-80}},
            color={0,0,0},
            thickness=1));
        connect(LeftGastricVein.q_out, IntestineVenule.q_out) annotation (Line(
            points={{-22,-80},{-14,-80},{-14,0},{-20,0}},
            color={0,0,0},
            thickness=1));
        connect(EsophagealVeins.q_in, LeftGastricVein.q_in) annotation (Line(
            points={{20,-90},{-52,-90},{-52,-80},{-42,-80}},
            color={0,0,0},
            thickness=1));
      end Ascites_Resistance_Shunts;

      model Ascites_Resistance_ShuntsWEso
        extends Ascites_Resistance(levittCase1SsSiIo(D=9.0007389101479e-06));
        ResistancePressureDep splenorenalShunt(
          R_nom(displayUnit="(mmHg.min)/ml") = 7999343244.9,
          Comp(displayUnit="ml/mmHg") = 7.5006157584566e-09,
          P_nom=1199.901486735,
          side=Lymphatics.Hemodynamics.Components.SideEnum.Central,
          useExternalCollapsingPressure=false)
          annotation (Placement(transformation(extent={{0,-72},{20,-52}})));
        output Physiolibrary.Types.VolumeFlowRate Q_shunt=splenorenalShunt.q_in.q;

        Integer phase(start=1);
        parameter Physiolibrary.Types.Fraction collateralFlowFr=0.1    "Fraction of collateral flow to portal vein flow, marking the existence of varices";
        parameter Physiolibrary.Types.Volume V_asc_min(displayUnit="l") = 0.005
          "Minimal observable ascites, marking the existence of ascites";
        parameter Physiolibrary.Types.Pressure HVPG_bleed=2266.480586055
                                                          "Minimal HVPG which results in variceal bleeding";
        replaceable ResistancePressureDep EsophagealVeins(
          useConductanceInput=true,
          L=0.05,
          d(displayUnit="mm") = 0.002,
          P_nom=1333.22387415,
          UsePrescribedDiameter=true) constrainedby
          Physiolibrary.Hydraulic.Components.Resistor(useConductanceInput=true)
          "Intestinal Arteriole resistance"
          annotation (Placement(transformation(extent={{20,-100},{40,-80}})));
        replaceable Physiolibrary.Hydraulic.Components.Resistor LeftGastricVein(
            useConductanceInput=true, Resistance=84*mmHg/Qnom)
          "Intestinal Arteriole resistance"
          annotation (Placement(transformation(extent={{-42,-90},{-22,-70}})));
        replaceable Physiolibrary.Hydraulic.Components.Resistor AzygousVein(
            Resistance(displayUnit="(mmHg.min)/l") = 7999343.2449)
          "Intestinal Arteriole resistance"
          annotation (Placement(transformation(extent={{52,-100},{72,-80}})));
        HagenPoiseulleConductance LGVhagenPoiseulleConductance(d_nominal=0.004,
            mu=6e-3)
          annotation (Placement(transformation(extent={{-54,-66},{-34,-46}})));
        HagenPoiseulleConductance LGVhagenPoiseulleConductance1(d_nominal=0.002,
            mu=6e-3)
          annotation (Placement(transformation(extent={{48,-80},{28,-60}})));
        replaceable Physiolibrary.Hydraulic.Components.Resistor gastricArt(
            Resistance=84*mmHg/Qnom) "Intestinal Arteriole resistance"
          annotation (Placement(transformation(extent={{-80,-90},{-60,-70}})));
        Physiolibrary.Hydraulic.Components.PumpPressureHead
                          HV1(useExternalCollapsingPressure=true, p_head0=0)
          "Hepatic vein (free)"
          annotation (Placement(transformation(extent={{-6,-100},{14,-80}})));
      equation

        when not Liver.dp > HVPG_bleed and splenorenalShunt.q_in.q >
            collateralFlowFr*IntestinesArt.q_in.q then
          phase = 2;
        elsewhen levittCase1SsSiIo.Av > V_asc_min then
          phase = 3;
        elsewhen Liver.dp > HVPG_bleed and splenorenalShunt.q_in.q >
            collateralFlowFr*IntestinesArt.q_in.q then
          phase = 4;
        end when;

        connect(splenorenalShunt.q_in, IntestineVenule.q_out) annotation (Line(
            points={{0,-62},{-14,-62},{-14,0},{-20,0}},
            color={0,0,0},
            thickness=1));
        connect(splenorenalShunt.q_out, HV.q_in) annotation (Line(
            points={{20,-62},{24,-62},{24,0},{28,0},{28,2.22045e-16},{52,
                2.22045e-16}},
            color={0,0,0},
            thickness=1));
        connect(splenorenalShunt.P_ext, levittCase1SsSiIo.Pa) annotation (Line(
              points={{0,-53},{0,-48},{52,-48},{52,16}},
                                                 color={0,0,127}));
        connect(EsophagealVeins.q_out, AzygousVein.q_in) annotation (Line(
            points={{40,-90},{52,-90}},
            color={0,0,0},
            thickness=1));
        connect(AzygousVein.q_out, q_out) annotation (Line(
            points={{72,-90},{78,-90},{78,0},{100,0}},
            color={0,0,0},
            thickness=1));
        connect(LGVhagenPoiseulleConductance.hydraulicconductance,
          LeftGastricVein.cond) annotation (Line(points={{-33.8,-56},{-28,-56},
                {-28,-74},{-32,-74}}, color={0,0,127}));
        connect(gastricArt.q_out, LeftGastricVein.q_in) annotation (Line(
            points={{-60,-80},{-42,-80}},
            color={0,0,0},
            thickness=1));
        connect(LeftGastricVein.q_out, IntestineVenule.q_out) annotation (Line(
            points={{-22,-80},{-14,-80},{-14,0},{-20,0}},
            color={0,0,0},
            thickness=1));
        connect(EsophagealVeins.q_in, HV1.q_out) annotation (Line(
            points={{20,-90},{14,-90}},
            color={0,0,0},
            thickness=1));
        connect(HV1.q_in, LeftGastricVein.q_in) annotation (Line(
            points={{-6,-90},{-52,-90},{-52,-80},{-42,-80}},
            color={0,0,0},
            thickness=1));
        connect(HV1.EP, levittCase1SsSiIo.Pa) annotation (Line(points={{-6,-82},
                {-6,-48},{52,-48},{52,16}},            color={0,0,127}));
      end Ascites_Resistance_ShuntsWEso;

      model Ascites_Resistance_EsoShunt "With esophageal shunt"
        extends Ascites_Resistance(levittCase1SsSiIo(D=9.0007389101479e-06));
        output Physiolibrary.Types.VolumeFlowRate Q_shunt=shunt.q_in.q;

        Integer phase(start=1);
        parameter Physiolibrary.Types.Fraction collateralFlowFr=0.1    "Fraction of collateral flow to portal vein flow, marking the existence of varices";
        parameter Physiolibrary.Types.Volume V_asc_min(displayUnit="l") = 0.005
          "Minimal observable ascites, marking the existence of ascites";
        parameter Physiolibrary.Types.Pressure HVPG_bleed=2266.480586055
                                                          "Minimal HVPG which results in variceal bleeding";
        ResistancePressureDep esophageal_abdominal(
          R_nom(displayUnit="(mmHg.min)/ml") = 799934324490,
          Comp(displayUnit="ml/mmHg") = 7.5006157584566e-09,
          P_nom=0,
          side=Lymphatics.Hemodynamics.Components.SideEnum.Central,
          useExternalCollapsingPressure=true)
          annotation (Placement(transformation(extent={{0,-98},{20,-78}})));
        Physiolibrary.Hydraulic.Components.PumpPressureHead
                          HV1(useExternalCollapsingPressure=true, p_head0=-266.64477483)
          "Hepatic vein (free)"
          annotation (Placement(transformation(extent={{40,-98},{60,-78}})));
        ResistancePressureDep esophageal_thoracic(
          R_nom(displayUnit="(mmHg.min)/ml") = 7999343244900,
          Comp(displayUnit="ml/mmHg") = 7.5006157584566e-10,
          P_nom=1333.22387415,
          side=Lymphatics.Hemodynamics.Components.SideEnum.Central,
          useExternalCollapsingPressure=false)
          annotation (Placement(transformation(extent={{72,-98},{92,-78}})));
      equation

        when shunt.q_in.q > collateralFlowFr*IntestinesArt.q_in.q then
          phase = 2;
        elsewhen levittCase1SsSiIo.Av > V_asc_min then
          phase = 3;
        elsewhen
          Liver.dp > HVPG_bleed then
          phase = 4;
        end when;

        connect(esophageal_abdominal.q_in, IntestineVenule.q_out) annotation (Line(
            points={{0,-88},{-14,-88},{-14,0},{-20,0}},
            color={0,0,0},
            thickness=1));
        connect(esophageal_abdominal.q_out, HV1.q_in) annotation (Line(
            points={{20,-88},{40,-88}},
            color={0,0,0},
            thickness=1));
        connect(HV.EP, esophageal_abdominal.P_ext) annotation (Line(points={{52,8},{52,
                -72},{32,-72},{32,-79},{0,-79}}, color={0,0,127}));
        connect(HV1.EP, esophageal_abdominal.P_ext) annotation (Line(points={{40,-80},
                {32,-80},{32,-79},{0,-79}}, color={0,0,127}));
        connect(esophageal_thoracic.q_in, HV1.q_out) annotation (Line(
            points={{72,-88},{60,-88}},
            color={0,0,0},
            thickness=1));
        connect(esophageal_thoracic.q_out, q_out) annotation (Line(
            points={{92,-88},{96,-88},{96,0},{100,0}},
            color={0,0,0},
            thickness=1));
      end Ascites_Resistance_EsoShunt;

      model Ascites_ResistanceDynamic
        "Levitt's model with resistances instead of fied pressure differences"
        extends partialAscites(
          redeclare replaceable Physiolibrary.Hydraulic.Components.Resistor Liver(
              useConductanceInput=true, Resistance=6*mmHg/Qnom),
          redeclare replaceable Physiolibrary.Hydraulic.Components.Resistor IntestinesArt(
              Resistance=84*mmHg/Qnom),
          redeclare replaceable Physiolibrary.Hydraulic.Components.Resistor IntestineVenule(
              Resistance=3*mmHg/Qnom),
          redeclare AscitesLevitt.LevitCase2DynamicsSiIo levittCase1SsSiIo(speedup=1));
        Modelica.Blocks.Sources.RealExpression liverConductance(y=1/((15)*mmHg/Qnom))
          annotation (Placement(transformation(extent={{-40,4},{-20,24}})));

        parameter Physiolibrary.Types.VolumeFlowRate Qnom=1.666666666666667e-08*(5000*
            0.2) "Nominal flow through the splanchnic circulation";
        parameter Physiolibrary.Types.Pressure mmHg=133.322;
        replaceable
        Physiolibrary.Hydraulic.Components.Resistor TIPSS(useConductanceInput=false,
            Resistance(displayUnit="(mmHg.min)/l") = 7999343.2449*(15/1.5))
          if useTIPPS
          constrainedby Physiolibrary.Hydraulic.Interfaces.OnePort
          annotation (Placement(transformation(extent={{0,-44},{20,-24}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a
                             q_in "Volume inflow" annotation (Placement(
              transformation(extent={{-114,-14},{-86,14}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_b
                             q_out "Volume outflow"
                               annotation (Placement(
              transformation(extent={{86,-14},{114,14}})));
        parameter Boolean useTIPPS=false;
        Modelica.Blocks.Sources.RealExpression tipssConductance(y=0)
          annotation (Placement(transformation(extent={{-40,-34},{-20,-14}})));

        output Physiolibrary.Types.Pressure PPV = Liver.q_in.pressure "Portal vein pressure";
        output Physiolibrary.Types.Pressure HVPG = Liver.dp "Hepatic venous pressure gradient";
        output Physiolibrary.Types.Volume V_asc = levittCase1SsSiIo.Va "Ascites volume";
        output Physiolibrary.Types.Pressure P_abd = levittCase1SsSiIo.Pa "Abdominal (ascites) pressure";
        output Physiolibrary.Types.VolumeFlowRate Q_liver = Liver.q_in.q;

      equation
        connect(Liver.cond, liverConductance.y)
          annotation (Line(points={{10,6},{10,14},{-19,14}}, color={0,0,127}));
        connect(TIPSS.q_in, IntestineVenule.q_out) annotation (Line(
            points={{0,-34},{-14,-34},{-14,0},{-20,0}},
            color={0,0,0},
            thickness=1));
        connect(TIPSS.q_out, HV.q_in) annotation (Line(
            points={{20,-34},{24,-34},{24,0},{52,0}},
            color={0,0,0},
            thickness=1));
        connect(HV.q_out, q_out) annotation (Line(
            points={{72,0},{100,0}},
            color={0,0,0},
            thickness=1));
        connect(q_in, IntestinesArt.q_in) annotation (Line(
            points={{-100,0},{-80,0}},
            color={0,0,0},
            thickness=1));
        connect(tipssConductance.y, TIPSS.cond) annotation (Line(points={{-19,
                -24},{10,-24},{10,-28}}, color={0,0,127}));
        annotation (Documentation(info="<html>
<p>The TIPS resistance taken from TIPS flow and PVP from Su et al (2012, PMID 22099870).</p>
</html>"));
      end Ascites_ResistanceDynamic;

      model ResistancePressureDep
        "Hydraulic resistor, dependent on the pressure at the inflow, outflow or centered"
       extends Physiolibrary.Hydraulic.Interfaces.OnePort;
       extends Physiolibrary.Icons.HydraulicResistor;

        parameter Boolean enable=true   "if false, no resistance is used"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="External inputs/outputs"));

        parameter Boolean useConductanceInput = false
          "=true, if external conductance value is used"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="External inputs/outputs"));

        Physiolibrary.Types.HydraulicResistance R "Actual resistance value";
        constant Real ni = 4e-3 "Blood dynamic viscosity";
        parameter Modelica.Units.SI.Length L = 1;
        parameter Physiolibrary.Types.HydraulicCompliance Comp(displayUnit="l/mmHg")=7.5006157584566e-05;
        Modelica.Units.SI.Radius r(start = 1e-3) "Radius";
        Modelica.Units.SI.Diameter d=2*r "Diameter" annotation (Dialog(group="External inputs/outputs",
              enable=UsePrescribedDiameter));
        Physiolibrary.Types.Volume V=Modelica.Constants.pi*(r^2)*L "Actual volume";

        parameter Physiolibrary.Types.HydraulicResistance R_nom(displayUnit="(mmHg.min)/l")=
             79993432.449 "Nominal resistance at P_nom";
        Physiolibrary.Types.Volume V_nom=Modelica.Constants.pi*(r_nom^2)*L "Zero-pressure (nominal) volume";
        parameter Physiolibrary.Types.Pressure P_nom=0
                                                     "Nominal end-point pressure";
        // initialization
        Modelica.Units.SI.Radius r_nom(start = 1e-3)
          "Nominal radius at given pressure, inferred from the P_nom";

        parameter SideEnum side=Lymphatics.Hemodynamics.Components.SideEnum.Central
          "Side at which the filling is relative to resistance - inflow, outflow, or averaged (central)";
        Physiolibrary.Types.Pressure P_inner;

        Physiolibrary.Types.RealIO.PressureInput P_ext=p_abd   if useExternalCollapsingPressure "ExternalPressure"
          annotation (Placement(transformation(extent={{-120,70},{-80,110}})));
          Physiolibrary.Types.Pressure P_transm = (P_inner - p_abd);
          parameter Boolean useExternalCollapsingPressure = false;

        // calculation from http://www.homepages.ucl.ac.uk/~uceseug/Fluids2/Notes_Viscosity.pdf
        Modelica.Units.SI.ShearStress shearStress = 4*ni*q_in.q/r^3 "Current wall shear stress";
        Modelica.Units.SI.ShearStress shearStress0= 4*ni*(P_nom/R_nom)/r_nom^3 "nominal shear stress";

      //   Modelica.Units.SI.ShearStress shearStress_ = dp*d/L/4;
      //   Modelica.Units.SI.ShearStress shearStress0_ = P_nom*(r_nom*2)/L/4;
        parameter RemodelingModel rm = RemodelingModel.Pd;
      //  parameter Boolean useShearRemodelling = false;

      parameter Boolean UsePrescribedDiameter = false;
      protected
        type Tension =   Real (final quantity = "Tension", final unit = "N/m");
        type Radius =   Modelica.Units.SI.Radius (
                                                final nominal =   1e-6);
        type CircumferentialLength =   Modelica.Units.SI.Length (
                                                               final nominal =   1e-3);

        Physiolibrary.Types.Pressure p_abd "Abdominal pressure";


        // TENSIONS
        Tension T_p = a*f_L + b*g_L "Passive vessel wall tension";
        Real f_L = L*(L - L_0)/L_0 "Tension function on circumferential wall length";
        Real g_L = L_0*(exp(d*(L - L_0)/L_0) - 1)
          "Tension function on circumferential wall length";
        CircumferentialLength L_0 = r_nom*(2*Modelica.Constants.pi)
          "Zero-pressure circumferential wall length";


        // TENSION PARAMETERS
      //  constant Real mmHg2Pa = 133.32;
        parameter Tension a = 108.5657 "identified by Matlab's cftool and fixed";
        parameter Tension b = 0.0026 "identified by Matlab's cftool and fixed";

      equation
        if not useExternalCollapsingPressure then
          p_abd = 0;
        end if;

      //   P_nom
        R_nom = 8*ni*L/(Modelica.Constants.pi*r_nom^4);


        if UsePrescribedDiameter then
          R = 8*ni*L/(Modelica.Constants.pi*((d/2)^4));
        else
          R = 8*ni*L/(Modelica.Constants.pi*(r^4));
        end if;


        if rm == RemodelingModel.Pd then
          // simple pressure difference - correlates with the shear
          max(dp - P_nom, 0)*Comp = V - V_nom;
        elseif rm == RemodelingModel.Pt then
          // transmural pressure difference
          max(P_transm - P_nom, 0)*Comp = V - V_nom;
        elseif rm == RemodelingModel.Wt then
          // pass
          P_inner * r  =  T_p;
        else
          der(r) = -1*(shearStress - shearStress0);
        end if;






        if side == SideEnum.Left then
          P_inner = q_in.pressure;
        elseif side == SideEnum.Central then
          P_inner = (q_in.pressure + q_out.pressure)/2;
        else // side == Side.Right then
          P_inner = q_out.pressure;
        end if;

        q_in.q =(q_in.pressure - q_out.pressure)/R;


        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
                  -100,-100},{100,100}}),
                         graphics={Text(
                extent={{-220,-40},{200,-80}},
                lineColor={0,0,255},
                fillColor={58,117,175},
                fillPattern=FillPattern.Solid,
                textString="%name"),
              Line(
                points={{-100,60},{-80,60},{-20,16},{22,16},{80,60},{100,60}},
                color={0,0,0},
                thickness=0.5,
                smooth=Smooth.Bezier,
                pattern=LinePattern.Dash),
              Line(
                points={{-100,-60},{-80,-60},{-20,-16},{22,-16},{80,-60},{100,-60}},
                color={0,0,0},
                thickness=0.5,
                smooth=Smooth.Bezier,
                pattern=LinePattern.Dash),
              Line(
                points={{-100,72},{-80,72},{-20,36},{20,36},{80,72},{100,72}},
                color={0,0,0},
                thickness=0.5,
                smooth=Smooth.Bezier),
              Line(
                points={{-100,-72},{-80,-72},{-20,-36},{20,-36},{80,-72},{100,-72}},
                color={0,0,0},
                thickness=0.5,
                smooth=Smooth.Bezier)}),
          Documentation(revisions="<html>
<p><i>2009-2010</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>This hydraulic conductance (resistance) element contains two connector sides. No hydraulic medium volume is changing in this element during simulation. That means that sum of flow in both connector sides is zero. The flow through element is determined by <b>Ohm&apos;s law</b>. It is used conductance (=1/resistance) because it could be numerical zero better then infinity in resistance. </p>
</html>"),Diagram(graphics={
              Ellipse(
                extent={{38,28},{58,-28}},
                lineColor={0,140,72},
                pattern=LinePattern.Dot,
                lineThickness=0.5),
              Text(
                extent={{40,18},{86,26}},
                textColor={0,0,0},
                textString="Model"),
              Ellipse(
                extent={{42,20},{54,-18}},
                lineColor={0,0,0},
                pattern=LinePattern.Dash,
                lineThickness=0.5),
              Ellipse(
                extent={{-58,26},{-46,-24}},
                lineColor={28,108,200},
                lineThickness=0.5),
              Ellipse(
                extent={{42,14},{52,-12}},
                lineColor={28,108,200},
                lineThickness=0.5),
              Line(
                points={{-52,26},{46,14}},
                color={28,108,200},
                thickness=0.5),
              Line(
                points={{-52,-24},{48,-12}},
                color={28,108,200},
                thickness=0.5),
              Text(
                extent={{-78,10},{-64,-8}},
                textColor={28,108,200},
                textString="Pin"),
              Text(
                extent={{60,10},{82,-8}},
                textColor={28,108,200},
                textString="Pout"),
              Ellipse(
                extent={{-12,20},{0,-18}},
                lineColor={0,0,0},
                pattern=LinePattern.Dash,
                lineThickness=0.5),
              Line(
                points={{-52,20},{48,20}},
                color={0,0,0},
                pattern=LinePattern.Dash,
                thickness=0.5),
              Line(
                points={{-52,-18},{48,-18}},
                color={0,0,0},
                pattern=LinePattern.Dash,
                thickness=0.5),
              Ellipse(
                extent={{-58,20},{-46,-18}},
                lineColor={0,0,0},
                pattern=LinePattern.Dash,
                lineThickness=0.5),
              Ellipse(
                extent={{-16,28},{4,-28}},
                lineColor={0,140,72},
                pattern=LinePattern.Dot,
                lineThickness=0.5),
              Line(
                points={{-52,28},{48,28}},
                color={0,140,72},
                pattern=LinePattern.Dot,
                thickness=0.5),
              Line(
                points={{-52,-28},{48,-28}},
                color={0,140,72},
                pattern=LinePattern.Dot,
                thickness=0.5),
              Ellipse(
                extent={{-62,28},{-42,-28}},
                lineColor={0,140,72},
                pattern=LinePattern.Dot,
                lineThickness=0.5),
              Text(
                extent={{50,28},{132,36}},
                textColor={0,140,72},
                textString="Higher midpoint pressure"),
              Line(
                points={{-60,-40},{-60,-100},{0,-100}},
                color={0,0,0},
                thickness=0.5),
              Line(
                points={{-60,-80},{0,-80}},
                color={28,108,200},
                thickness=0.5,
                pattern=LinePattern.Dash),
              Line(
                points={{-20,-40},{-40,-80},{-60,-80}},
                color={28,108,200},
                thickness=0.5),
              Line(
                points={{-40,-80},{-40,-100}},
                color={28,108,200},
                thickness=0.5,
                pattern=LinePattern.Dot),
              Text(
                extent={{-34,-104},{-18,-86}},
                textColor={28,108,200},
                textString="Vnom"),
              Text(
                extent={{-58,-84},{-42,-66}},
                textColor={28,108,200},
                textString="Pnom"),
              Text(
                extent={{-24,-6},{24,6}},
                textColor={28,108,200},
                textString="Transmural pressure",
                origin={-64,-74},
                rotation=90),
              Text(
                extent={{-24,-6},{24,6}},
                textColor={28,108,200},
                origin={-32,-106},
                rotation=360,
                textString="Shunts volume"),
              Line(
                points={{20,-40},{20,-100},{80,-100}},
                color={0,0,0},
                thickness=0.5),
              Text(
                extent={{-24,-6},{24,6}},
                textColor={28,108,200},
                origin={52,-108},
                rotation=360,
                textString="Shunts volume"),
              Text(
                extent={{-24,-6},{24,6}},
                textColor={28,108,200},
                origin={14,-72},
                rotation=90,
                textString="Resistance"),
              Line(
                points={{40,-72},{64,-94},{86,-98},{102,-98}},
                color={28,108,200},
                thickness=0.5,
                smooth=Smooth.Bezier),
              Text(
                extent={{42,-104},{58,-86}},
                textColor={28,108,200},
                textString="Vnom"),
              Line(
                points={{40,-72},{40,-100}},
                color={28,108,200},
                thickness=0.5,
                pattern=LinePattern.Dot)}));
      end ResistancePressureDep;

      type SideEnum = enumeration(
          Left                     "Left (inflow) side",
          Central                                                "Central",
          Right                                                                   "Right (outflow)")
        "Side of a resistance";
      model SplanchnicCirculation
        Physiolibrary.Hydraulic.Components.Resistor splanchnic(Resistance=
              7999343244.900001*(93/500))
          annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
        Physiolibrary.Hydraulic.Components.Resistor HV(Resistance=79993432.449)
          annotation (Placement(transformation(extent={{40,-10},{60,10}})));
        Physiolibrary.Hydraulic.Components.ElasticVessel Cliver(
          volume_start=0.001,
          ZeroPressureVolume=0.0003,
          Compliance=7.500615758456563e-09*(600/8))
          annotation (Placement(transformation(extent={{-10,10},{10,30}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_b q_out annotation (
            Placement(transformation(rotation=0, extent={{90,-50},{110,-30}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a q_in annotation (
            Placement(transformation(rotation=0, extent={{-110,-50},{-90,-30}})));
      equation
        connect(splanchnic.q_out,Cliver. q_in) annotation (Line(
            points={{-40,0},{2.22045e-16,0},{2.22045e-16,20}},
            color={0,0,0},
            thickness=1));
        connect(HV.q_in,Cliver. q_in) annotation (Line(
            points={{40,0},{2.22045e-16,0},{2.22045e-16,20},{0,20}},
            color={0,0,0},
            thickness=1));
        connect(q_out, HV.q_out) annotation (Line(points={{100,-40},{80,-40},{80,
                0},{60,0}}, color={0,0,0}));
        connect(q_in, splanchnic.q_in) annotation (Line(points={{-100,-40},{-80,
                -40},{-80,0},{-60,0}}, color={0,0,0}));
      end SplanchnicCirculation;

      model HagenPoiseulleConductance ""
        Physiolibrary.Types.RealIO.HydraulicConductanceOutput hydraulicconductance = 1/R
          annotation (Placement(transformation(extent={{92,-10},{112,10}})));

       Physiolibrary.Types.HydraulicResistance R = 8*mu*L/(Modelica.Constants.pi*(d/2)^4);
       Modelica.Units.SI.Diameter d;
       parameter Modelica.Units.SI.Diameter d_nominal(displayUnit="mm")=0.001;

       parameter Real mu = 3e-3;

       parameter Physiolibrary.Types.Length L=0.05;
       parameter Boolean useNominalDiameter = true;
       parameter Physiolibrary.Types.Pressure dp_nominal=1333.22387415;
       parameter Physiolibrary.Types.VolumeFlowRate q_nominal=
            3.3333333333333e-06;
      equation
        if useNominalDiameter then
          d = d_nominal;
        else
          dp_nominal = R*q_nominal;
        end if;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end HagenPoiseulleConductance;

      model ResistanceControlled "Integral-controlled resistance to maintain nominal flow rate"
       //extends Physiolibrary.Hydraulic.Interfaces.OnePort;

        extends Physiolibrary.Icons.Resistor;
        Physiolibrary.Types.VolumeFlowRate volumeFlowRate = q_in.q;
        Physiolibrary.Types.Pressure dp = q_in. pressure - q_out.pressure;
        parameter Physiolibrary.Types.HydraulicResistance Resistance = 1;
        parameter Physiolibrary.Types.VolumeFlowRate Q_nominal = 1;
        parameter Physiolibrary.Types.Time tau = 1;
        Physiolibrary.Hydraulic.Components.Resistor resistor(useConductanceInput=true,
            Resistance=7999343244.9)
          annotation (Placement(transformation(extent={{-22,-10},{-2,10}})));
        Physiolibrary.Hydraulic.Sensors.FlowMeasure flowMeasure
          annotation (Placement(transformation(extent={{40,-10},{60,10}})));
        Modelica.Blocks.Sources.RealExpression realExpression(y=c)
          annotation (Placement(transformation(extent={{-44,14},{-24,34}})));

        Physiolibrary.Types.HydraulicConductance c( start = 1/Resistance);

        Real action = (Q_nominal-flowMeasure.volumeFlow);
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a
                             q_in "Volume inflow" annotation (Placement(
              transformation(extent={{-114,-14},{-86,14}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_b
                             q_out "Volume outflow"
                               annotation (Placement(
              transformation(extent={{86,-14},{114,14}})));
      equation

        der(c)*tau = action;

        connect(resistor.q_out, flowMeasure.q_in) annotation (Line(
            points={{-2,0},{40,0}},
            color={0,0,0},
            thickness=1));
        connect(realExpression.y, resistor.cond)
          annotation (Line(points={{-23,24},{-12,24},{-12,6}}, color={0,0,127}));
        connect(flowMeasure.q_out, q_out) annotation (Line(
            points={{60,0},{100,0}},
            color={0,0,0},
            thickness=1));
        connect(resistor.q_in, q_in) annotation (Line(
            points={{-22,0},{-100,0}},
            color={0,0,0},
            thickness=1));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end ResistanceControlled;

      model LiverPortalVeinTree
        extends Physiolibrary.Hydraulic.Interfaces.partialOnePort;
        Physiolibrary.Hydraulic.Components.Resistor portalTrunk(Resistance(
              displayUnit="(mmHg.min)/l") = 799934.32449)
          annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
        Physiolibrary.Hydraulic.Components.Resistor leftPortalVein(Resistance(
              displayUnit="(mmHg.min)/l") = 799934.32449)
          annotation (Placement(transformation(extent={{-38,-70},{-18,-50}})));
        Physiolibrary.Hydraulic.Components.Resistor leftLateral(
            useConductanceInput=true,                           Resistance=0.17
              *R_tot)
          annotation (Placement(transformation(extent={{2,-50},{22,-30}})));
        Physiolibrary.Hydraulic.Components.Resistor leftMedial(
            useConductanceInput=true,                          Resistance=0.14*
              R_tot)
          annotation (Placement(transformation(extent={{2,-84},{22,-64}})));
        Physiolibrary.Hydraulic.Components.Resistor rightPosterior(
            useConductanceInput=true,                              Resistance=
              0.3*R_tot)
          annotation (Placement(transformation(extent={{2,20},{22,40}})));
        Physiolibrary.Hydraulic.Components.Resistor rightAnterior(
            useConductanceInput=true,                             Resistance=
              0.37*R_tot)
          annotation (Placement(transformation(extent={{2,-10},{22,10}})));
        parameter Physiolibrary.Types.HydraulicResistance R_tot=0
          "Total liver resistance";
        Modelica.Blocks.Math.Gain gain(k=1/0.3)
          annotation (Placement(transformation(extent={{48,26},{28,46}})));
        Modelica.Blocks.Math.Gain gain1(k=1/0.37)
          annotation (Placement(transformation(extent={{48,-4},{28,16}})));
        Modelica.Blocks.Math.Gain gain2(k=1/0.17)
          annotation (Placement(transformation(extent={{48,-42},{28,-22}})));
        Modelica.Blocks.Math.Gain gain3(k=1/0.14)
          annotation (Placement(transformation(extent={{48,-78},{28,-58}})));
        Physiolibrary.Types.RealIO.HydraulicConductanceInput
                                               cond(start=Conductance)=c if useConductanceInput
                                                         annotation (Placement(
              transformation(extent={{-20,-20},{20,20}},
              rotation=270,
              origin={2,60})));
      equation
        connect(q_in, portalTrunk.q_in) annotation (Line(
            points={{-100,0},{-80,0}},
            color={0,0,0},
            thickness=1));
        connect(leftPortalVein.q_out, leftLateral.q_in) annotation (Line(
            points={{-18,-60},{-4,-60},{-4,-40},{2,-40}},
            color={0,0,0},
            thickness=1));
        connect(leftPortalVein.q_out, leftMedial.q_in) annotation (Line(
            points={{-18,-60},{-4,-60},{-4,-74},{2,-74}},
            color={0,0,0},
            thickness=1));
        connect(portalTrunk.q_out, rightAnterior.q_in) annotation (Line(
            points={{-60,0},{2,0}},
            color={0,0,0},
            thickness=1));
        connect(portalTrunk.q_out, rightPosterior.q_in) annotation (Line(
            points={{-60,0},{-52,0},{-52,30},{2,30}},
            color={0,0,0},
            thickness=1));
        connect(leftPortalVein.q_in, rightAnterior.q_in) annotation (Line(
            points={{-38,-60},{-52,-60},{-52,0},{2,0}},
            color={0,0,0},
            thickness=1));
        connect(leftLateral.q_out, q_out) annotation (Line(
            points={{22,-40},{56,-40},{56,0},{62,0},{62,1.77636e-15},{100,
                1.77636e-15}},
            color={0,0,0},
            thickness=1));
        connect(leftMedial.q_out, q_out) annotation (Line(
            points={{22,-74},{56,-74},{56,0},{62,0},{62,1.77636e-15},{100,
                1.77636e-15}},
            color={0,0,0},
            thickness=1));
        connect(rightPosterior.q_out, q_out) annotation (Line(
            points={{22,30},{56,30},{56,1.77636e-15},{100,1.77636e-15}},
            color={0,0,0},
            thickness=1));
        connect(rightAnterior.q_out, q_out) annotation (Line(
            points={{22,0},{100,0}},
            color={0,0,0},
            thickness=1));
        connect(gain.y, rightPosterior.cond)
          annotation (Line(points={{27,36},{12,36}}, color={0,0,127}));
        connect(gain1.y, rightAnterior.cond)
          annotation (Line(points={{27,6},{12,6}}, color={0,0,127}));
        connect(gain2.y, leftLateral.cond) annotation (Line(points={{27,-32},{
                27,-34},{12,-34}}, color={0,0,127}));
        connect(gain3.y, leftMedial.cond)
          annotation (Line(points={{27,-68},{12,-68}}, color={0,0,127}));
        connect(gain.u, cond) annotation (Line(points={{50,36},{54,36},{54,60},
                {2,60}}, color={0,0,127}));
        connect(gain1.u, cond) annotation (Line(points={{50,6},{66,6},{66,60},{
                2,60}}, color={0,0,127}));
        connect(gain2.u, cond) annotation (Line(points={{50,-32},{58,-32},{58,
                -34},{66,-34},{66,60},{2,60}}, color={0,0,127}));
        connect(gain3.u, cond) annotation (Line(points={{50,-68},{66,-68},{66,
                60},{2,60}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end LiverPortalVeinTree;

      model ResistancePressureDepTensions
        "Hydraulic resistor, dependent on the pressure at the inflow, outflow or centered"
       extends Physiolibrary.Hydraulic.Interfaces.OnePort;
       extends Physiolibrary.Icons.HydraulicResistor;

        parameter Boolean enable=true   "if false, no resistance is used"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="External inputs/outputs"));

        parameter Boolean useConductanceInput = false
          "=true, if external conductance value is used"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="External inputs/outputs"));

        Physiolibrary.Types.HydraulicResistance R "Actual resistance value";
        constant Real ni = 4e-3 "Blood dynamic viscosity";
        constant Modelica.Units.SI.Length L = 1;
        parameter Physiolibrary.Types.HydraulicCompliance Comp(displayUnit="l/mmHg")=7.5006157584566e-05;
        Modelica.Units.SI.Radius r "Radius";
        Physiolibrary.Types.Volume V = Modelica.Constants.pi * (r - R0)^2*L;

        Modelica.Units.SI.Radius R0 "Zero resistance, inferred from the P_nom at R0_nom";
        Modelica.Units.SI.Radius R_nom "Zero resistance, inferred from the P_nom at R0_nom";
        parameter Physiolibrary.Types.HydraulicResistance r0_nom(displayUnit="(mmHg.min)/l")=
           79993432.449                                          "Nominal resistance at P_nom";
        parameter Physiolibrary.Types.Pressure P_nom=1333.22387415
                                                     "Nominal end-point pressure";

        parameter SideEnum side=Lymphatics.Hemodynamics.Components.SideEnum.Central
          "Side at which the filling is relative to resistance - inflow, outflow, or averaged (central)";
        Physiolibrary.Types.Pressure P_inner;
        Physiolibrary.Hydraulic.Components.Resistor resistor
          annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
        Physiolibrary.Hydraulic.Components.Resistor resistor1
          annotation (Placement(transformation(extent={{40,-10},{60,10}})));
        ComplianceTube complianceTube annotation (Placement(transformation(
                rotation=0, extent={{-20,60},{0,80}})));
      equation
        // to calculate the R_0 parameter
        Modelica.Constants.pi * (R_nom - R0)^2*L = Comp*P_nom;
        r0_nom = 8*ni*L/(Modelica.Constants.pi*R_nom^4);

        R = 8*ni*L/(Modelica.Constants.pi*r^4);

        V =Comp*P_inner;

        if side == SideEnum.Left then
          P_inner = q_in.pressure;
        elseif side == SideEnum.Central then
          P_inner = (q_in.pressure + q_out.pressure)/2;
        else // side == Side.Right then
          P_inner = q_out.pressure;
        end if;

        q_in.q =(q_in.pressure - q_out.pressure)/R;
        connect(q_in, resistor.q_in) annotation (Line(
            points={{-100,0},{-60,0}},
            color={0,0,0},
            thickness=1));
        connect(q_out, resistor1.q_out) annotation (Line(
            points={{100,0},{60,0}},
            color={0,0,0},
            thickness=1));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
                  -100,-100},{100,100}}),
                         graphics={Text(
                extent={{-220,-40},{200,-80}},
                lineColor={0,0,255},
                fillColor={58,117,175},
                fillPattern=FillPattern.Solid,
                textString="%name"),
              Line(
                points={{-100,60},{-80,60},{-20,16},{22,16},{80,60},{100,60}},
                color={0,0,0},
                thickness=0.5,
                smooth=Smooth.Bezier,
                pattern=LinePattern.Dash),
              Line(
                points={{-100,-60},{-80,-60},{-20,-16},{22,-16},{80,-60},{100,-60}},
                color={0,0,0},
                thickness=0.5,
                smooth=Smooth.Bezier,
                pattern=LinePattern.Dash),
              Line(
                points={{-100,72},{-80,72},{-20,36},{20,36},{80,72},{100,72}},
                color={0,0,0},
                thickness=0.5,
                smooth=Smooth.Bezier),
              Line(
                points={{-100,-72},{-80,-72},{-20,-36},{20,-36},{80,-72},{100,-72}},
                color={0,0,0},
                thickness=0.5,
                smooth=Smooth.Bezier)}),
          Documentation(revisions="<html>
<p><i>2009-2010</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>This hydraulic conductance (resistance) element contains two connector sides. No hydraulic medium volume is changing in this element during simulation. That means that sum of flow in both connector sides is zero. The flow through element is determined by <b>Ohm&apos;s law</b>. It is used conductance (=1/resistance) because it could be numerical zero better then infinity in resistance. </p>
</html>"));
      end ResistancePressureDepTensions;

      model ComplianceTube
        extends Physiolibrary.Icons.ElasticBalloon;
        ADAN_main.Components.Subsystems.Systemic.Vessel_modules.Auxiliary.Compliances.compliance_tensionDataFit
          compliance_tensionDataFit(
          V = volume,
          phi = 0,
          l=l,
          p0=1333.22387415,
          r_n=r_n)   annotation (Placement(transformation(extent={{-20,60},{0,80}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a port_a
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Physiolibrary.Types.RealIO.HydraulicConductanceOutput hydraulicconductance "Hagen-poiseulle model of resistance"
          annotation (Placement(transformation(extent={{80,-10},{100,10}})));

      Physiolibrary.Types.Volume volume(start = Modelica.Constants.pi*r_n^2*l);
      parameter Physiolibrary.Types.Length l = 0.02;

       parameter Real mu = 3e-3;

        parameter Modelica.Units.SI.Radius r_n=0.005 "nominal vessel radius";
      equation
        hydraulicconductance = (Modelica.Constants.pi*(compliance_tensionDataFit.r)^4)/(8*mu*l);

        port_a.q = der(volume);

        port_a.pressure = compliance_tensionDataFit.p;
        annotation (Icon(graphics={Line(
                points={{-98,-116},{-96,-90},{-82,-60},{14,-62},{58,24},{70,86}},
                color={217,67,180},
                thickness=1,
                smooth=Smooth.Bezier)}));
      end ComplianceTube;

      package Tests

        model Ascites_PG
          extends Components.partialAscites(
                                       Liver(useExternalPressureHead=true));
          Modelica.Blocks.Sources.RealExpression realExpression(y=-time*133.322)
            annotation (Placement(transformation(extent={{-52,-52},{-32,-32}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
            annotation (Placement(transformation(extent={{100,-10},{80,10}})));
        equation
          connect(Liver.p_headIn,realExpression. y) annotation (Line(points={{10,10},
                  {-24,10},{-24,-42},{-31,-42}},
                                           color={0,0,127}));
          connect(CVP.y, HV.q_out) annotation (Line(
              points={{80,0},{72,0}},
              color={0,0,0},
              thickness=1));
        end Ascites_PG;

        model Ascites_Resistance
          extends Components.partialAscites(
            redeclare Physiolibrary.Hydraulic.Components.Resistor Liver(
                useConductanceInput=true, Resistance=6*mmHg/Qnom),
            redeclare Physiolibrary.Hydraulic.Components.Resistor IntestinesArt(
                Resistance=84*mmHg/Qnom),
            redeclare Physiolibrary.Hydraulic.Components.Resistor IntestineVenule(
                Resistance=3*mmHg/Qnom));
          Modelica.Blocks.Sources.RealExpression realExpression1(y=1/(time*mmHg/Qnom))
            annotation (Placement(transformation(extent={{-32,18},{-12,38}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume unlimitedVolume2(P=666.611937075)
            annotation (Placement(transformation(extent={{-176,-20},{-196,0}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume unlimitedVolume1(P=13332.2387415)
            annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

          parameter Physiolibrary.Types.VolumeFlowRate Qnom=1.666666666666667e-08*(5000*
              0.2) "Nominal flow through the splanchnic circulation";
          parameter Physiolibrary.Types.Pressure mmHg=133.322;
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
            annotation (Placement(transformation(extent={{110,-10},{90,10}})));
          Physiolibrary.Hydraulic.Components.Resistor TIPSS(
            enable=true,
            useConductanceInput=false,
            Resistance(displayUnit="(mmHg.min)/l") = 7999343.2449*(15/1.5*1e15))
            annotation (Placement(transformation(extent={{-2,-50},{18,-30}})));
          Components.Ascites_Resistance ascites_Resistance
            annotation (Placement(transformation(extent={{-2,-86},{18,-66}})));
          Components.SplanchnicCirculation splanchnicCirculation
            annotation (Placement(transformation(extent={{-2,-130},{18,-110}})));
        equation
          connect(Liver.cond, realExpression1.y) annotation (Line(points={{10,6},{
                  10,14},{-4,14},{-4,28},{-11,28}},
                                             color={0,0,127}));
          connect(IntestinesArt.q_in, unlimitedVolume1.y) annotation (Line(
              points={{-80,0},{-90,0}},
              color={0,0,0},
              thickness=1));
          connect(CVP.y, HV.q_out) annotation (Line(
              points={{90,0},{72,0}},
              color={0,0,0},
              thickness=1));
          connect(TIPSS.q_in, IntestineVenule.q_out) annotation (Line(
              points={{-2,-40},{-14,-40},{-14,2.22045e-16},{-20,2.22045e-16}},
              color={0,0,0},
              thickness=1));
          connect(TIPSS.q_out, HV.q_in) annotation (Line(
              points={{18,-40},{24,-40},{24,2.22045e-16},{52,2.22045e-16}},
              color={0,0,0},
              thickness=1));
          connect(ascites_Resistance.q_in, unlimitedVolume1.y) annotation (Line(
              points={{-2,-76},{-86,-76},{-86,0},{-90,0}},
              color={0,0,0},
              thickness=1));
          connect(ascites_Resistance.q_out, HV.q_out) annotation (Line(
              points={{18,-76},{86,-76},{86,0},{72,0}},
              color={0,0,0},
              thickness=1));
          connect(splanchnicCirculation.q_out, HV.q_out) annotation (Line(
              points={{18,-124},{34,-124},{34,-76},{86,-76},{86,0},{72,0}},
              color={0,0,0},
              thickness=1));
          connect(splanchnicCirculation.q_in, unlimitedVolume1.y) annotation (
              Line(
              points={{-2,-124},{-22,-124},{-22,-74},{-86,-74},{-86,0},{-90,0}},
              color={0,0,0},
              thickness=1));
          annotation (Documentation(info="<html>
<p>The TIPS resistance taken from TIPS flow and PVP from Su et al (2012, PMID 22099870).</p>
</html>"));
        end Ascites_Resistance;

        model TestResistance
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume unlimitedVolume1(
              usePressureInput=true, P=13332.2387415)
            annotation (Placement(transformation(extent={{-104,-10},{-84,10}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(usePressureInput=
                true, P=666.611937075)
            annotation (Placement(transformation(extent={{116,-10},{96,10}})));
          Components.ResistancePressureDep resistancePressureDep_trans(
            L=0.01,
            R_nom=799934324.49,
            Comp(displayUnit="ml/mmHg") = 7.5006157584566e-07,
            P_nom=13332.2387415,
            side=Lymphatics.Hemodynamics.Components.SideEnum.Right,
            useExternalCollapsingPressure=true,
            rm=Lymphatics.Hemodynamics.Components.RemodelingModel.Pt,
            r(start=0.0005973604195040681),
            r_nom(start=0.0005973604195040681))
            annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
          Modelica.Blocks.Sources.RealExpression realExpression(y=max(time, 5)*
                133.322)
            annotation (Placement(transformation(extent={{92,30},{112,50}})));
          Modelica.Blocks.Sources.RealExpression realExpression1(y=(10*sin(6.28*
                time) + 100)*133.322)
            annotation (Placement(transformation(extent={{-144,34},{-124,54}})));
          ResistancePressureDep resistancePressureDep_dp(
            L=0.01,
            Comp(displayUnit="ml/mmHg") = 7.5006157584566e-12,
            R_nom=799934324.49,
            P_nom=13332.2387415,
            r(start=0.0005973604195040681),
            r_nom(start=0.0005973604195040681))
            annotation (Placement(transformation(extent={{-10,22},{10,42}})));
          ResistancePressureDep resistancePressureDep_shear(
            L=0.01,
            R_nom=799934324.49,
            P_nom=13332.2387415,
            rm=Lymphatics.Hemodynamics.Components.RemodelingModel.Shear,
            r_nom(start=0.0005973604195040681))
            annotation (Placement(transformation(extent={{-12,-40},{8,-20}})));
          Modelica.Blocks.Sources.RealExpression realExpression2(y=10*133.322)
            annotation (Placement(transformation(extent={{-42,4},{-22,24}})));
        equation
          connect(unlimitedVolume1.y, resistancePressureDep_trans.q_in)
            annotation (Line(
              points={{-84,0},{-10,0}},
              color={0,0,0},
              thickness=1));
          connect(CVP.y, resistancePressureDep_trans.q_out) annotation (Line(
              points={{96,0},{10,0}},
              color={0,0,0},
              thickness=1));
          connect(realExpression.y, CVP.pressure) annotation (Line(points={{113,
                  40},{150,40},{150,0},{116,0}}, color={0,0,127}));
          connect(unlimitedVolume1.pressure, realExpression1.y) annotation (Line(
                points={{-104,0},{-123,0},{-123,44}}, color={0,0,127}));
          connect(unlimitedVolume1.y, resistancePressureDep_dp.q_in)
            annotation (Line(
              points={{-84,0},{-48,0},{-48,32},{-10,32}},
              color={0,0,0},
              thickness=1));
          connect(CVP.y, resistancePressureDep_dp.q_out) annotation (Line(
              points={{96,0},{54,0},{54,32},{10,32}},
              color={0,0,0},
              thickness=1));
          connect(unlimitedVolume1.y, resistancePressureDep_shear.q_in)
            annotation (Line(
              points={{-84,0},{-48,0},{-48,-30},{-12,-30}},
              color={0,0,0},
              thickness=1));
          connect(CVP.y, resistancePressureDep_shear.q_out) annotation (Line(
              points={{96,0},{52,0},{52,-30},{8,-30}},
              color={0,0,0},
              thickness=1));
          connect(realExpression2.y, resistancePressureDep_trans.P_ext)
            annotation (Line(points={{-21,14},{-10,14},{-10,9}}, color={0,0,127}));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)),
            experiment(
              StartTime=5,
              StopTime=10,
              Tolerance=1e-06,
              __Dymola_Algorithm="Cvode"));
        end TestResistance;

        model AscitesResistanceBased
          "Levitts ascites with nominal hemodynamics break down"
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume unlimitedVolume(P=666.611937075)
            annotation (Placement(transformation(extent={{100,-10},{80,10}})));
          AscitesLevitt.LevittCase1SsSiIo levittCase1SsSiIo
            annotation (Placement(transformation(extent={{30,54},{50,74}})));
          Physiolibrary.Hydraulic.Components.PumpPressureHead
                            HV(useExternalCollapsingPressure=true, p_head0=-266.64477483)
            "Hepatic vein (free)"
            annotation (Placement(transformation(extent={{20,-10},{40,10}})));
          Physiolibrary.Hydraulic.Sensors.PressureMeasure PRa "Right atrial pressure"
            annotation (Placement(transformation(extent={{80,40},{60,60}})));
          Physiolibrary.Hydraulic.Components.PumpPressureHead Liver(
              useExternalPressureHead=true, p_head0=-799.93432449)
            annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
          Physiolibrary.Hydraulic.Components.PumpPressureHead IntestineVenule(p_head0=-399.967162245)
            "Intestinal venule pressure drop"
            annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
          Physiolibrary.Hydraulic.Sensors.PressureMeasure Ppv
            annotation (Placement(transformation(extent={{-40,58},{-20,78}})));
          Physiolibrary.Hydraulic.Components.PumpPressureHead IntestinesArt(p_head0=-11199.08054286)
            "Intestinal Arteriole resistance"
            annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
          Physiolibrary.Hydraulic.Sensors.PressureMeasure Pc
            annotation (Placement(transformation(extent={{-68,66},{-48,86}})));
          Physiolibrary.Hydraulic.Sensors.PressureMeasure Phv
            annotation (Placement(transformation(extent={{-2,50},{18,70}})));
          Modelica.Blocks.Sources.RealExpression realExpression(y=-time*133.322)
            annotation (Placement(transformation(extent={{-94,12},{-74,32}})));
          parameter Physiolibrary.Types.VolumeFlowRate Qnom=1.666666666666667e-08*(
              5000*0.2) "Nominal flow through the splanchnic circulation";
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume unlimitedVolume1(P=
                13332.2387415)
            annotation (Placement(transformation(extent={{-130,-30},{-110,-10}})));
          Physiolibrary.Hydraulic.Components.Resistor R_intestArt(Resistance=84*
                mmHg/Qnom)
            annotation (Placement(transformation(extent={{-100,-30},{-80,-10}})));
          parameter Physiolibrary.Types.Pressure mmHg=133.322;
          Physiolibrary.Hydraulic.Components.Resistor R_IntestVen(Resistance=3*mmHg
                /Qnom)
            annotation (Placement(transformation(extent={{-60,-30},{-40,-10}})));
          Physiolibrary.Hydraulic.Components.Resistor R_liver(useConductanceInput=
                true, Resistance=6*mmHg/Qnom)
            annotation (Placement(transformation(extent={{-20,-10},{0,-30}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume unlimitedVolume2(P=
                666.611937075)
            annotation (Placement(transformation(extent={{100,-30},{80,-10}})));
          Physiolibrary.Hydraulic.Components.PumpPressureHead PG_HV(
              useExternalCollapsingPressure=true, p_head0=-266.64477483)
            "Hepatic vein (free)"
            annotation (Placement(transformation(extent={{20,-30},{40,-10}})));
          Modelica.Blocks.Sources.RealExpression realExpression1(y=1/(time*mmHg/
                Qnom))
            annotation (Placement(transformation(extent={{-72,-70},{-52,-50}})));
          Physiolibrary.Hydraulic.Components.Resistor R_shunt(Resistance(
                displayUnit="(mmHg.min)/l") = 7999343.2449*(15/1.5))
            annotation (Placement(transformation(extent={{-18,-84},{2,-64}})));
        equation
          connect(unlimitedVolume.y, PRa.q_in) annotation (Line(
              points={{80,0},{74,0},{74,44}},
              color={0,0,0},
              thickness=1));
          connect(PRa.pressure, levittCase1SsSiIo.Pra) annotation (Line(points={{64,46},
                  {56,46},{56,64},{50,64}}, color={0,0,127}));
          connect(HV.q_out, unlimitedVolume.y) annotation (Line(
              points={{40,0},{80,0}},
              color={0,0,0},
              thickness=1));
          connect(Liver.q_out, HV.q_in) annotation (Line(
              points={{0,0},{20,0}},
              color={0,0,0},
              thickness=1));
          connect(IntestineVenule.q_out, Liver.q_in) annotation (Line(
              points={{-40,0},{-20,0}},
              color={0,0,0},
              thickness=1));
          connect(Ppv.q_in, IntestineVenule.q_out) annotation (Line(
              points={{-34,62},{-34,0},{-40,0}},
              color={0,0,0},
              thickness=1));
          connect(levittCase1SsSiIo.Pp, Ppv.pressure)
            annotation (Line(points={{30,64},{-24,64}}, color={0,0,127}));
          connect(IntestinesArt.q_out, IntestineVenule.q_in) annotation (Line(
              points={{-80,0},{-60,0}},
              color={0,0,0},
              thickness=1));
          connect(IntestinesArt.q_out, Pc.q_in) annotation (Line(
              points={{-80,0},{-62,0},{-62,70}},
              color={0,0,0},
              thickness=1));
          connect(Liver.q_out,Phv. q_in) annotation (Line(
              points={{0,0},{4,0},{4,54}},
              color={0,0,0},
              thickness=1));
          connect(Phv.pressure, levittCase1SsSiIo.Phv)
            annotation (Line(points={{14,56},{30,56}}, color={0,0,127}));
          connect(Pc.pressure, levittCase1SsSiIo.Pc)
            annotation (Line(points={{-52,72},{30,72}}, color={0,0,127}));
          connect(Liver.p_headIn, realExpression.y) annotation (Line(points={{-10,10},
                  {-10,22},{-73,22}},      color={0,0,127}));
          connect(HV.EP, levittCase1SsSiIo.Pa) annotation (Line(points={{20,8},{20,
                  40},{40,40},{40,56}}, color={0,0,127}));
          connect(unlimitedVolume1.y, R_intestArt.q_in) annotation (Line(
              points={{-110,-20},{-100,-20}},
              color={0,0,0},
              thickness=1));
          connect(R_IntestVen.q_in, R_intestArt.q_out) annotation (Line(
              points={{-60,-20},{-80,-20}},
              color={0,0,0},
              thickness=1));
          connect(R_liver.q_in, R_IntestVen.q_out) annotation (Line(
              points={{-20,-20},{-40,-20}},
              color={0,0,0},
              thickness=1));
          connect(PG_HV.q_in, R_liver.q_out) annotation (Line(
              points={{20,-20},{0,-20}},
              color={0,0,0},
              thickness=1));
          connect(PG_HV.q_out, unlimitedVolume2.y) annotation (Line(
              points={{40,-20},{80,-20}},
              color={0,0,0},
              thickness=1));
          connect(PG_HV.EP, HV.EP) annotation (Line(points={{20,-12},{14,-12},{14,8},
                  {20,8}}, color={0,0,127}));
          connect(realExpression1.y, R_liver.cond) annotation (Line(points={{-51,
                  -60},{-10,-60},{-10,-26}}, color={0,0,127}));
          connect(R_shunt.q_in, R_IntestVen.q_out) annotation (Line(
              points={{-18,-74},{-30,-74},{-30,-20},{-40,-20}},
              color={0,0,0},
              thickness=1));
          connect(R_shunt.q_out, R_liver.q_out) annotation (Line(
              points={{2,-74},{6,-74},{6,-20},{0,-20}},
              color={0,0,0},
              thickness=1));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)),
            experiment(
              StartTime=6,
              StopTime=25,
              __Dymola_Algorithm="Dassl"));
        end AscitesResistanceBased;

        model ResistanceControlTest
          Components.ResistanceControlled resistanceControlled(
            Resistance(displayUnit="(mmHg.min)/l") = 799934324.49,
            Q_nominal(displayUnit="l/min") = 8.3333333333333e-05,
            tau=60)
            annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume unlimitedVolume1(P=
                13332.2387415)
            annotation (Placement(transformation(extent={{-110,-12},{-90,8}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
            annotation (Placement(transformation(extent={{110,-12},{90,8}})));
          Components.Ascites_Resistance ascites_Resistance(
            liverConductance(y=1),
            Qnom(displayUnit="l/min") = 1.6666666666667e-05,
            TIPSS(enable=true, Resistance=Modelica.Constants.inf),
            Liver(
              enable=true,
              useConductanceInput=false,
              Resistance(displayUnit="(mmHg.min)/l") = 119990148.6735),
            redeclare Components.ResistanceControlled IntestinesArt(
              Resistance(displayUnit="(mmHg.min)/l") = 399967162.245,
              Q_nominal(displayUnit="l/min") = 1.6666666666667e-05,
              tau(displayUnit="s") = 1))
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
        equation
          connect(unlimitedVolume1.y, resistanceControlled.q_in) annotation (Line(
              points={{-90,-2},{-10,-2}},
              color={0,0,0},
              thickness=1));
          connect(resistanceControlled.q_out, CVP.y) annotation (Line(
              points={{10,-2},{90,-2}},
              color={0,0,0},
              thickness=1));
          connect(ascites_Resistance.q_in, resistanceControlled.q_in) annotation (
             Line(
              points={{-10,-40},{-44,-40},{-44,-42},{-82,-42},{-82,-2},{-10,-2}},
              color={0,0,0},
              thickness=1));

          connect(ascites_Resistance.q_out, CVP.y) annotation (Line(
              points={{10,-40},{74,-40},{74,-2},{90,-2}},
              color={0,0,0},
              thickness=1));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));
        end ResistanceControlTest;

        model TestComplianceTube
          Components.ComplianceTube complianceTube(l=0.2, r_n(displayUnit="mm")=
                 0.003)
            annotation (Placement(transformation(extent={{-8,-10},{12,10}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedPump unlimitedPump(SolutionFlow(
                displayUnit="ml/min")=1E-06)
            annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
          Physiolibrary.Hydraulic.Components.Resistor resistor(useConductanceInput=true)
            annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
          inner ADAN_main.Components.Settings
                         settings(
            phi0=0,
            veins_UsePhiEffect=false,
            veins_delayed_activation=false)
            annotation (Placement(transformation(extent={{-88,46},{-68,66}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedOutflowPump
            unlimitedOutflowPump(SolutionFlow=0)
            annotation (Placement(transformation(extent={{64,-40},{84,-20}})));
          Components.ResistancePressureDep resistancePressureDep(
            R_nom(displayUnit="(mmHg.min)/ml") = 7999340.0,
            L=0.2,
            side=Lymphatics.Hemodynamics.Components.SideEnum.Left,
            Comp(displayUnit="ml/mmHg") = 7.50062E-09,
            r(start=0.001))
            annotation (Placement(transformation(extent={{18,-40},{38,-20}})));
        equation
          connect(resistor.q_out, complianceTube.port_a) annotation (Line(
              points={{-20,0},{-12,0},{2,0}},
              color={0,0,0},
              thickness=1,
              smooth=Smooth.Bezier));
          connect(resistor.q_in, unlimitedPump.q_out) annotation (Line(
              points={{-40,0},{-60,0}},
              color={0,0,0},
              thickness=1,
              smooth=Smooth.Bezier));
          connect(complianceTube.hydraulicconductance, resistor.cond) annotation (Line(
              points={{11,0},{20,0},{20,16},{-30,16},{-30,6}},
              color={0,0,127},
              smooth=Smooth.Bezier));
          connect(resistancePressureDep.q_in, complianceTube.port_a) annotation (
              Line(
              points={{18,-30},{2,-30},{2,0}},
              color={0,0,0},
              thickness=1,
              smooth=Smooth.Bezier));
          connect(unlimitedOutflowPump.q_in, resistancePressureDep.q_out)
            annotation (Line(
              points={{64,-30},{38,-30}},
              color={0,0,0},
              thickness=1,
              smooth=Smooth.Bezier));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));
        end TestComplianceTube;
      end Tests;

      package Schematics
        model ModelSchematicsSmall
          AscitesLevitt.LevittCase1SsSiIo levittCase1SsSiIo
            annotation (Placement(transformation(extent={{24,8},{44,28}})));
          Physiolibrary.Hydraulic.Components.PumpPressureHead hepatic_vein(
              useExternalCollapsingPressure=true, p_head0=-266.64477483)
            "Hepatic vein (free)"
            annotation (Placement(transformation(extent={{34,-12},{54,8}})));
          Physiolibrary.Hydraulic.Sensors.PressureMeasure PRa "Right atrial pressure"
            annotation (Placement(transformation(extent={{68,14},{48,34}})));
          Physiolibrary.Hydraulic.Components.Resistor         Liver(
              useConductanceInput=true, Resistance=6*mmHg/Qnom)
            annotation (Placement(transformation(extent={{-24,-14},{-4,6}})));
          Physiolibrary.Hydraulic.Components.Resistor Intest_Venules(Resistance=3
                *mmHg/Qnom) "Intestinal venule pressure drop"
            annotation (Placement(transformation(extent={{-54,-12},{-34,8}})));
          Physiolibrary.Hydraulic.Sensors.PressureMeasure Ppv
            annotation (Placement(transformation(extent={{-40,12},{-20,32}})));
          Physiolibrary.Hydraulic.Components.Resistor Intest_Art(Resistance=84*
                mmHg/Qnom) "Intestinal Arteriole resistance"
            annotation (Placement(transformation(extent={{-82,-12},{-62,8}})));
          Physiolibrary.Hydraulic.Sensors.PressureMeasure Pc
            annotation (Placement(transformation(extent={{-62,6},{-42,26}})));
          Physiolibrary.Hydraulic.Sensors.PressureMeasure Phv
            annotation (Placement(transformation(extent={{6,-2},{26,18}})));
          Modelica.Blocks.Sources.RealExpression HVPG_nominal(y=1/(time*mmHg/Qnom))
            annotation (Placement(transformation(extent={{30,2},{10,-18}})));
          replaceable Components.ResistancePressureDep Shunt(enable=true,
              useConductanceInput=false) if useTIPPS constrainedby
            Components.ResistancePressureDep
            annotation (Placement(transformation(extent={{-24,2},{-4,22}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedPump pump(SolutionFlow(
                displayUnit="l/min") = 1.6666666666667e-05)
            annotation (Placement(transformation(extent={{-106,-12},{-86,8}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
            annotation (Placement(transformation(extent={{84,-12},{64,8}})));
        equation
          connect(PRa.pressure,levittCase1SsSiIo. Pra) annotation (Line(points={{52,20},
                  {48,20},{48,18},{44,18}}, color={0,0,127}));
          connect(Liver.q_out, hepatic_vein.q_in) annotation (Line(
              points={{-4,-4},{6,-4},{6,-2},{34,-2}},
              color={0,0,0},
              thickness=1));
          connect(Intest_Venules.q_out, Liver.q_in) annotation (Line(
              points={{-34,-2},{-28,-2},{-28,-4},{-24,-4}},
              color={0,0,0},
              thickness=1));
          connect(Ppv.q_in, Intest_Venules.q_out) annotation (Line(
              points={{-34,16},{-34,-2}},
              color={0,0,0},
              thickness=1));
          connect(levittCase1SsSiIo.Pp,Ppv. pressure)
            annotation (Line(points={{24,18},{-24,18}}, color={0,0,127}));
          connect(Intest_Art.q_out, Intest_Venules.q_in) annotation (Line(
              points={{-62,-2},{-54,-2}},
              color={0,0,0},
              thickness=1));
          connect(Intest_Art.q_out, Pc.q_in) annotation (Line(
              points={{-62,-2},{-58,-2},{-58,10},{-56,10}},
              color={0,0,0},
              thickness=1));
          connect(Liver.q_out,Phv. q_in) annotation (Line(
              points={{-4,-4},{6,-4},{6,-2},{12,-2},{12,2}},
              color={0,0,0},
              thickness=1));
          connect(Phv.pressure,levittCase1SsSiIo. Phv)
            annotation (Line(points={{22,4},{22,10},{24,10}},
                                                       color={0,0,127}));
          connect(Pc.pressure,levittCase1SsSiIo. Pc)
            annotation (Line(points={{-46,12},{-46,26},{24,26}},
                                                        color={0,0,127}));
          connect(hepatic_vein.EP, levittCase1SsSiIo.Pa)
            annotation (Line(points={{34,6},{34,10}}, color={0,0,127}));
          connect(PRa.q_in, hepatic_vein.q_out) annotation (Line(
              points={{62,18},{62,-2},{54,-2}},
              color={0,0,0},
              thickness=1));
          connect(Liver.cond, HVPG_nominal.y) annotation (Line(points={{-14,2},{4,
                  2},{4,-8},{9,-8}}, color={0,0,127}));
          connect(Shunt.q_in, Intest_Venules.q_out) annotation (Line(
              points={{-24,12},{-28,12},{-28,-2},{-34,-2}},
              color={0,0,0},
              thickness=1));
          connect(Shunt.q_out, hepatic_vein.q_in) annotation (Line(
              points={{-4,12},{6,12},{6,-2},{34,-2}},
              color={0,0,0},
              thickness=1));
          connect(CVP.y, hepatic_vein.q_out) annotation (Line(
              points={{64,-2},{54,-2}},
              color={0,0,0},
              thickness=1));
          connect(pump.q_out, Intest_Art.q_in) annotation (Line(
              points={{-86,-2},{-82,-2}},
              color={0,0,0},
              thickness=1));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false), graphics={
                Text(
                  extent={{-112,2},{-82,14}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="1 L/min"),
                Text(
                  extent={{-92,-22},{-54,-10}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="(84 mmHg)"),
                Text(
                  extent={{-62,-22},{-24,-10}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="(3 mmHg)"),
                Text(
                  extent={{-36,-22},{16,-10}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="(4-30 mmHg)"),
                Text(
                  extent={{54,-24},{92,-12}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="5 mmHg")}));
        end ModelSchematicsSmall;

        model HVPG_shuntsSchematics "Model for view not for simulation"
          extends Components.partialAscites(
            redeclare Physiolibrary.Hydraulic.Components.Resistor Liver(
                useConductanceInput=true, Resistance=6*mmHg/Qnom),
            redeclare Physiolibrary.Hydraulic.Components.Resistor IntestinesArt(
                Resistance=84*mmHg/Qnom),
            redeclare Physiolibrary.Hydraulic.Components.Resistor IntestineVenule(
                Resistance=3*mmHg/Qnom));
          Modelica.Blocks.Sources.RealExpression HVPG_nominal(y=1/(time*mmHg/Qnom))
            annotation (Placement(transformation(extent={{-80,4},{-60,24}})));

          parameter Physiolibrary.Types.VolumeFlowRate Qnom=1.666666666666667e-08*(5000*
              0.2) "Nominal flow through the splanchnic circulation";
          parameter Physiolibrary.Types.Pressure mmHg=133.322;
          replaceable Components.ResistancePressureDep collateralShunt(enable=
                true, useConductanceInput=false) if useTIPPS constrainedby
            Physiolibrary.Hydraulic.Interfaces.OnePort
            annotation (Placement(transformation(extent={{0,-28},{20,-8}})));
          parameter Boolean useTIPPS=false;
          Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump(
              SolutionFlow(displayUnit="l/min") = 1.6666666666667e-05)
            annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
            annotation (Placement(transformation(extent={{108,-10},{88,10}})));
        equation
          connect(Liver.cond, HVPG_nominal.y)
            annotation (Line(points={{10,6},{10,14},{-59,14}}, color={0,0,127}));
          connect(collateralShunt.q_in, IntestineVenule.q_out) annotation (Line(
              points={{0,-18},{-14,-18},{-14,0},{-20,0}},
              color={0,0,0},
              thickness=1));
          connect(collateralShunt.q_out, HV.q_in) annotation (Line(
              points={{20,-18},{24,-18},{24,0},{52,0}},
              color={0,0,0},
              thickness=1));
          connect(CVP.y, HV.q_out) annotation (Line(
              points={{88,0},{72,0}},
              color={0,0,0},
              thickness=1));
          connect(unlimitedPump.q_out, IntestinesArt.q_in) annotation (Line(
              points={{-88,0},{-80,0}},
              color={0,0,0},
              thickness=1));
          annotation (Documentation(info="<html>
<p>The TIPS resistance taken from TIPS flow and PVP from Su et al (2012, PMID 22099870).</p>
</html>"),   Diagram(graphics={
                Rectangle(
                  extent={{-52,48},{48,-28}},
                  lineColor={217,67,180},
                  lineThickness=1),
                Text(
                  extent={{-114,4},{-84,16}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="1 L/min"),
                Text(
                  extent={{-90,-20},{-52,-8}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="(84 mmHg)"),
                Text(
                  extent={{-50,-20},{-12,-8}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="(3 mmHg)"),
                Text(
                  extent={{-96,22},{-44,34}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="(4-30 mmHg)"),
                Text(
                  extent={{78,-24},{116,-12}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="5 mmHg")}));
        end HVPG_shuntsSchematics;

        model ModelSchematicsLarge
          Physiolibrary.Hydraulic.Components.Resistor Liver(useConductanceInput=
               false, Resistance=6*mmHg/Qnom)
            annotation (Placement(transformation(extent={{2,-10},{22,10}})));
          AscitesLevitt.LevittCase1SsSiIo AscitesVolume
            annotation (Placement(transformation(extent={{44,18},{64,38}})));
          Physiolibrary.Hydraulic.Components.PumpPressureHead HV(
              useExternalCollapsingPressure=true, p_head0=-266.64477483)
            "Hepatic vein (free)"
            annotation (Placement(transformation(extent={{60,-10},{80,10}})));
          Physiolibrary.Hydraulic.Components.Resistor IntestineVenule(
              Resistance=3*mmHg/Qnom) "Intestinal venule pressure drop"
            annotation (Placement(transformation(extent={{-36,-10},{-16,10}})));
          Physiolibrary.Hydraulic.Components.Resistor IntestinesArt(Resistance=
                84*mmHg/Qnom) "Intestinal Arteriole resistance"
            annotation (Placement(transformation(extent={{-76,-10},{-56,10}})));
          replaceable ResistancePressureDep collateralShunt(enable=true,
              useConductanceInput=false) if useTIPPS constrainedby
            ResistancePressureDep
            annotation (Placement(transformation(extent={{2,-40},{22,-20}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedPump unlimitedPump(
              SolutionFlow(displayUnit="l/min") = 1.6666666666667e-05)
            annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));
          Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
            annotation (Placement(transformation(extent={{108,10},{88,-10}})));
          parameter Physiolibrary.Types.VolumeFlowRate Qnom=
              1.666666666666667e-08*(5000*0.2)
                   "Nominal flow through the splanchnic circulation";
          parameter Physiolibrary.Types.Pressure mmHg=133.322;
        equation
          connect(collateralShunt.q_out, HV.q_in) annotation (Line(
              points={{22,-30},{34,-30},{34,0},{60,0}},
              color={0,0,0},
              thickness=1));
          connect(CVP.y, HV.q_out) annotation (Line(
              points={{88,0},{80,0}},
              color={0,0,0},
              thickness=1));
          connect(unlimitedPump.q_out, IntestinesArt.q_in) annotation (Line(
              points={{-88,0},{-76,0}},
              color={0,0,0},
              thickness=1));
          connect(Liver.q_out,HV. q_in) annotation (Line(
              points={{22,0},{60,0}},
              color={0,0,0},
              thickness=1));
          connect(IntestineVenule.q_out,Liver. q_in) annotation (Line(
              points={{-16,0},{2,0}},
              color={0,0,0},
              thickness=1));
          connect(IntestinesArt.q_out,IntestineVenule. q_in) annotation (Line(
              points={{-56,0},{-36,0}},
              color={0,0,0},
              thickness=1));
          connect(HV.EP, AscitesVolume.Pa) annotation (Line(points={{60,8},{60,
                  20},{54,20}}, color={0,0,127}));
          connect(collateralShunt.q_in,IntestineVenule. q_out) annotation (Line(
              points={{2,-30},{-8,-30},{-8,0},{-16,0}},
              color={0,0,0},
              thickness=1));
          connect(CVP.y,HV. q_out) annotation (Line(
              points={{88,0},{80,0}},
              color={0,0,0},
              thickness=1));
          connect(unlimitedPump.q_out,IntestinesArt. q_in) annotation (Line(
              points={{-88,0},{-76,0}},
              color={0,0,0},
              thickness=1));
          annotation (Diagram(graphics={
                Text(
                  extent={{-114,8},{-84,20}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="1 L/min"),
                Text(
                  extent={{-86,8},{-48,20}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="(84 mmHg)"),
                Text(
                  extent={{-44,8},{-6,20}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="(3 mmHg)"),
                Text(
                  extent={{-12,8},{40,20}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="(4-35 mmHg)"),
                Text(
                  extent={{78,8},{116,20}},
                  textColor={28,108,200},
                  fontSize=14,
                  textString="5 mmHg"),
                Line(
                  points={{46,36},{-42,36},{-42,0},{-42,36}},
                  color={28,108,200},
                  pattern=LinePattern.Dash),
                Line(
                  points={{46,28},{-8,28},{-8,0}},
                  color={28,108,200},
                  pattern=LinePattern.Dash),
                Line(
                  points={{44,20},{28,20},{28,0}},
                  color={28,108,200},
                  pattern=LinePattern.Dash),
                Line(
                  points={{66,28},{86,28},{86,0}},
                  color={28,108,200},
                  pattern=LinePattern.Dash),
                Line(
                  points={{4,-10},{20,10}},
                  color={0,0,0},
                  pattern=LinePattern.Dash,
                  thickness=0.5,
                  arrow={Arrow.None,Arrow.Filled}),
                Text(
                  extent={{-8,30},{40,34}},
                  textColor={28,108,200},
                  textStyle={TextStyle.Italic},
                  textString="Pressure inputs",
                  fontSize=12)}));
        end ModelSchematicsLarge;
      end Schematics;

      type RemodelingModel = enumeration(
          Pd                                                                             "Pressure difference",
          Pt                                "Transmural pressure",
          Shear                                                                                                       "Shear stress",
          Wt                                                          "Wall tension")
            "Remodelling assumption for the shunt";
    end Components;

    package Experiments
      model HVPG_shunts "Evaluation of shunts"
        extends Components.Ascites_Resistance_ShuntsWEso(splenorenalShunt(P_nom=1066.57909932));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump(
            SolutionFlow(displayUnit="l/min") = 1.6666666666667e-05)
          annotation (Placement(transformation(extent={{-120,-70},{-100,-50}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
          annotation (Placement(transformation(extent={{120,-70},{100,-50}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume MAP(P=13332.2387415)
          annotation (Placement(transformation(extent={{-122,-90},{-102,-70}})));
      equation
        connect(HV.q_out, q_out) annotation (Line(
            points={{72,0},{100,0}},
            color={0,0,0},
            thickness=1));
        connect(q_in, IntestinesArt.q_in) annotation (Line(
            points={{-100,0},{-80,0}},
            color={0,0,0},
            thickness=1));
        connect(CVP.y, HV.q_out) annotation (Line(
            points={{100,-60},{80,-60},{80,0},{72,0}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump.q_out, IntestinesArt.q_in) annotation (Line(
            points={{-100,-60},{-88,-60},{-88,0},{-80,0}},
            color={0,0,0},
            thickness=1));
        connect(MAP.y, gastricArt.q_in) annotation (Line(
            points={{-102,-80},{-80,-80}},
            color={0,0,0},
            thickness=1));
        annotation (Documentation(info="<html>
<p>The TIPS resistance taken from TIPS flow and PVP from Su et al (2012, PMID 22099870).</p>
</html>"), experiment(
            StartTime=4,
            StopTime=35,
            __Dymola_Algorithm="Dassl"));
      end HVPG_shunts;

      model HVPGShuntsComparison
        Physiolibrary.Hydraulic.Sources.UnlimitedPump unlimitedPump1(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,100},{-40,120}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump2(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,70},{-40,90}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump5(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
        Components.Ascites_Resistance ascites_NoShunts(liverConductance(y=Lc))
          annotation (Placement(transformation(extent={{-20,100},{0,120}})));
        Components.Ascites_Resistance_Shunts ascites_Shunts(
          splenorenalShunt(
            Comp=7.50062E-09,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom),
          useTIPPS=false,
          liverConductance(y=Lc))
          annotation (Placement(transformation(extent={{-20,40},{0,60}})));
        Components.Ascites_Resistance_Shunts ascites_ShuntStiff(
          splenorenalShunt(
            Comp=3.75031E-09,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom),
          useTIPPS=false,
          liverConductance(y=Lc))
          annotation (Placement(transformation(extent={{-20,70},{0,90}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
          annotation (Placement(transformation(extent={{100,-10},{80,10}})));
        parameter Physiolibrary.Types.VolumeFlowRate Inflow(displayUnit="l/min")
          =1.6666666666667e-05   "Splanchnic perfusion";
        parameter Physiolibrary.Types.HydraulicCompliance Shunt_Compliance=
            7.50062e-08;
        parameter Physiolibrary.Types.Pressure Shunt_Pnom(displayUnit="mmHg")=
          1066.58    "Nominal end-point pressure";
        parameter Physiolibrary.Types.HydraulicResistance Shunt_R0nom(
            displayUnit="(mmHg.min)/l")=7999340000
                         "Nominal resistance at P_nom";
        parameter Physiolibrary.Types.HydraulicResistance R_TIPSS(displayUnit=
              "(mmHg.min)/l")=39996700;
        parameter Physiolibrary.Types.Pressure intestinalPressureDrop(
            displayUnit="mmHg")=11199.08054286
          "Pressure drop at intestinal arteries";
        Physiolibrary.Types.HydraulicConductance Lc=1/((time)*
            ascites_NoShunts.mmHg/ascites_NoShunts.Qnom) "liver conductance";
        parameter Physiolibrary.Types.VolumeFlowRate EmbolizedInflow(
            displayUnit="l/min") = 1e-05
          "Volumetric flow of solution if useSolutionFlowInput=false";
        Components.Ascites_Resistance_Shunts ascites_ShuntsEmb(
          splenorenalShunt(
            Comp=7.50062e-09,
            d=ascites_Shunts.splenorenalShunt.d,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom,
            UsePrescribedDiameter=true),
          useTIPPS=false,
          liverConductance(y=Lc))
          annotation (Placement(transformation(extent={{-26,-18},{-6,2}})));
        Components.Ascites_Resistance_Shunts ascites_ShuntStiffEmb(
          splenorenalShunt(
            Comp=3.75031e-09,
            d=ascites_ShuntStiff.splenorenalShunt.d,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom,
            UsePrescribedDiameter=true),
          useTIPPS=false,
          liverConductance(y=Lc))
          annotation (Placement(transformation(extent={{-26,12},{-6,32}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump3(
            SolutionFlow(displayUnit="l/min") = EmbolizedInflow)
          annotation (Placement(transformation(extent={{-66,12},{-46,32}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump4(
            SolutionFlow(displayUnit="l/min") = EmbolizedInflow)
          annotation (Placement(transformation(extent={{-66,-18},{-46,2}})));
        Components.Ascites_Resistance_Shunts ascites_ShuntStiffEmbSS(
          splenorenalShunt(
            Comp=3.75031E-09,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom),
          useTIPPS=false,
          liverConductance(y=Lc))
          annotation (Placement(transformation(extent={{-26,-52},{-6,-32}})));
        Components.Ascites_Resistance_Shunts ascites_ShuntsEmbSS(
          splenorenalShunt(
            Comp=7.50062E-09,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom),
          useTIPPS=false,
          liverConductance(y=Lc))
          annotation (Placement(transformation(extent={{-26,-82},{-6,-62}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump6(
            SolutionFlow(displayUnit="l/min") = EmbolizedInflow)
          annotation (Placement(transformation(extent={{-64,-54},{-44,-34}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump7(
            SolutionFlow(displayUnit="l/min") = EmbolizedInflow)
          annotation (Placement(transformation(extent={{-64,-84},{-44,-64}})));
        Components.Ascites_Resistance ascites_NoShuntsEmb(liverConductance(y=Lc))
          annotation (Placement(transformation(extent={{-26,-114},{-6,-94}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump8(
            SolutionFlow(displayUnit="l/min") = EmbolizedInflow)
          annotation (Placement(transformation(extent={{-64,-116},{-44,-96}})));
      equation
        connect(ascites_NoShunts.q_out, CVP.y) annotation (Line(
            points={{0,110},{70,110},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_NoShunts.q_in, unlimitedPump1.q_out) annotation (Line(
            points={{-20,110},{-40,110}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Shunts.q_out, CVP.y) annotation (Line(
            points={{0,50},{70,50},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump5.q_out, ascites_Shunts.q_in) annotation (Line(
            points={{-40,50},{-20,50}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump2.q_out, ascites_ShuntStiff.q_in) annotation (Line(
            points={{-40,80},{-20,80}},
            color={0,0,0},
            thickness=1));
        connect(ascites_ShuntStiff.q_out, CVP.y) annotation (Line(
            points={{0,80},{70,80},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump3.q_out, ascites_ShuntStiffEmb.q_in) annotation (
            Line(
            points={{-46,22},{-26,22}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump4.q_out, ascites_ShuntsEmb.q_in) annotation (Line(
            points={{-46,-8},{-26,-8}},
            color={0,0,0},
            thickness=1));
        connect(ascites_ShuntStiffEmb.q_out, CVP.y) annotation (Line(
            points={{-6,22},{70,22},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_ShuntsEmb.q_out, CVP.y) annotation (Line(
            points={{-6,-8},{34,-8},{34,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump6.q_out, ascites_ShuntStiffEmbSS.q_in) annotation
          (Line(
            points={{-44,-44},{-44,-42},{-26,-42}},
            color={0,0,0},
            thickness=1));
        connect(ascites_ShuntsEmbSS.q_in, unlimitedPump7.q_out) annotation (
            Line(
            points={{-26,-72},{-26,-74},{-44,-74}},
            color={0,0,0},
            thickness=1));
        connect(ascites_ShuntsEmbSS.q_out, CVP.y) annotation (Line(
            points={{-6,-72},{32,-72},{32,-68},{72,-68},{72,6},{70,6},{70,0},{
                80,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_ShuntStiffEmbSS.q_out, CVP.y) annotation (Line(
            points={{-6,-42},{72,-42},{72,6},{70,6},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump8.q_out, ascites_NoShuntsEmb.q_in) annotation (
            Line(
            points={{-44,-106},{-44,-104},{-26,-104}},
            color={0,0,0},
            thickness=1));
        connect(ascites_NoShuntsEmb.q_out, CVP.y) annotation (Line(
            points={{-6,-104},{36,-104},{36,-106},{72,-106},{72,6},{70,6},{70,0},
                {80,0}},
            color={0,0,0},
            thickness=1));
        annotation (
          Icon(coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(
            StartTime=5,
            StopTime=30,
            __Dymola_NumberOfIntervals=200,
            Tolerance=1e-06,
            __Dymola_Algorithm="Cvode"),
          __Dymola_Commands(file(ensureSimulated=true) = "TIPPS.mos"
              "Plot HVPG for shunt and non-shunt", file=
                "Fig2 Treatment with TIPSS.mos" "Treatment with TIPSS",
            file="ClinicalPhases-Compliant Stiff.mos" "ClinicalPhases",
            file="ClinicalPhases.mos" "ClinicalPhases"));
      end HVPGShuntsComparison;

      model HVPGShuntsComparison_IncreasingInflow
        extends HVPGShuntsComparison(
          unlimitedPump1(useSolutionFlowInput=true),
          unlimitedPump2(useSolutionFlowInput=true),
          unlimitedPump5(useSolutionFlowInput=true));
        Modelica.Blocks.Sources.Ramp ramp(
          height=(height - 1)*Inflow,
          duration=30.0,
          offset=Inflow,
          startTime=5)
          annotation (Placement(transformation(extent={{-130,70},{-110,90}})));
        parameter Real height=2.5 "Height of ramps";
      equation
        connect(ramp.y, unlimitedPump1.solutionFlow) annotation (Line(points={{
                -109,80},{-76,80},{-76,117},{-50,117}}, color={0,0,127}));
        connect(ramp.y, unlimitedPump2.solutionFlow) annotation (Line(points={{
                -109,80},{-76,80},{-76,96},{-50,96},{-50,87}}, color={0,0,127}));
        connect(ramp.y, unlimitedPump5.solutionFlow) annotation (Line(points={{
                -109,80},{-66,80},{-66,64},{-50,64},{-50,57}}, color={0,0,127}));
        annotation (experiment(
            StartTime=5,
            StopTime=30,
            __Dymola_NumberOfIntervals=2000,
            Tolerance=1e-06,
            __Dymola_Algorithm="Cvode"));
      end HVPGShuntsComparison_IncreasingInflow;

      model HVPGShuntsComparison_extended
        extends HVPGShuntsComparison(R_TIPSS(displayUnit="(mmHg.min)/l"));
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume PA(P=13332.2387415)
          "Arterial pressure"
          annotation (Placement(transformation(extent={{-60,-110},{-40,-90}})));
        Components.Ascites_Resistance_Shunts ascites_Shunts_FixedPressure(
          splenorenalShunt(
            Comp(displayUnit="ml/mmHg") = Shunt_Compliance,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom),
          useTIPPS=false,
          IntestinesArt(Resistance=intestinalPressureDrop/
                ascites_Shunts_FixedPressure.Qnom))
          annotation (Placement(transformation(extent={{-20,-110},{0,-90}})));
        parameter Physiolibrary.Types.VolumeFlowRate Inflow(displayUnit="l/min")=1.6666666666667e-05
                                 "Splanchnic perfusion";
        parameter Physiolibrary.Types.VolumeFlowRate InflowReduced(displayUnit=
              "l/min")=1.3333333333333e-05
                                 "Splanchnic perfusion";
      //  parameter Physiolibrary.Types.HydraulicCompliance Shunt_Compliance(
      //      displayUnit="ml/mmHg")=7.50062e-07;
        parameter Physiolibrary.Types.Pressure Shunt_Pnom(displayUnit="mmHg")=1066.58
                     "Nominal end-point pressure";
        parameter Physiolibrary.Types.HydraulicResistance Shunt_R0nom(displayUnit="(mmHg.min)/l")=
           7999340000    "Nominal resistance at P_nom";
        parameter Physiolibrary.Types.HydraulicResistance R_TIPSS(displayUnit="(mmHg.min)/l")=
           39996700;
        parameter Physiolibrary.Types.Pressure intestinalPressureDrop(displayUnit="mmHg")=
           11199.08054286
          "Pressure drop at intestinal arteries";
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump3(SolutionFlow(
              displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,-50},{-40,-30}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump4(SolutionFlow(
              displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,-80},{-40,-60}})));
        Components.Ascites_Resistance_Shunts ascites_Shunts_TIPS_acute(
            splenorenalShunt(
            Comp=7.50062e-09,
            d=ascites_Shunts.shunt.d,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom,
            UsePrescribedDiameter=true), useTIPPS=true)
          annotation (Placement(transformation(extent={{-20,-80},{0,-60}})));
        Components.Ascites_Resistance_Shunts ascites_ShuntStiff_TIPS_acute(
            splenorenalShunt(
            Comp=3.75031e-09,
            d=ascites_ShuntStiff.shunt.d,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom,
            UsePrescribedDiameter=true), useTIPPS=true)
          annotation (Placement(transformation(extent={{-20,-50},{0,-30}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump6(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,10},{-40,30}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump7(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,-20},{-40,0}})));
        Components.Ascites_Resistance_Shunts ascites_Shunts_TIPS_remodel(useTIPPS=
              true)
          annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
        Components.Ascites_Resistance_Shunts ascites_TIPS(splenorenalShunt(
            Comp(displayUnit="m3/Pa") = 1e-12,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom(displayUnit="(mmHg.min)/l") = 7.9993432449e+16), useTIPPS=
              true)
          annotation (Placement(transformation(extent={{-20,10},{0,30}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume PA1(usePressureInput=
              true, P=13332.2387415)
          "Arterial pressure"
          annotation (Placement(transformation(extent={{-60,-140},{-40,-120}})));
        Components.Ascites_Resistance_Shunts
          ascites_Shunts_TIPS_acute_presssure(splenorenalShunt(
            Comp=7.50062e-09,
            d=ascites_Shunts.shunt.d,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom,
            UsePrescribedDiameter=true), useTIPPS=true)
          annotation (Placement(transformation(extent={{-20,-140},{0,-120}})));
        Modelica.Blocks.Sources.RealExpression ShuntsPressure(y=ascites_Shunts.q_in.pressure)
          annotation (Placement(transformation(extent={{-100,-140},{-80,-120}})));
        Components.Ascites_Resistance_Shunts ascites_TIPS_pressure(
            splenorenalShunt(
            Comp(displayUnit="m3/Pa") = 1e-12,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom(displayUnit="(mmHg.min)/l") = 7.9993432449e+16), useTIPPS=
              true)
          annotation (Placement(transformation(extent={{-20,-170},{0,-150}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume PA2(usePressureInput=
              true, P=13332.2387415)
          "Arterial pressure"
          annotation (Placement(transformation(extent={{-60,-170},{-40,-150}})));
        Modelica.Blocks.Sources.RealExpression NoShuntsPressure(y=
              ascites_NoShunts.q_in.pressure) annotation (Placement(
              transformation(extent={{-100,-170},{-80,-150}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump9(
            SolutionFlow(displayUnit="l/min") = InflowReduced)
          annotation (Placement(transformation(extent={{-60,-200},{-40,-180}})));
        Components.Ascites_Resistance_Shunts ascites_Shunts_spleenEmb(
          splenorenalShunt(
            Comp=7.50062e-09,
            d=ascites_Shunts.shunt.d,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom,
            UsePrescribedDiameter=true),
          useTIPPS=false,
          liverConductance(y=Lc))
          annotation (Placement(transformation(extent={{-20,-200},{0,-180}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump8(
            SolutionFlow(displayUnit="l/min") = InflowReduced)
          annotation (Placement(transformation(extent={{-60,-232},{-40,-212}})));
        Components.Ascites_Resistance_Shunts ascites_Shunts_spleenEmb1(
          splenorenalShunt(
            Comp=7.5006157584566e-13,
            d(displayUnit="mm") = 1e-06,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom,
            UsePrescribedDiameter=true),
          useTIPPS=false,
          liverConductance(y=Lc))
          annotation (Placement(transformation(extent={{-20,-232},{0,-212}})));
      equation
        connect(PA.y, ascites_Shunts_FixedPressure.q_in) annotation (Line(
            points={{-40,-100},{-20,-100}},
            color={0,0,0},
            thickness=1,
            smooth=Smooth.Bezier));
        connect(ascites_Shunts_FixedPressure.q_out, CVP.y) annotation (Line(
            points={{0,-100},{70,-100},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Shunts_TIPS_acute.q_out, CVP.y) annotation (Line(
            points={{0,-70},{70,-70},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump4.q_out, ascites_Shunts_TIPS_acute.q_in)
          annotation (Line(
            points={{-40,-70},{-20,-70}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump3.q_out, ascites_ShuntStiff_TIPS_acute.q_in)
          annotation (Line(
            points={{-40,-40},{-20,-40}},
            color={0,0,0},
            thickness=1));
        connect(ascites_ShuntStiff_TIPS_acute.q_out, CVP.y) annotation (Line(
            points={{0,-40},{70,-40},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Shunts_TIPS_remodel.q_out, CVP.y) annotation (Line(
            points={{0,-10},{70,-10},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump7.q_out, ascites_Shunts_TIPS_remodel.q_in)
          annotation (Line(
            points={{-40,-10},{-20,-10}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump6.q_out, ascites_TIPS.q_in) annotation (Line(
            points={{-40,20},{-20,20}},
            color={0,0,0},
            thickness=1));
        connect(ascites_TIPS.q_out, CVP.y) annotation (Line(
            points={{0,20},{70,20},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(PA1.y, ascites_Shunts_TIPS_acute_presssure.q_in) annotation (
            Line(
            points={{-40,-130},{-20,-130}},
            color={0,0,0},
            thickness=1));
        connect(CVP.y, ascites_Shunts_TIPS_acute_presssure.q_out) annotation (
            Line(
            points={{80,0},{70,0},{70,-130},{0,-130}},
            color={0,0,0},
            thickness=1));
        connect(ShuntsPressure.y, PA1.pressure)
          annotation (Line(points={{-79,-130},{-60,-130}}, color={0,0,127}));
        connect(NoShuntsPressure.y, PA2.pressure)
          annotation (Line(points={{-79,-160},{-60,-160}}, color={0,0,127}));
        connect(PA2.y, ascites_TIPS_pressure.q_in) annotation (Line(
            points={{-40,-160},{-20,-160}},
            color={0,0,0},
            thickness=1));
        connect(ascites_TIPS_pressure.q_out,
          ascites_Shunts_TIPS_acute_presssure.q_out) annotation (Line(
            points={{0,-160},{70,-160},{70,-130},{0,-130}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Shunts_spleenEmb.q_out, CVP.y) annotation (Line(
            points={{0,-190},{70,-190},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump9.q_out, ascites_Shunts_spleenEmb.q_in)
          annotation (Line(
            points={{-40,-190},{-20,-190}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Shunts_spleenEmb1.q_out, CVP.y) annotation (Line(
            points={{0,-222},{70,-222},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump8.q_out, ascites_Shunts_spleenEmb1.q_in)
          annotation (Line(
            points={{-40,-222},{-20,-222}},
            color={0,0,0},
            thickness=1));
        annotation (
          Icon(coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(
            StartTime=4,
            StopTime=35,
            __Dymola_NumberOfIntervals=200,
            Tolerance=1e-06,
            __Dymola_Algorithm="Cvode"),
          __Dymola_Commands(file(ensureSimulated=true) = "TIPPS.mos"
              "Plot HVPG for shunt and non-shunt", file=
                "Fig2 Treatment with TIPSS.mos" "Treatment with TIPSS",
            file="ClinicalPhases-Compliant Stiff.mos" "ClinicalPhases",
            file="ClinicalPhases.mos" "ClinicalPhases"));
      end HVPGShuntsComparison_extended;

      model HVPGShuntsComparison_Size
        extends HVPGShuntsComparison(ascites_NoShunts(Liver(useConductanceInput=
                 true, Resistance(displayUnit="(mmHg.min)/l") = 199983581.1225)));

        Components.HagenPoiseulleConductance hagenPoiseulleConductance(
          d_nominal(displayUnit="mm") = 0.0045,
          mu=5e-3,
          L=0.08,
          useNominalDiameter=true,
          q_nominal=1.6666666666667e-06)
          annotation (Placement(transformation(extent={{-14,-68},{6,-48}})));
        parameter Physiolibrary.Types.HydraulicResistance LiverResistance(displayUnit=
             "(mmHg.min)/l") = 199983581.1225
          "Hydraulic conductance if useConductanceInput=false";
      equation

        annotation (experiment(
            StartTime=5,
            StopTime=35,
            __Dymola_NumberOfIntervals=200,
            Tolerance=1e-06,
            __Dymola_Algorithm="Euler"));
      end HVPGShuntsComparison_Size;

      model CardiovascularSystem_GCG_Asc
        "Cardiovascular part of Guyton-Coleman-Granger's model from 1972"
         extends Modelica.Icons.Example;
         import Physiolibrary.Hydraulic;
        Hydraulic.Components.ElasticVessel pulmonaryVeinsAndLeftAtrium(
          volume_start(displayUnit="l") = 0.0004,
          ZeroPressureVolume(displayUnit="l") = 0.0004,
          Compliance(displayUnit="l/mmHg") = 7.5006157584566e-08)
          annotation (Placement(transformation(extent={{4,74},{24,94}})));
        Hydraulic.Components.ElasticVessel pulmonaryArteries(
          ZeroPressureVolume(displayUnit="l") = 0.00030625,
          Compliance(displayUnit="l/mmHg") = 3.6002955640592e-08,
          volume_start(displayUnit="l") = 0.00038)
          annotation (Placement(transformation(extent={{-62,74},{-42,94}})));
        Hydraulic.Components.Conductor
                 pulmonary(Conductance(displayUnit="l/(mmHg.min)") = 4.1665920538226e-08)
          annotation (Placement(transformation(extent={{-30,74},{-10,94}})));
        Hydraulic.Components.ElasticVessel arteries(
          volume_start(displayUnit="l") = 0.00085,
          ZeroPressureVolume(displayUnit="l") = 0.000495,
          Compliance(displayUnit="l/mmHg") = 2.6627185942521e-08)
          annotation (Placement(transformation(extent={{14,-46},{34,-26}})));
        Hydraulic.Components.ElasticVessel veins(
          Compliance(displayUnit="l/mmHg") = 6.1880080007267e-07,
          volume_start(displayUnit="l") = 0.00325,
          ZeroPressureVolume(displayUnit="l") = 0.00295)
          annotation (Placement(transformation(extent={{-64,-46},{-44,-26}})));
        Hydraulic.Components.Conductor
                 nonMuscle(Conductance(displayUnit="m3/(Pa.s)") = 2.6e-9)
          annotation (Placement(transformation(extent={{-24,-46},{-4,-26}})));
        Hydraulic.Sensors.PressureMeasure pressureMeasure
          annotation (Placement(transformation(extent={{-78,26},{-58,46}})));
        Hydraulic.Components.Pump rightHeart(useSolutionFlowInput=true)
          annotation (Placement(transformation(extent={{-56,8},{-36,28}})));
        Hydraulic.Sensors.PressureMeasure pressureMeasure1
          annotation (Placement(transformation(extent={{-8,26},{12,46}})));
        Hydraulic.Components.Pump leftHeart(useSolutionFlowInput=true)
          annotation (Placement(transformation(extent={{16,6},{36,26}})));
        Physiolibrary.Types.Constants.VolumeFlowRateConst LNormalCO(k(displayUnit="l/min")=
               8.3333333333333e-05)
          annotation (Placement(transformation(extent={{12,42},{20,50}})));
        Hydraulic.Components.Conductor
                 kidney(Conductance(displayUnit="l/(mmHg.min)") = 1.4126159678427e-09)
          annotation (Placement(transformation(extent={{-24,-64},{-4,-44}})));
        Hydraulic.Components.Conductor
                 muscle(Conductance(displayUnit="l/(mmHg.min)") = 1.3001067314658e-09)
          annotation (Placement(transformation(extent={{-24,-28},{-4,-8}})));
        Hydraulic.Components.Conductor
                 largeVeins(Conductance(displayUnit="l/(mmHg.min)") = 1.6888886482791e-07)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-84,-8})));
        Hydraulic.Components.ElasticVessel rightAtrium(
          volume_start(displayUnit="l") = 0.0001,
          ZeroPressureVolume(displayUnit="l") = 0.0001,
          Compliance(displayUnit="l/mmHg") = 3.7503078792283e-08)
          annotation (Placement(transformation(extent={{-82,8},{-62,28}})));
        Physiolibrary.Blocks.Factors.Spline rightStarling(data={{-6,0,0},{-3,0.15,0.104},
              {-1,0.52,0.48},{2,1.96,0.48},{4,2.42,0.123},{8,2.7,0}}, Xscale=101325/760)
          "At filling pressure 0mmHg (because external thorax pressure is -4mmHg) is normal cardiac output (effect=1)."
          annotation (Placement(transformation(extent={{-56,22},{-36,42}})));
        Physiolibrary.Blocks.Factors.Spline leftStarling(data={{-4,0,0},{-1,0.72,0.29},
              {0,1.01,0.29},{3,1.88,0.218333},{10,2.7,0}}, Xscale=101325/760)
          "At filling pressure -0.0029mmHg (because external thorax pressure is -4mmHg) is normal cardiac output (effect=1)."
          annotation (Placement(transformation(extent={{16,22},{36,42}})));
        replaceable Components.Ascites_Resistance ascites_Resistance(
          liverConductance(y=1),
          TIPSS(enable=true, Resistance=39996716.2245),
          Liver(
            enable=true,
            useConductanceInput=false,
            Resistance(displayUnit="(mmHg.min)/l") = 199983581.1225),
          redeclare replaceable Components.ResistanceControlled IntestinesArt(
            Resistance(displayUnit="(mmHg.min)/l") = 639947459.592,
            Q_nominal(displayUnit="l/min") = 1.6666666666667e-05,
            tau=500000.0),
          useTIPPS=false) constrainedby Components.partialAscites
          annotation (Placement(transformation(extent={{-4,-84},{-24,-64}})));
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(table=[0.0,8.3333333333333e-05;
              500,8.3333333333333e-05; 501,5.3333333333333e-05; 1000,5.3333333333333e-05;
              1001,8.3333333333333e-05], extrapolation=Modelica.Blocks.Types.Extrapolation.HoldLastPoint)
          annotation (Placement(transformation(extent={{92,42},{72,62}})));
      equation
        connect(pulmonaryArteries.q_in,pulmonary. q_in) annotation (Line(
            points={{-52,84},{-30,84}},
            thickness=1));
        connect(pulmonary.q_out, pulmonaryVeinsAndLeftAtrium.q_in) annotation (
            Line(
            points={{-10,84},{14,84}},
            thickness=1));
        connect(veins.q_in, nonMuscle.q_in)  annotation (Line(
            points={{-54,-36},{-24,-36}},
            thickness=1));
        connect(nonMuscle.q_out, arteries.q_in)  annotation (Line(
            points={{-4,-36},{24,-36}},
            thickness=1));
        connect(rightHeart.q_out,pulmonaryArteries. q_in) annotation (Line(
            points={{-36,18},{-28,18},{-28,60},{-70,60},{-70,84},{-52,84}},
            thickness=1));
        connect(leftHeart.q_in, pulmonaryVeinsAndLeftAtrium.q_in) annotation (
            Line(
            points={{16,16},{-4,16},{-4,60},{32,60},{32,84},{14,84}},
            thickness=1));
        connect(leftHeart.q_out,arteries. q_in) annotation (Line(
            points={{36,16},{44,16},{44,-36},{24,-36}},
            thickness=1));
        connect(pressureMeasure1.q_in, pulmonaryVeinsAndLeftAtrium.q_in)
          annotation (Line(
            points={{-2,30},{-4,30},{-4,60},{32,60},{32,84},{14,84}},
            thickness=1));
        connect(muscle.q_out, arteries.q_in) annotation (Line(
            points={{-4,-18},{10,-18},{10,-36},{24,-36}},
            thickness=1));
        connect(kidney.q_out, arteries.q_in) annotation (Line(
            points={{-4,-54},{10,-54},{10,-36},{24,-36}},
            thickness=1));
        connect(kidney.q_in, nonMuscle.q_in) annotation (Line(
            points={{-24,-54},{-34,-54},{-34,-36},{-24,-36}},
            thickness=1));
        connect(muscle.q_in, nonMuscle.q_in) annotation (Line(
            points={{-24,-18},{-34,-18},{-34,-36},{-24,-36}},
            thickness=1));
        connect(veins.q_in, largeVeins.q_out) annotation (Line(
            points={{-54,-36},{-84,-36},{-84,-18}},
            thickness=1));
        connect(largeVeins.q_in, rightAtrium.q_in) annotation (Line(
            points={{-84,2},{-84,18},{-72,18}},
            thickness=1));
        connect(rightAtrium.q_in, rightHeart.q_in) annotation (Line(
            points={{-72,18},{-56,18}},
            thickness=1));
        connect(rightHeart.solutionFlow, rightStarling.y) annotation (Line(
            points={{-46,25},{-46,28},{-46,28}},
            color={0,0,127}));
        connect(leftStarling.y, leftHeart.solutionFlow) annotation (Line(
            points={{26,28},{26,23}},
            color={0,0,127}));
        connect(pressureMeasure.pressure, rightStarling.u) annotation (Line(
            points={{-62,32},{-54,32}},
            color={0,0,127}));
        connect(pressureMeasure1.pressure, leftStarling.u) annotation (Line(
            points={{8,32},{18,32}},
            color={0,0,127}));
        connect(pressureMeasure.q_in, rightAtrium.q_in) annotation (Line(
            points={{-72,30},{-72,18}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Resistance.q_in, arteries.q_in) annotation (Line(
            points={{-4,-74},{10,-74},{10,-36},{24,-36}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Resistance.q_out, nonMuscle.q_in) annotation (Line(
            points={{-24,-74},{-34,-74},{-34,-36},{-24,-36}},
            color={0,0,0},
            thickness=1));
        connect(combiTimeTable.y[1], leftStarling.yBase)
          annotation (Line(points={{71,52},{26,52},{26,34}}, color={0,0,127}));
        connect(combiTimeTable.y[1], rightStarling.yBase)
          annotation (Line(points={{71,52},{-46,52},{-46,34}}, color={0,0,127}));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics={Text(
                extent={{-82,-80},{80,-100}},
                lineColor={175,175,175},
                textString=
                    "Circulation part of Guyton-Coleman-Granger's model from 1972")}),
                                                Documentation(info="<html>
<p>Cardiovascular subsystem in famous Guyton-Coleman-Granger model from 1972. </p>
<p><br/>Model, all parameters and all initial values are from article: </p>
<p>A.C. Guyton, T.G. Coleman, H.J. Granger (1972). &quot;Circulation: overall regulation.&quot; Annual review of physiology 34(1): 13-44.</p>
</html>",   revisions="<html>
<p><i>2014</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(
            StopTime=1500,
            Tolerance=1e-06,
            __Dymola_Algorithm="Cvode"));
      end CardiovascularSystem_GCG_Asc;

      model CardiovascularSystem_GCG_Asc_NoReg
        extends CardiovascularSystem_GCG_Asc(ascites_Resistance(redeclare
              Physiolibrary.Hydraulic.Components.Resistor IntestinesArt(
                Resistance=84*ascites_Resistance.mmHg/ascites_Resistance.Qnom)));
      end CardiovascularSystem_GCG_Asc_NoReg;

      model CardiovascularSystem_GCG_Asc_NoReg_NoCOreg
        extends CardiovascularSystem_GCG_Asc_NoReg(
          rightStarling(enabled=false),
          leftStarling(enabled=false),
          combiTimeTable(extrapolation=Modelica.Blocks.Types.Extrapolation.HoldLastPoint, table=[
                0.0,8.3333333333333E-05; 500.0,8.3333333333333E-05; 501.0,6.66666E-05;
                1000.0,6.66666E-05; 1001.0,8.3333333333333E-05]),
          ascites_Resistance(
            Liver(Resistance(displayUnit="(mmHg.min)/l") = 199984000.0)));
        annotation (__Dymola_Commands(file="Reducing CO.mos" "Reducing CO"));
      end CardiovascularSystem_GCG_Asc_NoReg_NoCOreg;

      model CardiovascularSystem_GCG_Asc_NoReg_NoCOreg_TIPSS
        extends CardiovascularSystem_GCG_Asc_NoReg_NoCOreg(redeclare
            Components.Ascites_Resistance_Shunts ascites_Resistance(
            useTIPPS=true,
            Liver(useConductanceInput=false),
            TIPSS(Resistance=39996716.2245),
            liverConductance(y=1),
            splenorenalShunt(
              Comp(displayUnit="ml/mmHg") = 1.5001231516913e-07,
              R_nom(displayUnit="(mmHg.min)/l"),
              P_nom=1466.546261565)));
      end CardiovascularSystem_GCG_Asc_NoReg_NoCOreg_TIPSS;

      model HVPGCasesComparison
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
          annotation (Placement(transformation(extent={{100,-10},{80,10}})));
        Components.Ascites_Resistance_Shunts ascites_Shunts(splenorenalShunt(
            Comp(displayUnit="ml/mmHg") = Shunt_Compliance,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom), useTIPPS=false)
          annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump2(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
        parameter Physiolibrary.Types.VolumeFlowRate Inflow(displayUnit="l/min")=
           1.6666666666667e-05   "Splanchnic perfusion";
        parameter Physiolibrary.Types.HydraulicCompliance Shunt_Compliance(
            displayUnit="ml/mmHg")=1.50012e-07;
        parameter Physiolibrary.Types.Pressure Shunt_Pnom(displayUnit="mmHg")=
          1466.55    "Nominal end-point pressure";
        parameter Physiolibrary.Types.HydraulicResistance Shunt_R0nom=
            7999343244.9 "Nominal resistance at P_nom";
        parameter Physiolibrary.Types.HydraulicResistance R_TIPSS(displayUnit=
              "(mmHg.min)/l") = 39996700.0;
        Components.Ascites_Resistance_Shunts ascites_ShuntsStiff(
            splenorenalShunt(
            Comp(displayUnit="ml/mmHg") = 7.5006157584566e-09,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom), useTIPPS=false)
          annotation (Placement(transformation(extent={{-20,-40},{0,-20}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump1(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,-40},{-40,-20}})));
      equation
        connect(ascites_Shunts.q_out, CVP.y) annotation (Line(
            points={{0,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Shunts.q_in, unlimitedPump2.q_out) annotation (Line(
            points={{-20,0},{-40,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_ShuntsStiff.q_in, unlimitedPump1.q_out) annotation (
            Line(
            points={{-20,-30},{-40,-30}},
            color={0,0,0},
            thickness=1));
        connect(ascites_ShuntsStiff.q_out, CVP.y) annotation (Line(
            points={{0,-30},{70,-30},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        annotation (
          Icon(coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(
            StartTime=4,
            StopTime=35,
            __Dymola_NumberOfIntervals=200,
            Tolerance=1e-06,
            __Dymola_Algorithm="Cvode"),
          __Dymola_Commands(file(ensureSimulated=true) = "TIPPS.mos"
              "Plot HVPG for shunt and non-shunt", file=
                "Fig2 Treatment with TIPSS.mos" "Treatment with TIPSS"));
      end HVPGCasesComparison;

      model AscitesPhases
        Components.Ascites_Resistance_Shunts ascites_Shunts_stiff(
          splenorenalShunt(
            Comp(displayUnit="ml/mmHg") = 1.50012E-08,
            P_nom(displayUnit="mmHg") = 1466.55,
            R_nom=7999343244.9),
          useTIPPS=false,
          levittCase1SsSiIo(D=D))
          annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump2(SolutionFlow(
              displayUnit="l/min") = 1.6666666666667e-05)
          annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
          annotation (Placement(transformation(extent={{80,-10},{60,10}})));

        Components.Ascites_Resistance_Shunts ascites_Shunts_compliant(
          splenorenalShunt(
            Comp(displayUnit="ml/mmHg") = 7.5006157584566e-08,
            P_nom(displayUnit="mmHg") = 1466.55,
            R_nom=7999343244.9),
          useTIPPS=false,
          levittCase1SsSiIo(D=D))
          annotation (Placement(transformation(extent={{-40,-42},{-20,-22}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump1(SolutionFlow(
              displayUnit="l/min") = 1.6666666666667e-05)
          annotation (Placement(transformation(extent={{-80,-42},{-60,-22}})));
        parameter Physiolibrary.Types.HydraulicCompliance D=9.00074E-06;
        Components.Ascites_Resistance_Shunts ascites_NoShunt(
          splenorenalShunt(
            Comp(displayUnit="l/mmHg") = 7.5006157584566e-15,
            P_nom(displayUnit="mmHg") = 1466.55,
            R_nom=7.9993432449e+18),
          useTIPPS=false,
          levittCase1SsSiIo(D=D))
          annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump3(
            SolutionFlow(displayUnit="l/min") = 1.6666666666667e-05)
          annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
      equation

        connect(ascites_Shunts_stiff.q_out, CVP.y) annotation (Line(
            points={{-20,0},{60,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Shunts_stiff.q_in, unlimitedPump2.q_out) annotation (Line(
            points={{-40,0},{-60,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Shunts_compliant.q_out, CVP.y) annotation (Line(
            points={{-20,-32},{20,-32},{20,0},{60,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Shunts_compliant.q_in, unlimitedPump1.q_out) annotation (Line(
            points={{-40,-32},{-60,-32}},
            color={0,0,0},
            thickness=1));
        connect(ascites_NoShunt.q_out, CVP.y) annotation (Line(
            points={{-20,30},{20,30},{20,0},{60,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_NoShunt.q_in, unlimitedPump3.q_out) annotation (Line(
            points={{-40,30},{-60,30}},
            color={0,0,0},
            thickness=1));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end AscitesPhases;

      model HVPGShuntsComparison_Shear
        extends HVPGShuntsComparison(ascites_Shunts(splenorenalShunt(
                                                          rm=Lymphatics.Hemodynamics.Components.RemodelingModel.Pt)));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump3(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,-60},{-40,-40}})));
        Components.Ascites_Resistance_Shunts ascites_ShuntStiff_dp(
            splenorenalShunt(
            Comp(displayUnit="ml/mmHg") = 7.50062e-10,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom,
            rm=Lymphatics.Hemodynamics.Components.RemodelingModel.Pd), useTIPPS=
             false)
          annotation (Placement(transformation(extent={{-20,-60},{0,-40}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump4(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,-100},{-40,-80}})));
        Components.Ascites_Resistance_Shunts ascites_Shunt_shear(
            splenorenalShunt(
            Comp(displayUnit="ml/mmHg") = 7.50062e-10,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom,
            rm=Lymphatics.Hemodynamics.Components.RemodelingModel.Shear),
            useTIPPS=false)
          annotation (Placement(transformation(extent={{-20,-100},{0,-80}})));
      equation
        connect(unlimitedPump3.q_out, ascites_ShuntStiff_dp.q_in) annotation (
            Line(
            points={{-40,-50},{-20,-50}},
            color={0,0,0},
            thickness=1));
        connect(ascites_ShuntStiff_dp.q_out, CVP.y) annotation (Line(
            points={{0,-50},{70,-50},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump4.q_out, ascites_Shunt_shear.q_in) annotation (
            Line(
            points={{-40,-90},{-20,-90}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Shunt_shear.q_out, CVP.y) annotation (Line(
            points={{0,-90},{70,-90},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
      end HVPGShuntsComparison_Shear;

      model TipsFlow
        Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium =
              Medium,
          p=Medium.p_default,
          nPorts=3)             annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={70,30})));
        inner Modelica.Fluid.System system
          annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
        Modelica.Fluid.Pipes.StaticPipe pipe(
          redeclare package Medium = Medium,
          length(displayUnit="cm") = L,
          diameter=diameter,
          redeclare model FlowModel =
              Modelica.Fluid.Pipes.BaseClasses.FlowModels.DetailedPipeFlow (show_Res=true))
          annotation (Placement(transformation(extent={{-18,20},{2,40}})));
        replaceable package Medium =
            Modelica.Media.Incompressible.Examples.Glycol47
          constrainedby Modelica.Media.Interfaces.PartialMedium annotation (
            choicesAllMatching=true);
        Modelica.Fluid.Sources.MassFlowSource_T boundary2(
          redeclare package Medium = Medium,
          use_m_flow_in=false,
          m_flow(displayUnit="kg/min") = 0.016666666666667,
          nPorts=1) annotation (Placement(transformation(extent={{-80,20},{-60,
                  40}})));

      Physiolibrary.Types.Pressure dp = pipe.port_a.p - pipe.port_b.p;
      Real dp_pipe_mmHg = (pipe.port_a.p - pipe.port_b.p)/133.322;
      Real dp_cond_mmHg = conductor.dp/133.322;
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump7(SolutionFlow(
              displayUnit="l/min") = 1.6666666666667e-05)
          annotation (Placement(transformation(extent={{-78,-50},{-58,-30}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=Medium.p_default)
          annotation (Placement(transformation(extent={{82,-40},{62,-20}})));
        Components.HagenPoiseulleConductance hagenPoiseulleConductance(
          d_nominal(displayUnit="mm") = diameter,
          mu=5e-3,
          L=L,
          useNominalDiameter=true,
          q_nominal=1.6666666666667e-06)
          annotation (Placement(transformation(extent={{-34,-76},{-14,-56}})));
        Physiolibrary.Hydraulic.Components.Conductor conductor(useConductanceInput=true)
          annotation (Placement(transformation(extent={{-32,-42},{-12,-22}})));
        parameter Modelica.Units.SI.Diameter diameter=d_mm/1000
          "Diameter of circular pipe";
          parameter Real d_mm = 1;

        Modelica.Fluid.Sources.MassFlowSource_T boundary3(
          redeclare package Medium = Medium,
          use_m_flow_in=false,
          m_flow(displayUnit="kg/min") = 0.016666666666667,
          nPorts=1) annotation (Placement(transformation(extent={{-80,48},{-60,
                  68}})));
        Modelica.Fluid.Fittings.Orifices.ThickEdgedOrifice thickEdgedOrifice1(
            redeclare package Medium = Medium, geometry=
              Modelica.Fluid.Fittings.BaseClasses.Orifices.ThickEdgedOrifice.Choices.circular(
                    diameter*2,
                    diameter,
                    L))
          annotation (Placement(transformation(extent={{-20,48},{0,68}})));
        Modelica.Fluid.Sources.MassFlowSource_T boundary4(
          redeclare package Medium = Medium,
          use_m_flow_in=false,
          m_flow(displayUnit="kg/min") = 0.016666666666667,
          nPorts=1) annotation (Placement(transformation(extent={{-80,-10},{-60,
                  10}})));
        Modelica.Fluid.Fittings.Bends.CurvedBend curvedBend(redeclare package
            Medium = Medium, geometry(
            d_hyd=diameter,
            R_0(displayUnit="m") = sqrt(360/60*L/3.14),
            delta(displayUnit="deg") = 1.0471975511966))
          annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
        parameter Physiolibrary.Types.Length L=0.01;
      equation
        connect(unlimitedPump7.q_out, conductor.q_in) annotation (Line(
            points={{-58,-40},{-38,-40},{-38,-32},{-32,-32}},
            color={0,0,0},
            thickness=1));
        connect(CVP.y, conductor.q_out) annotation (Line(
            points={{62,-30},{60,-30},{60,-32},{-12,-32}},
            color={0,0,0},
            thickness=1));
        connect(hagenPoiseulleConductance.hydraulicconductance, conductor.cond)
          annotation (Line(points={{-13.8,-66},{-4,-66},{-4,-64},{10,-64},{10,-18},{-22,
                -18},{-22,-26}}, color={0,0,127}));
        connect(boundary3.ports[1], thickEdgedOrifice1.port_a)
          annotation (Line(points={{-60,58},{-20,58}}, color={0,127,255}));
        connect(boundary2.ports[1], pipe.port_a)
          annotation (Line(points={{-60,30},{-18,30}}, color={0,127,255}));
        connect(thickEdgedOrifice1.port_b, boundary1.ports[1]) annotation (Line(
              points={{0,58},{50,58},{50,31.3333},{60,31.3333}}, color={0,127,
                255}));
        connect(pipe.port_b, boundary1.ports[2]) annotation (Line(points={{2,30},
                {32,30},{32,30},{60,30}}, color={0,127,255}));
        connect(boundary4.ports[1], curvedBend.port_a)
          annotation (Line(points={{-60,0},{-20,0}}, color={0,127,255}));
        connect(boundary1.ports[3], curvedBend.port_b) annotation (Line(points={{60,
                28.6667},{8,28.6667},{8,0},{0,0}},      color={0,127,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TipsFlow;

      model HVPGShuntsComparison_transmural
        "Comparison remodeling sensitivity to dP and P_transmural"
        extends HVPGShuntsComparison;
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump3(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,-50},{-40,-30}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump4(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,-80},{-40,-60}})));
        Components.Ascites_Resistance_Shunts ascites_Shunts_transm(
            splenorenalShunt(
            Comp=7.50062e-09,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom,
            useExternalCollapsingPressure=true,
            rm=Lymphatics.Hemodynamics.Components.RemodelingModel.Pt), useTIPPS=
             false)
          annotation (Placement(transformation(extent={{-20,-80},{0,-60}})));
        Components.Ascites_Resistance_Shunts ascites_ShuntStiff_transm(
            splenorenalShunt(
            useConductanceInput=true,
            Comp=3.75031e-09,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom,
            useExternalCollapsingPressure=true,
            rm=Lymphatics.Hemodynamics.Components.RemodelingModel.Pt), useTIPPS=
             false)
          annotation (Placement(transformation(extent={{-20,-50},{0,-30}})));
      equation
        connect(ascites_Shunts_transm.q_out, CVP.y) annotation (Line(
            points={{0,-70},{70,-70},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump4.q_out, ascites_Shunts_transm.q_in) annotation (
            Line(
            points={{-40,-70},{-20,-70}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump3.q_out, ascites_ShuntStiff_transm.q_in)
          annotation (Line(
            points={{-40,-40},{-20,-40}},
            color={0,0,0},
            thickness=1));
        connect(ascites_ShuntStiff_transm.q_out, CVP.y) annotation (Line(
            points={{0,-40},{70,-40},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
      end HVPGShuntsComparison_transmural;

      model HVPGShuntsForSimulator
        Physiolibrary.Hydraulic.Sources.UnlimitedPump unlimitedPump1(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,100},{-40,120}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump2(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,70},{-40,90}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump   unlimitedPump5(
            SolutionFlow(displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
        Components.Ascites_Resistance ascites_NoShunts(liverConductance(y=Lc))
          annotation (Placement(transformation(extent={{-20,100},{0,120}})));
        Components.Ascites_Resistance_Shunts ascites_Shunts(
          splenorenalShunt(
            Comp=7.50062E-09,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom),
          useTIPPS=true,
          liverConductance(y=Lc),
          TIPSS(Resistance=TipsResistance))
          annotation (Placement(transformation(extent={{-20,40},{0,60}})));
        Components.Ascites_Resistance_Shunts ascites_ShuntDefault(
          splenorenalShunt(
            Comp=7.50062E-09,
            P_nom(displayUnit="mmHg") = Shunt_Pnom,
            R_nom=Shunt_R0nom),
          useTIPPS=false,
          liverConductance(y=Lc),
          TIPSS(Resistance(displayUnit="(mmHg.min)/l")))
          annotation (Placement(transformation(extent={{-20,70},{0,90}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
          annotation (Placement(transformation(extent={{100,-10},{80,10}})));
        parameter Physiolibrary.Types.VolumeFlowRate Inflow(displayUnit="l/min")=
           1.6666666666667e-05   "Splanchnic perfusion";
        parameter Physiolibrary.Types.HydraulicCompliance Shunt_Compliance=7.50062E-08;
        parameter Physiolibrary.Types.Pressure Shunt_Pnom(displayUnit="mmHg")=
          1066.58    "Nominal end-point pressure";
        parameter Physiolibrary.Types.HydraulicResistance Shunt_R0nom(
            displayUnit="(mmHg.min)/l")=7999340000
                         "Nominal resistance at P_nom";
        parameter Physiolibrary.Types.HydraulicResistance R_TIPSS(displayUnit=
              "(mmHg.min)/l")=39996700;
        parameter Physiolibrary.Types.Pressure intestinalPressureDrop(
            displayUnit="mmHg")=11199.08054286
          "Pressure drop at intestinal arteries";
        Physiolibrary.Types.HydraulicConductance Lc=1/((time + 1e-3)*
            ascites_NoShunts.mmHg/ascites_NoShunts.Qnom) "liver conductance";
        parameter Physiolibrary.Types.HydraulicResistance TipsResistance=if TipsOn > 0.5 then 7999343.2449*
            (15/1.5) else 1e12 "Hydraulic conductance if useConductanceInput=false";
        parameter Real TipsOn=0;
      equation
        connect(ascites_NoShunts.q_out, CVP.y) annotation (Line(
            points={{0,110},{70,110},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_NoShunts.q_in, unlimitedPump1.q_out) annotation (Line(
            points={{-20,110},{-40,110}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Shunts.q_out, CVP.y) annotation (Line(
            points={{0,50},{70,50},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump5.q_out, ascites_Shunts.q_in) annotation (Line(
            points={{-40,50},{-20,50}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump2.q_out, ascites_ShuntDefault.q_in) annotation (Line(
            points={{-40,80},{-20,80}},
            color={0,0,0},
            thickness=1));
        connect(ascites_ShuntDefault.q_out, CVP.y) annotation (Line(
            points={{0,80},{70,80},{70,0},{80,0}},
            color={0,0,0},
            thickness=1));
        annotation (
          Icon(coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(
            StopTime=30,
            __Dymola_NumberOfIntervals=200,
            Tolerance=1e-06,
            __Dymola_Algorithm="Cvode"),
          __Dymola_Commands(file(ensureSimulated=true) = "TIPPS.mos"
              "Plot HVPG for shunt and non-shunt", file=
                "Fig2 Treatment with TIPSS.mos" "Treatment with TIPSS",
            file="ClinicalPhases-Compliant Stiff.mos" "ClinicalPhases",
            file="ClinicalPhases.mos" "ClinicalPhases"));
      end HVPGShuntsForSimulator;

      model DynamicHVPG
        Physiolibrary.Hydraulic.Sources.UnlimitedPump unlimitedPump1(SolutionFlow(
              displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-50,20},{-30,40}})));
        Components.Ascites_Resistance ascites_NoShunts(liverConductance(y=Lc))
          annotation (Placement(transformation(extent={{-10,20},{10,40}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
          annotation (Placement(transformation(extent={{110,-30},{90,-10}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedPump unlimitedPump2(SolutionFlow(
              displayUnit="l/min") = Inflow)
          annotation (Placement(transformation(extent={{-50,-12},{-30,8}})));
        Components.Ascites_ResistanceDynamic ascites_dynamic(liverConductance(y=
               Lc))
          annotation (Placement(transformation(extent={{-10,-12},{10,8}})));

      parameter Physiolibrary.Types.HydraulicConductance Lc=1/((15)*
            ascites_NoShunts.mmHg/ascites_NoShunts.Qnom) "liver conductance";
        parameter Physiolibrary.Types.VolumeFlowRate Inflow=1.6666666666667e-05
          "Volumetric flow of solution";
      equation
        connect(ascites_NoShunts.q_out,CVP. y) annotation (Line(
            points={{10,30},{80,30},{80,-20},{90,-20}},
            color={0,0,0},
            thickness=1));
        connect(ascites_NoShunts.q_in,unlimitedPump1. q_out) annotation (Line(
            points={{-10,30},{-30,30}},
            color={0,0,0},
            thickness=1));
        connect(unlimitedPump2.q_out, ascites_dynamic.q_in) annotation (Line(
            points={{-30,-2},{-10,-2}},
            color={0,0,0},
            thickness=1));
        connect(ascites_dynamic.q_out, CVP.y) annotation (Line(
            points={{10,-2},{80,-2},{80,-20},{90,-20}},
            color={0,0,0},
            thickness=1));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=1728000, __Dymola_Algorithm="Dassl"));
      end DynamicHVPG;

      model HVPG_shunts_EsophagealParams
        extends HVPG_shunts(
          MAP(P(displayUnit="mmHg") = 13332.2387415),
          splenorenalShunt(Comp(displayUnit="ml/mmHg") = 3.75031E-09),
          gastricArt(Resistance(displayUnit="(mmHg.min)/l") = 5599540000.0),
          EsophagealVeins(
            rm=Lymphatics.Hemodynamics.Components.RemodelingModel.Pd,
            UsePrescribedDiameter=false,
            Comp(displayUnit="ml/mmHg") = 7.50062e-11,
            P_nom(displayUnit="mmHg") = 1999.84,
            R_nom(displayUnit="(mmHg.min)/l") = 79993400000));
      end HVPG_shunts_EsophagealParams;

      model HVPG_shunts_EsophagealParams_TestParams
        extends HVPG_shunts_EsophagealParams(EsophagealVeins(Comp(displayUnit=
                  "ml/mmHg") = 3.75031E-10), splenorenalShunt(Comp=7.50062E-11,
              L=0.05));
      end HVPG_shunts_EsophagealParams_TestParams;
    end Experiments;
  end Hemodynamics;

  package Deprecated
    model AbdominalCompliance
      extends Modelica.Blocks.Interfaces.SISO(u(unit="m3"), y(unit="Pa"));

      parameter Physiolibrary.Types.Volume V0 "Nominal volume";
      parameter Physiolibrary.Types.Pressure P0 "Nominal pressure";
      parameter Physiolibrary.Types.HydraulicCompliance D = 800 "Compliance";
    equation
      u = (y - P0)*D + V0;

    end AbdominalCompliance;

    model Ascites
      parameter Physiolibrary.Types.OsmoticPermeability L_Y=1.250102626409427e-10
          *(7.86/60);
          parameter Physiolibrary.Types.Pressure P_min=266.64477483
                                                       "nominal abodminal cavity pressure";
      Physiolibrary.Hydraulic.Sources.UnlimitedPump SplanchnicInflow(SolutionFlow=8.3333333333333e-06)
        annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));
      Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=266.64477483)
        annotation (Placement(transformation(extent={{100,-100},{80,-80}})));
      Physiolibrary.Hydraulic.Components.PumpPressureHead HV(p_head0=
            266.64477483) "Hepatic vein pressure gradient"
        annotation (Placement(transformation(extent={{60,-100},{40,-80}})));
      Physiolibrary.Hydraulic.Components.PumpPressureHead HPVG(p_head0=
            533.28954966) "Hepatic vein pressure gradient"
        annotation (Placement(transformation(extent={{20,-100},{0,-80}})));
      Physiolibrary.Hydraulic.Sensors.PressureMeasure P_HV
        "Pressure in hepatic vein"
        annotation (Placement(transformation(extent={{34,-78},{54,-58}})));
      Physiolibrary.Hydraulic.Sensors.PressureMeasure P_PV "portal vein pressure"
        annotation (Placement(transformation(extent={{-8,-78},{12,-58}})));
      Physiolibrary.Hydraulic.Components.PumpPressureHead P_IC(p_head0=
            399.967162245) "Intestinal capillaries pressure drop"
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
      Deprecated.AbdominalCompliance abdominalCompliance(
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
  end Deprecated;
  annotation (uses(Physiolibrary(version="2.4.1"), Modelica(version="4.0.0"),
      ADAN_main(version="1.1")));
end Lymphatics;
