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
  annotation (uses(Physiolibrary(version="2.4.1"), Modelica(version="4.0.0")));
end Lymphatics;
