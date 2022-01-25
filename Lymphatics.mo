within ;
package Lymphatics

  package OncoticPressures
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
  end OncoticPressures;

  package Hemodynamics
    model CVS_GCG
      extends Physiolibrary.Hydraulic.Examples.CardiovascularSystem_GCG(
          nonMuscle(Conductance(displayUnit="m3/(Pa.s)") = 2.6e-9));
      Components.Ascites_Resistance ascites_Resistance(
        TIPSS(enable=true, Resistance=Modelica.Constants.inf),
        Liver(enable=true, useConductanceInput=true),
        realExpression1(y=1/(max(time/100, 5)*ascites_Resistance.mmHg/
              ascites_Resistance.Qnom)))
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
      extends CVS_GCG(ascites_Resistance(TIPSS(Resistance=7999343.2449*(15/1.5))));
      annotation (experiment(StopTime=3000, __Dymola_Algorithm="Dassl"));
    end CVS_GCG_TIPS;

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

    package Tests
      model partialAscites
        "Levitts ascites with nominal hemodynamics break down"
        AscitesLevitt.LevittCase1SsSiIo levittCase1SsSiIo
          annotation (Placement(transformation(extent={{42,54},{62,74}})));
        Physiolibrary.Hydraulic.Components.PumpPressureHead
                          HV(useExternalCollapsingPressure=true, p_head0=-266.64477483)
          "Hepatic vein (free)"
          annotation (Placement(transformation(extent={{40,-10},{60,10}})));
        Physiolibrary.Hydraulic.Sensors.PressureMeasure PRa "Right atrial pressure"
          annotation (Placement(transformation(extent={{86,58},{66,78}})));
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
          annotation (Placement(transformation(extent={{-20,58},{0,78}})));
        replaceable
        Physiolibrary.Hydraulic.Components.PumpPressureHead IntestinesArt(p_head0=-2666.4477483)
          constrainedby Physiolibrary.Hydraulic.Interfaces.OnePort
          "Intestinal Arteriole resistance"
          annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
        Physiolibrary.Hydraulic.Sensors.PressureMeasure Pc
          annotation (Placement(transformation(extent={{-48,66},{-28,86}})));
        Physiolibrary.Hydraulic.Sensors.PressureMeasure Phv
          annotation (Placement(transformation(extent={{18,50},{38,70}})));
      equation
        connect(PRa.pressure, levittCase1SsSiIo.Pra) annotation (Line(points={{70,64},
                {62,64}},                 color={0,0,127}));
        connect(Liver.q_out, HV.q_in) annotation (Line(
            points={{20,0},{40,0}},
            color={0,0,0},
            thickness=1));
        connect(IntestineVenule.q_out, Liver.q_in) annotation (Line(
            points={{-20,0},{0,0}},
            color={0,0,0},
            thickness=1));
        connect(Ppv.q_in, IntestineVenule.q_out) annotation (Line(
            points={{-14,62},{-14,0},{-20,0}},
            color={0,0,0},
            thickness=1));
        connect(levittCase1SsSiIo.Pp, Ppv.pressure)
          annotation (Line(points={{42,64},{-4,64}},  color={0,0,127}));
        connect(IntestinesArt.q_out, IntestineVenule.q_in) annotation (Line(
            points={{-60,0},{-40,0}},
            color={0,0,0},
            thickness=1));
        connect(IntestinesArt.q_out, Pc.q_in) annotation (Line(
            points={{-60,0},{-42,0},{-42,70}},
            color={0,0,0},
            thickness=1));
        connect(Liver.q_out,Phv. q_in) annotation (Line(
            points={{20,0},{24,0},{24,54}},
            color={0,0,0},
            thickness=1));
        connect(Phv.pressure, levittCase1SsSiIo.Phv)
          annotation (Line(points={{34,56},{42,56}}, color={0,0,127}));
        connect(Pc.pressure, levittCase1SsSiIo.Pc)
          annotation (Line(points={{-32,72},{42,72}}, color={0,0,127}));
        connect(HV.EP, levittCase1SsSiIo.Pa) annotation (Line(points={{40,8},{
                40,40},{52,40},{52,56}}, color={0,0,127}));
        connect(PRa.q_in, HV.q_out) annotation (Line(
            points={{80,62},{80,0},{60,0}},
            color={0,0,0},
            thickness=1));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          experiment(
            StartTime=6,
            StopTime=25,
            __Dymola_Algorithm="Dassl"));
      end partialAscites;

      model Ascites_PG
        extends Tests.partialAscites(Liver(useExternalPressureHead=true));
        Modelica.Blocks.Sources.RealExpression realExpression(y=-time*133.322)
          annotation (Placement(transformation(extent={{-52,-52},{-32,-32}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(P=666.611937075)
          annotation (Placement(transformation(extent={{100,-10},{80,10}})));
      equation
        connect(Liver.p_headIn,realExpression. y) annotation (Line(points={{10,10},
                {-24,10},{-24,-42},{-31,-42}},
                                         color={0,0,127}));
        connect(CVP.y, HV.q_out) annotation (Line(
            points={{80,0},{60,0}},
            color={0,0,0},
            thickness=1));
      end Ascites_PG;

      model Ascites_Resistance
        extends Tests.partialAscites(
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
          enable=false,
          useConductanceInput=false,
          Resistance(displayUnit="(mmHg.min)/l") = 7999343.2449*(15/1.5))
          annotation (Placement(transformation(extent={{-2,-50},{18,-30}})));
        Components.Ascites_Resistance ascites_Resistance
          annotation (Placement(transformation(extent={{-2,-86},{18,-66}})));
        SplanchnicCirculation splanchnicCirculation
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
            points={{90,0},{60,0}},
            color={0,0,0},
            thickness=1));
        connect(TIPSS.q_in, IntestineVenule.q_out) annotation (Line(
            points={{-2,-40},{-14,-40},{-14,2.22045e-16},{-20,2.22045e-16}},
            color={0,0,0},
            thickness=1));
        connect(TIPSS.q_out, HV.q_in) annotation (Line(
            points={{18,-40},{24,-40},{24,2.22045e-16},{40,2.22045e-16}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Resistance.q_in, unlimitedVolume1.y) annotation (Line(
            points={{-2,-76},{-86,-76},{-86,0},{-90,0}},
            color={0,0,0},
            thickness=1));
        connect(ascites_Resistance.q_out, HV.q_out) annotation (Line(
            points={{18,-76},{86,-76},{86,0},{60,0}},
            color={0,0,0},
            thickness=1));
        connect(splanchnicCirculation.q_out, HV.q_out) annotation (Line(
            points={{18,-124},{34,-124},{34,-76},{86,-76},{86,0},{60,0}},
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
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume unlimitedVolume1(P=
              13332.2387415)
          annotation (Placement(transformation(extent={{-104,-10},{-84,10}})));
        Physiolibrary.Hydraulic.Sources.UnlimitedVolume CVP(usePressureInput=
              true, P=666.611937075)
          annotation (Placement(transformation(extent={{116,-10},{96,10}})));
        Components.ResistancePressureDep resistancePressureDep(
          Comp(displayUnit="l/mmHg") = 6.0004926067653e-06,
          r0_nom(displayUnit="(mmHg.min)/l") = 719940892.041,
          P_nom=666.611937075)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Modelica.Blocks.Sources.RealExpression realExpression(y=max(time, 5)*
              133.322)
          annotation (Placement(transformation(extent={{92,30},{112,50}})));
      equation
        connect(unlimitedVolume1.y, resistancePressureDep.q_in) annotation (
            Line(
            points={{-84,0},{-10,0}},
            color={0,0,0},
            thickness=1));
        connect(CVP.y, resistancePressureDep.q_out) annotation (Line(
            points={{96,0},{10,0}},
            color={0,0,0},
            thickness=1));
        connect(realExpression.y, CVP.pressure) annotation (Line(points={{113,
                40},{150,40},{150,0},{116,0}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TestResistance;
    end Tests;

    package Components
      model Ascites_Resistance
        extends Tests.partialAscites(
          redeclare Physiolibrary.Hydraulic.Components.Resistor Liver(
              useConductanceInput=true, Resistance=6*mmHg/Qnom),
          redeclare Physiolibrary.Hydraulic.Components.Resistor IntestinesArt(
              Resistance=84*mmHg/Qnom),
          redeclare Physiolibrary.Hydraulic.Components.Resistor IntestineVenule(
              Resistance=3*mmHg/Qnom));
        Modelica.Blocks.Sources.RealExpression realExpression1(y=1/(time*mmHg/Qnom))
          annotation (Placement(transformation(extent={{-30,18},{-10,38}})));

        parameter Physiolibrary.Types.VolumeFlowRate Qnom=1.666666666666667e-08*(5000*
            0.2) "Nominal flow through the splanchnic circulation";
        parameter Physiolibrary.Types.Pressure mmHg=133.322;
        replaceable
        Physiolibrary.Hydraulic.Components.Resistor TIPSS(useConductanceInput=false,
            Resistance(displayUnit="(mmHg.min)/l") = 7999343.2449*(15/1.5))
          constrainedby Physiolibrary.Hydraulic.Interfaces.OnePort
          annotation (Placement(transformation(extent={{-2,-44},{18,-24}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a
                             q_in "Volume inflow" annotation (Placement(
              transformation(extent={{-114,-14},{-86,14}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_b
                             q_out "Volume outflow"
                               annotation (Placement(
              transformation(extent={{86,-14},{114,14}})));
      equation
        connect(Liver.cond, realExpression1.y) annotation (Line(points={{10,6},{10,14},
                {-2,14},{-2,28},{-9,28}},  color={0,0,127}));
        connect(TIPSS.q_in, IntestineVenule.q_out) annotation (Line(
            points={{-2,-34},{-14,-34},{-14,0},{-20,0}},
            color={0,0,0},
            thickness=1));
        connect(TIPSS.q_out, HV.q_in) annotation (Line(
            points={{18,-34},{24,-34},{24,0},{40,0}},
            color={0,0,0},
            thickness=1));
        connect(HV.q_out, q_out) annotation (Line(
            points={{60,0},{100,0}},
            color={0,0,0},
            thickness=1));
        connect(q_in, IntestinesArt.q_in) annotation (Line(
            points={{-100,0},{-80,0}},
            color={0,0,0},
            thickness=1));
        annotation (Documentation(info="<html>
<p>The TIPS resistance taken from TIPS flow and PVP from Su et al (2012, PMID 22099870).</p>
</html>"));
      end Ascites_Resistance;

      model ResistancePressureDep
        "Hydraulic resistor, dependent on the filling"
       extends Physiolibrary.Hydraulic.Interfaces.OnePort;
       extends Physiolibrary.Icons.HydraulicResistor;

        parameter Boolean enable=true   "if false, no resistance is used"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="External inputs/outputs"));

        parameter Boolean useConductanceInput = false
          "=true, if external conductance value is used"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="External inputs/outputs"));

        Physiolibrary.Types.HydraulicResistance r
          "Informative resistance value";
        constant Real ni = 4e-3 "Blood dynamic viscosity";
        constant Modelica.Units.SI.Length L = 1;
        parameter Physiolibrary.Types.HydraulicCompliance Comp;
        Modelica.Units.SI.Radius R;
        Physiolibrary.Types.Volume V = Modelica.Constants.pi * (R - R0)^2*L;

        Modelica.Units.SI.Radius R0 "Zero resistance, inferred from the P_nom at R0_nom";
        Modelica.Units.SI.Radius R_nom "Zero resistance, inferred from the P_nom at R0_nom";
        parameter Physiolibrary.Types.HydraulicResistance r0_nom "Nominal resistance at P_nom";
        parameter Physiolibrary.Types.Pressure P_nom "Nominal end-point pressure";
      equation
        // to calculate the R_0 parameter
        Modelica.Constants.pi * (R_nom - R0)^2*L = Comp*P_nom;
        r0_nom = 8*ni*L/(Modelica.Constants.pi*R_nom^4);

        r = 8*ni*L/(Modelica.Constants.pi*R^4);

        V = Comp*q_out.pressure;

        q_in.q = (q_in.pressure - q_out.pressure)/r;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
                  -100,-100},{100,100}}),
                         graphics={Text(
                extent={{-220,-40},{200,-80}},
                lineColor={0,0,255},
                fillColor={58,117,175},
                fillPattern=FillPattern.Solid,
                textString="%name")}),
          Documentation(revisions="<html>
<p><i>2009-2010</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>This hydraulic conductance (resistance) element contains two connector sides. No hydraulic medium volume is changing in this element during simulation. That means that sum of flow in both connector sides is zero. The flow through element is determined by <b>Ohm&apos;s law</b>. It is used conductance (=1/resistance) because it could be numerical zero better then infinity in resistance. </p>
</html>"));
      end ResistancePressureDep;
    end Components;

    model StarlingCurve
      Physiolibrary.Blocks.Factors.Spline
                                  rightStarling(data={{-6,0,0},{-3,0.15,0.104},{-1,0.52,
            0.48},{2,1.96,0.48},{4,2.42,0.123},{8,2.7,0}}, Xscale=101325/760)
        "At filling pressure 0mmHg (because external thorax pressure is -4mmHg) is normal cardiac output (effect=1)."
        annotation (Placement(transformation(extent={{-50,-6},{-30,14}})));
      Physiolibrary.Blocks.Factors.Spline
                                  leftStarling(data={{-4,0,0},{-1,0.72,0.29},{0,1.01,
            0.29},{3,1.88,0.218333},{10,2.7,0}}, Xscale=101325/760)
        "At filling pressure -0.0029mmHg (because external thorax pressure is -4mmHg) is normal cardiac output (effect=1)."
        annotation (Placement(transformation(extent={{-50,20},{-30,40}})));
      Physiolibrary.Types.Constants.VolumeFlowRateConst volumeFlowRate(k(
            displayUnit="l/min") = 8.3333333333333e-05)
        annotation (Placement(transformation(extent={{-66,60},{-58,68}})));
      Modelica.Blocks.Sources.RealExpression realExpression(y=time*mmHg)
        annotation (Placement(transformation(extent={{-100,20},{-80,40}})));
      constant Real mmHg=133.322;

    Physiolibrary.Types.Pressure preload = realExpression.y + p_th;
    parameter Physiolibrary.Types.Pressure p_th=533.28954966
                                                    "Thoracic cavity pressure";
    Physiolibrary.Types.VolumeFlowRate CO_LV = leftStarling.y;
    Physiolibrary.Types.VolumeFlowRate CO_RV = rightStarling.y;
    equation
      connect(realExpression.y, leftStarling.u)
        annotation (Line(points={{-79,30},{-48,30}}, color={0,0,127}));
      connect(realExpression.y, rightStarling.u) annotation (Line(points={{-79,30},{
              -70,30},{-70,4},{-48,4}}, color={0,0,127}));
      connect(volumeFlowRate.y, leftStarling.yBase)
        annotation (Line(points={{-57,64},{-40,64},{-40,32}}, color={0,0,127}));
      connect(volumeFlowRate.y, rightStarling.yBase)
        annotation (Line(points={{-57,64},{-40,64},{-40,6}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StartTime=-5,
          StopTime=25,
          __Dymola_Algorithm="Dassl"));
    end StarlingCurve;

    model Kofranek2014
      extends
        Physiolibrary.Hydraulic.Examples.Kofranek2014.NonPulsatileCirculation(
          RT(k(displayUnit="(mmHg.min)/l") = Modelica.Constants.inf));
      Physiolibrary.Hydraulic.Components.Conductor
               nonMuscle(Conductance(displayUnit="m3/(Pa.s)") = 2.6e-9)
        annotation (Placement(transformation(extent={{-4,-112},{16,-92}})));
      Physiolibrary.Hydraulic.Components.Conductor
               muscle(Conductance(displayUnit="l/(mmHg.min)") =
          1.3001067314658e-09)
        annotation (Placement(transformation(extent={{-4,-94},{16,-74}})));
      Components.Ascites_Resistance ascites_Resistance(
        TIPSS(enable=true, Resistance=Modelica.Constants.inf),
        Liver(enable=true, useConductanceInput=true),
        realExpression1(y=1/(max(time/200, 5)*ascites_Resistance.mmHg/
              ascites_Resistance.Qnom)))
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
  end Hemodynamics;

  package AscitesLevitt "Ascites model by Levitt and Levitt"
    model LevittCase1SS
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
    end LevittCase1SS;

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

    model LD_concs
      extends AscitesLevitt.LevitCase2Dynamics(
        Ap=22.0,
        Pgrad2=5.0,
        Pvg=0.0);
      constant Physiolibrary.Types.Pressure mmHg=133.322387415;
      OncoticPressures.POnc2Conc pOnc2Conc(
        inputIsPressure=false,
        inputValue=4,
        AGf=4/3)
        annotation (Placement(transformation(extent={{-40,-20},{-20,0}})));
      OncoticPressures.POnc2Conc protein_asictes(
        inputIsPressure=true,
        inputValue=Aa,
        AGf=4/3) annotation (Placement(transformation(extent={{0,-20},{20,0}})));
      OncoticPressures.POnc2Conc protein_intes(
        inputIsPressure=true,
        inputValue=Ai,
        AGf=4/3)
        annotation (Placement(transformation(extent={{40,-20},{60,0}})));
    end LD_concs;

    model LSS_concs
      extends AscitesLevitt.LevittCase1SS;
      OncoticPressures.POnc2Conc protein_plasma(inputIsPressure=true,
          inputValue=Ap)
        annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
      OncoticPressures.POnc2Conc protein_asc(inputIsPressure=true, inputValue=
            Aa) annotation (Placement(transformation(extent={{-20,20},{0,40}})));
    end LSS_concs;

    package Test
      model Healthy
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

    model LevittCase1SS_SI
      "Model of dynamic ascites build-up by Levitt and Levitt (2012), PMID 22453061. Converted from Maple code."
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
    end LevittCase1SS_SI;

    model LevittCase1SsSiIo
      "Model of static ascites build-up by Levitt and Levitt (2012), PMID 22453061. Converted from Maple code into SI units and IO"
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
      Physiolibrary.Types.RealIO.PressureInput Pp "portal vein pressure input"
         annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
      Physiolibrary.Types.RealIO.PressureInput Pra "Right atrial pressure input"
        annotation (Placement(transformation(extent={{80,-20},{120,20}})));
      Physiolibrary.Types.RealIO.PressureOutput Pa "Ascites pressure output"
        annotation (Placement(transformation(extent={{-20,-100},{20,-60}})));
      Physiolibrary.Types.RealIO.PressureInput Phv "Hepatic vein pressure input"
        annotation (Placement(transformation(extent={{-120,-100},{-80,-60}})));
      Physiolibrary.Types.RealIO.PressureInput Pc "Intestinal capillary pressure"
        annotation (Placement(transformation(extent={{-120,60},{-80,100}})));
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
          __Dymola_Algorithm="Dassl"));
    end LevittCase1SsSiIo;

    model LevittCase1SS_SIInverted
      "Model of dynamic ascites build-up by Levitt and Levitt (2012), PMID 22453061. Converted from Maple code. Inverted the causality of HVPG as an output."
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

    Physiolibrary.Types.Pressure Pa=time*133.322
                                    annotation (Dialog(tab="General", group="Inputs"));
    Physiolibrary.Types.Pressure Aa;
    Physiolibrary.Types.Pressure Pgrad annotation (Dialog(tab="General", group="Inputs"));

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
          StopTime=25,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Dassl"));
    end LevittCase1SS_SIInverted;
  end AscitesLevitt;

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
