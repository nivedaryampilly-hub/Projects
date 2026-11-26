%% =====================================================================
%  REGENERATIVE BRAKING WITH SUPERCAPACITOR  -  SIMULINK PROJECT
%  E-SCOOTER (S1 Pro class)  -  REAL KATPADI DRIVE CYCLE
%  ---------------------------------------------------------------------
%  Run this ONE script. It will
%    1. read the real Katpadi drive cycle from the Excel file
%    2. build the Simulink model  RegenBraking_SC.slx  programmatically
%       (motor/inverter model, 14S16P battery model with SOC, EMS,
%        PID-controlled bidirectional DC-DC converter, supercapacitor)
%    3. auto-tune the PID current controller (fminsearch, ITAE cost)
%    4. sweep the supercapacitor size, re-tuning the PID for each size
%    5. pick the optimum size (min battery Wh/km + weighted SC cost)
%    6. run the optimum design, print the results and draw all graphs
%       (including the battery SOC and SC state-of-energy plot)
%
%  Simulink model structure (open it after the run to inspect / edit):
%
%   Drive cycle (v,a) --> [EMS: motor model + rules] --> i_ref --> (+/-) --> [PID Controller block]
%                              |  P_bus, P_fric                      ^               |
%                              v                                     | iL            v
%                        [Battery: OCV+R0+RC, SOC] <-- D, iL -- [Duty = Vc/Vbus + u, limited]
%                              |  Vbus (filtered)                    |
%             Bidirectional DC-DC (averaged):   L*diL/dt = D*Vbus - Vc - (ESR+RL)*iL
%             Supercapacitor:                    C*dVc/dt = iL - Vc/Rleak
%
%  Put Indian_Drive_Cycle_Full_1180s_1s.xlsx in the same folder.
%  Requires: MATLAB + Simulink (R2019b or newer recommended). No other toolboxes.
%  Runtime : long (one stiff-solver simulation of 1180 s per SC size).
%            Shorten Cgrid below if you want a quicker run.
%  =====================================================================
clear; clc; close all;
mdl = 'RegenBraking_SC';

%% 1. PARAMETERS ---------------------------------------------------------
% ---- vehicle (scooter + rider)
P.m0 = 200;  P.k_aero = 0.5*1.2*0.54;  P.Crr = 0.012;  P.g = 9.81;
P.rw = 0.233;  P.G = 4.0;  P.eta_g = 0.97;
% ---- motor & inverter (MCU) from datasheet screenshot
P.T_pk = 58;  P.P_pk = 11e3;  P.I_pk = 65;  P.Kt = P.T_pk/P.I_pk;
P.ke_LL = 12.4e-3*60/(2*pi);                      % 12.4 V/krpm in V.s/rad
P.Rs = 0.045;  P.kfe = 0.05;  P.kw = 2e-4;  P.R_inv = 0.004;  P.k_sw = 0.06;  P.Ifw = 25;
% ---- battery pack from datasheet screenshot (14S16P NMC811, 21700)
P.Ns = 14;  P.Np = 16;  P.Q = 76.8;  P.Vbus = 51.8;
P.R0 = 0.026;  P.R1 = 0.010;  P.tau1 = 30;  P.C1 = P.tau1/P.R1;
P.soc0 = 0.80;  P.tau_bus = 5e-3;
P.soc_tab = [0 .05 .10 .20 .30 .40 .50 .60 .70 .80 .90 1];
P.ocv_tab = [3.00 3.35 3.50 3.60 3.66 3.72 3.78 3.85 3.93 4.03 4.12 4.20];
% ---- energy management
P.Pb_chg = 2e3;  P.Pb_share = 0.6;  P.P_th = 3e3;  P.Paux = 100;  P.v_lo = 1;  P.v_hi = 3;
% ---- supercapacitor bank
P.Vmax = 40;  P.Vmin = 20;  P.Vc0 = 30;  P.Imax = 250;
P.ESRC = 1.134;  P.RleakC = 2e4;  P.Whkg = 2.3;  P.cost_kJ = 3;  P.w_cost = 0.004;
% ---- converter + controller
P.L = 1e-4;  P.RL = 0.003;  P.dt_c = 5e-4;  P.dt_e = 1e-2;  P.Dmin = 0.02;  P.Dmax = 0.98;

%% 2. DRIVE CYCLE (real Katpadi drive cycle from Excel) -------------------
cycFile = 'Indian_Drive_Cycle_Full_1180s_1s.xlsx';
if ~isfile(cycFile)
    error('Drive-cycle file "%s" not found. Put it in the current folder.', cycFile);
end
Dc = readmatrix(cycFile);                         % Time(s), Velocity(km/hr)
Dc = Dc(all(isfinite(Dc),2), 1:2);
cyc.t = (0:P.dt_e:Dc(end,1))';
cyc.v = interp1(Dc(:,1), Dc(:,2)/3.6, cyc.t, 'linear');
cyc.a = gradient(cyc.v, P.dt_e);
T = cyc.t(end);
setBase('ts_v', timeseries(cyc.v, cyc.t));
setBase('ts_a', timeseries(cyc.a, cyc.t));
dist_km = trapz(cyc.t, cyc.v)/1000;
fprintf('Katpadi drive cycle: %.0f s, %.2f km, mean %.1f km/h, max %.0f km/h\n', ...
        T, dist_km, 3.6*dist_km*1000/T, 3.6*max(cyc.v));

%% 3. BUILD THE SIMULINK MODEL ------------------------------------------
setBase('pvec', makePvec(50, P));  setBase('Vc0', P.Vc0);
setBase('soc0', P.soc0);  setBase('Vbf0', ocvPack(P.soc0, P));
g0 = tunePID(50, P);  setBase('Kp', g0(1)); setBase('Ki', g0(2)); setBase('Kd', g0(3));
buildModel(mdl, T);
fprintf('Model "%s" built.\n', mdl);

%% 4. BASELINES ----------------------------------------------------------
B_fr = baselineRun(P, cyc, false);                % friction-only braking
B_bt = baselineRun(P, cyc, true);                 % battery-only regen
Wh_fric = B_fr.WhKm;  Wh_batt = B_bt.WhKm;
Pw0 = wheelPower(cyc.v, cyc.a, P.m0, P);
E_brake_Wh = trapz(cyc.t, max(-Pw0,0))/3600;

%% 5. SC SIZE SWEEP (PID re-tuned for every size) -------------------------
Cgrid = [10 20 35 50 75 100 150];
nC = numel(Cgrid);  WhKm = zeros(nC,1);  Jobj = zeros(nC,1);  G = zeros(nC,3);  EscChg = zeros(nC,1);
fprintf('\nSweeping supercapacitor size with Simulink...\n');
for i = 1:nC
    [M, ~] = runCase(mdl, Cgrid(i), P, cyc);
    G(i,:) = M.gains;  WhKm(i) = M.WhPerKm;  EscChg(i) = M.EscChg_Wh;
    Jobj(i) = M.WhPerKm + P.w_cost*scCost(Cgrid(i), P);
    fprintf('  C = %5.1f F | %6.2f Wh/km | J = %6.2f | Kp=%.5f Ki=%.3f\n', ...
            Cgrid(i), WhKm(i), Jobj(i), G(i,1), G(i,2));
end

%% 6. OPTIMUM (interpolated minimum, then verified by a final simulation) -
Cfine = linspace(min(Cgrid), max(Cgrid), 500);
Jfine = interp1(Cgrid, Jobj, Cfine, 'pchip');
[J_opt, ik] = min(Jfine);  C_opt = Cfine(ik);
fprintf('\nOptimum capacitance (interpolated): C* = %.2f F -> verifying...\n', C_opt);
[Mopt, S] = runCase(mdl, C_opt, P, cyc);
J_opt = Mopt.WhPerKm + P.w_cost*scCost(C_opt, P);

%% 7. PRINT RESULTS -------------------------------------------------------
save_fr = Wh_fric*dist_km - Mopt.Enet_Wh;
save_bt = Wh_batt*dist_km - Mopt.Enet_Wh;
fprintf('\n================ FINAL RESULTS (optimal design) ================\n');
fprintf(' Optimal SC capacitance          : %8.2f F\n',    C_opt);
fprintf(' SC ESR                          : %8.1f mOhm\n', 1e3*P.ESRC/C_opt);
fprintf(' SC bank mass                    : %8.1f kg\n',   scMass(C_opt,P));
fprintf(' SC bank rated energy            : %8.1f kJ\n',   0.5*C_opt*P.Vmax^2/1e3);
fprintf(' SC bank cost (estimate)         : %8.0f $\n',    scCost(C_opt,P));
fprintf(' PID gains  Kp / Ki / Kd         : %.5f / %.3f / %.2e\n', Mopt.gains);
fprintf(' Current-loop RMS tracking error : %8.2f A\n',    Mopt.rmsErr);
fprintf('----------------------------------------------------------------\n');
fprintf(' Braking energy at wheels        : %8.1f Wh\n',  E_brake_Wh);
fprintf(' Friction-brake energy wasted    : %8.1f Wh\n',  Mopt.Efric_Wh);
fprintf(' Energy stored into SC           : %8.1f Wh\n',  Mopt.EscChg_Wh);
fprintf(' SC + converter losses           : %8.1f Wh\n',  Mopt.Eloss_Wh);
fprintf('----------------------------------------------------------------\n');
fprintf(' Consumption, friction-only      : %8.2f Wh/km\n', Wh_fric);
fprintf(' Consumption, battery-only regen : %8.2f Wh/km\n', Wh_batt);
fprintf(' Consumption, battery + SC       : %8.2f Wh/km\n', Mopt.WhPerKm);
fprintf(' Saving vs friction-only         : %8.1f %%\n', 100*save_fr/(Wh_fric*dist_km));
fprintf(' Saving vs battery-only regen    : %8.1f %%\n', 100*save_bt/(Wh_batt*dist_km));
fprintf(' Braking-energy recovery ratio   : %8.1f %%\n', 100*save_fr/E_brake_Wh);
fprintf(' Peak battery power (batt only)  : %8.1f kW\n', max(abs(B_bt.Pbatt))/1e3);
fprintf(' Peak battery power (batt + SC)  : %8.1f kW\n', max(abs(S.Pbatt))/1e3);
fprintf(' Battery RMS current (batt only) : %8.1f A\n',  sqrt(mean(B_bt.Ib.^2)));
fprintf(' Battery RMS current (batt + SC) : %8.1f A\n',  sqrt(mean(S.Ib.^2)));
fprintf('----------------------------------------------------------------\n');
fprintf(' Battery SOC start               : %8.2f %%\n', 100*P.soc0);
fprintf(' Battery SOC end, friction-only  : %8.2f %%\n', 100*B_fr.soc(end));
fprintf(' Battery SOC end, battery regen  : %8.2f %%\n', 100*B_bt.soc(end));
fprintf(' Battery SOC end, battery + SC   : %8.2f %%\n', 100*S.soc(end));
fprintf('================================================================\n');

%% 8. GRAPHS ---------------------------------------------------------------
t = S.t;
c1 = [0.00 0.45 0.74]; c2 = [0.85 0.33 0.10]; c3 = [0.47 0.67 0.19]; c4 = [0.49 0.18 0.56];

figure('Name','Speed profile','Color','w'); plot(t, S.v*3.6, 'Color', c1, 'LineWidth', 1.4); grid on;
xlabel('Time [s]'); ylabel('Speed [km/h]'); title('Real Katpadi drive cycle');
figure('Name','Wheel power','Color','w'); Pw = wheelPower(S.v, gradient(S.v, t), P.m0 + scMass(C_opt,P), P)/1e3;
area(t, max(Pw,0), 'FaceColor', c2, 'EdgeColor','none'); hold on;
area(t, min(Pw,0), 'FaceColor', c3, 'EdgeColor','none'); grid on;
title('Wheel power (traction / braking)'); ylabel('Wheel power [kW]'); xlabel('Time [s]'); legend('Traction','Braking');

figure('Name','SC voltage','Color','w'); plot(t, S.Vc, 'Color', c1, 'LineWidth', 1.4); hold on;
yline(P.Vmax,'r--','V_{max}'); yline(P.Vmin,'k--','V_{min}'); grid on;
xlabel('Time [s]'); ylabel('V_{SC} [V]'); title(sprintf('Supercapacitor voltage (C = %.1f F)', C_opt));
figure('Name','SC current tracking (PID)','Color','w'); plot(t, S.iref, 'Color', c2, 'LineWidth', 1.6); hold on;
plot(t, S.iL, 'Color', c1, 'LineWidth', 1.0); grid on;
xlabel('Time [s]'); ylabel('SC current [A]'); legend('Reference','Measured (PID output)');
title('PID current tracking (+ = charging, regen)');
figure('Name','Converter duty cycle','Color','w'); plot(t, S.D, 'Color', c4); grid on;
ylabel('Duty cycle'); xlabel('Time [s]'); title('Converter duty cycle');

figure('Name','Battery power','Color','w'); plot(cyc.t, B_bt.Pbatt/1e3, 'Color', [.6 .6 .6], 'LineWidth', 1.2); hold on;
plot(t, S.Pbatt/1e3, 'Color', c1, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Battery power [kW]'); legend('Battery (no SC)','Battery (with SC)');
title('Battery power: SC shaves peaks and absorbs regen');
figure('Name','Supercapacitor power','Color','w'); plot(t, S.Psc/1e3, 'Color', c3, 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('[kW]'); title('Supercapacitor power (+ charging)');
figure('Name','Friction brake power','Color','w'); area(t, S.Pfric/1e3, 'FaceColor', c2, 'EdgeColor','none'); grid on;
ylabel('[kW]'); xlabel('Time [s]'); title('Power dissipated in friction brakes');

figure('Name','Energy vs SC size','Color','w'); plot(Cgrid, WhKm, 'o-', 'Color', c1, 'LineWidth', 1.5, 'MarkerFaceColor', c1); hold on;
yline(Wh_batt,'--','Battery-only regen','Color',[.4 .4 .4]);
yline(Wh_fric,'--','Friction only','Color',c2); grid on;
xlabel('C [F]'); ylabel('Wh/km'); title('Energy consumption vs. SC size');
figure('Name','Objective vs SC size','Color','w'); plot(Cfine, Jfine, 'Color', c4, 'LineWidth', 1.5); hold on;
plot(Cgrid, Jobj, 'o', 'Color', c4);
plot(C_opt, J_opt, 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k'); grid on;
xlabel('C [F]'); ylabel('Objective J'); title(sprintf('Objective (optimum C* = %.1f F)', C_opt));
figure('Name','Captured energy vs SC mass','Color','w'); yyaxis left; plot(Cgrid, EscChg, 's-', 'LineWidth', 1.4); ylabel('Energy into SC [Wh]');
yyaxis right; plot(Cgrid, arrayfun(@(c) scMass(c,P), Cgrid), '^-', 'LineWidth', 1.4);
ylabel('SC mass [kg]'); grid on; xlabel('C [F]'); title('Captured energy vs. added mass');
figure('Name','PID gains vs SC size','Color','w'); yyaxis left; semilogy(Cgrid, G(:,1), 'o-', 'LineWidth', 1.4); ylabel('K_p');
yyaxis right; semilogy(Cgrid, G(:,2), 's-', 'LineWidth', 1.4); ylabel('K_i');
grid on; xlabel('C [F]'); title('Auto-tuned PID gains vs. SC size');

Rt = P.ESRC/C_opt + P.RL;  Kp0 = 2*pi*150*P.L/P.Vbus;
[ts0, i0] = stepSim([Kp0, Kp0*Rt/P.L, Kp0*1e-4], Rt, P);  [ts1, i1] = stepSim(Mopt.gains, Rt, P);
figure('Name','PID step response','Color','w','Position',[170 90 900 500]);
plot(ts0*1e3, i0, '--', 'Color', c2, 'LineWidth', 1.4); hold on;
plot(ts1*1e3, i1, 'Color', c1, 'LineWidth', 1.8); yline(100,'k:'); grid on;
xlabel('Time [ms]'); ylabel('Inductor current [A]');
legend('Initial analytic gains','Optimised gains','Reference','Location','southeast');
title('Current-loop step response (100 A step)');

figure('Name','Energy consumption comparison','Color','w'); b = bar([Wh_fric Wh_batt Mopt.WhPerKm], 0.6); b.FaceColor = 'flat';
b.CData = [c2; [.6 .6 .6]; c3]; grid on;
set(gca,'XTickLabel',{'Friction only','Battery regen','Battery + SC'}); ylabel('Wh/km');
title('Energy consumption comparison');
text(1:3, [Wh_fric Wh_batt Mopt.WhPerKm], compose('%.1f',[Wh_fric Wh_batt Mopt.WhPerKm]), ...
     'HorizontalAlignment','center','VerticalAlignment','bottom');
figure('Name','Braking energy split','Color','w'); pie([Mopt.Efric_Wh, Mopt.EscChg_Wh, max(E_brake_Wh - Mopt.Efric_Wh - Mopt.EscChg_Wh,0)+eps]);
legend({'Friction brakes (lost)','Into supercapacitor','Battery regen + losses'}, 'Location','southoutside');
title('Where the braking energy goes');

% ---- SOC graph (battery SOC, effective SOC, SOC rate, SC state of energy)
SoE    = 100*max(min((S.Vc.^2 - P.Vmin^2)/(P.Vmax^2 - P.Vmin^2), 1), 0);
% effective SOC = battery SOC + SC stored energy (relative to start) expressed in battery SOC
socEff = S.soc + (0.5*C_opt*(S.Vc.^2 - P.Vc0^2)/3600)/(P.Q*P.Vbus);
win    = t <= 300;                                     % zoom window with several stop-go events
rate   = @(s) 6000*movmean(gradient(s, t(win)), 300);  % dSOC/dt in %/min (s = windowed SOC), 3 s moving average
figure('Name','SOC full cycle','Color','w'); plot(t, 100*B_fr.soc, 'Color', c2, 'LineWidth', 1.4); hold on;
plot(t, 100*B_bt.soc, 'Color', [.5 .5 .5], 'LineWidth', 1.4);
plot(t, 100*S.soc, 'Color', c3, 'LineWidth', 1.8);
plot(t, 100*socEff, '--', 'Color', c1, 'LineWidth', 1.6); grid on;
xlabel('Time [s]'); ylabel('SOC [%]'); legend('Friction only','Battery regen','Battery + SC (battery only)', ...
       'Battery + SC (effective, incl. SC energy)','Location','northeast');
title('Battery state of charge over the Katpadi cycle');
figure('Name','SOC zoom 0-300 s','Color','w'); ax2 = gca; plot(t(win), 100*B_fr.soc(win), 'Color', c2, 'LineWidth', 1.4); hold on;
plot(t(win), 100*B_bt.soc(win), 'Color', [.5 .5 .5], 'LineWidth', 1.6);
plot(t(win), 100*S.soc(win), 'Color', c3, 'LineWidth', 1.8);
plot(t(win), 100*socEff(win), '--', 'Color', c1, 'LineWidth', 1.8); grid on;
yl = [min([B_fr.soc(win); B_bt.soc(win); S.soc(win); socEff(win)]), ...
      max([B_fr.soc(win); B_bt.soc(win); S.soc(win); socEff(win)])]*100;
ylim(yl + 0.05*diff(yl)*[-1 1]);  shadeBraking(t(win), Pw0(win), ax2);
xlabel('Time [s]'); ylabel('SOC [%]'); title('Zoom (0-300 s): SOC falls while accelerating, rises while braking (green bands)');
figure('Name','SOC rate','Color','w'); ax3 = gca; h1 = plot(t(win), rate(B_bt.soc(win)), 'Color', [.5 .5 .5], 'LineWidth', 1.6); hold on;
h2 = plot(t(win), rate(S.soc(win)), 'Color', c3, 'LineWidth', 1.6);
h3 = plot(t(win), rate(socEff(win)), '--', 'Color', c1, 'LineWidth', 1.6);
yline(0, 'k-'); grid on;
shadeBraking(t(win), Pw0(win), ax3);
xlabel('Time [s]'); ylabel('dSOC/dt [%/min]'); title('SOC rate: negative = discharging (traction), positive = charging (regen)');
legend([h1 h2 h3], {'Battery regen','Battery + SC (battery only)','Battery + SC (effective)'}, 'Location','best');
figure('Name','SOC gain from SC','Color','w'); plot(t, 100*(S.soc - B_bt.soc), 'Color', c1, 'LineWidth', 1.4); grid on;
xlabel('Time [s]'); ylabel('\DeltaSOC [%]'); title('SOC gained by adding the supercapacitor (vs. battery-only regen)');
figure('Name','SC state of energy','Color','w'); plot(t, SoE, 'Color', c4, 'LineWidth', 1.4); grid on;
ylabel('SC state of energy [%]'); xlabel('Time [s]');
title(sprintf('Supercapacitor state of energy (C = %.1f F)', C_opt));

open_system(mdl);  save_system(mdl);
fprintf('\nSimulink model saved as %s.slx\n', mdl);

%% =====================================================================
%                          LOCAL FUNCTIONS
%  =====================================================================
function setBase(name, val), assignin('base', name, val); end

function [M, S] = runCase(mdl, C, P, cyc)
% Tune PID for this C, push parameters to base workspace, run Simulink.
g = tunePID(C, P);
setBase('pvec', makePvec(C, P));  setBase('Vc0', P.Vc0);
setBase('soc0', P.soc0);  setBase('Vbf0', ocvPack(P.soc0, P));
setBase('Kp', g(1));  setBase('Ki', g(2));  setBase('Kd', g(3));
so = sim(mdl);
T  = cyc.t(end);  tg = (0:0.01:T)';
names = {'v','Vc','iL','iref','D','Pbus','Pfric','Pbatt','Psc','soc','Ib','Vbus'};
S.t = tg;
for k = 1:numel(names), S.(names{k}) = getSig(so, ['out_' names{k}], tg); end
Eb    = getSig(so, 'out_Eb',    tg);  Efric = getSig(so, 'out_Efric', tg);
Eloss = getSig(so, 'out_Eloss', tg);
dEsc  = 0.5*C*(P.Vc0^2 - S.Vc(end)^2)/3600;
M.gains     = g;
M.Enet_Wh   = Eb(end) + dEsc;
M.WhPerKm   = M.Enet_Wh/(trapz(tg, S.v)/1000);
M.Efric_Wh  = Efric(end);
M.Eloss_Wh  = Eloss(end);
M.EscChg_Wh = trapz(tg, max(S.Psc,0))/3600;
M.rmsErr    = sqrt(mean((S.iref - S.iL).^2));
end

function y = getSig(so, name, tg)
ts = so.get(name);  t = ts.Time(:);  y = squeeze(ts.Data);  y = y(:);
[t, iu] = unique(t);  y = y(iu);
y = interp1(t, y, tg, 'linear', 'extrap');
end

function p = makePvec(C, P)
% Parameter vector read by every MATLAB Function block in the model.
p = [P.m0 + scMass(C,P), P.k_aero, P.Crr, P.g, P.rw, P.G, P.eta_g, ...
     P.T_pk, P.P_pk, P.Kt, P.ke_LL, P.Rs, P.kfe, P.kw, P.R_inv, P.k_sw, P.Ifw, ...
     P.Vmax, P.Vmin, P.Imax, P.P_th, P.Pb_chg, P.v_lo, P.v_hi, P.Paux, ...
     P.ESRC/C, P.RL, C, P.L, P.RleakC/C, ...
     P.R0, P.R1, P.C1, P.Q, P.tau_bus, P.Ns, P.Pb_share];
%    1 m  2 kaero 3 Crr 4 g 5 rw 6 G 7 eta_g 8 T_pk 9 P_pk 10 Kt 11 ke 12 Rs 13 kfe
%    14 kw 15 R_inv 16 k_sw 17 Ifw 18 Vmax 19 Vmin 20 Imax 21 P_th 22 Pb_chg 23 v_lo
%    24 v_hi 25 Paux 26 ESR 27 RL 28 C 29 L 30 Rleak 31 R0 32 R1 33 C1 34 Q
%    35 tau_bus 36 Ns 37 Pb_share
end

% ----------------------------------------------------------------------
function buildModel(mdl, T)
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl);
set_param(mdl, 'Solver','ode23tb', 'StopTime',num2str(T), 'MaxStep','0.01', ...
    'RelTol','1e-4', 'AbsTol','1e-6', 'ReturnWorkspaceOutputs','on', 'SaveOutput','off');
L = @(s,d) add_line(mdl, s, d, 'autorouting', 'on');

% ---- sources
addB(mdl,'simulink/Sources/From Workspace','v_ref',[30 40 120 70], ...
     'VariableName','ts_v','OutputAfterFinalValue','Holding final value');
addB(mdl,'simulink/Sources/From Workspace','a_ref',[30 110 120 140], ...
     'VariableName','ts_a','OutputAfterFinalValue','Holding final value');
addB(mdl,'simulink/Sources/Constant','pvec',[30 560 120 590],'Value','pvec');

% ---- Energy management system: motor model + rule-based supervisory controller
emsCode = strjoin({
'function [P_bus, i_ref, P_fric] = ems(v, a, Vsc, Vb, socB, p)'
'% Energy management: + i_ref charges the supercapacitor (regen).'
'Vmax = p(18); Vmin = p(19); Imax = p(20); Pth0 = p(21); Pbchg = p(22);'
'Paux = p(25); eta = p(7); Pbsh = p(37);'
'P_bus = 0; i_ref = 0; P_fric = 0;'
'[Pw, Psh, Pe] = motorpoint(v, a, Vb, p);'
'Vs = max(Vsc, 1);'
'socS = (Vs - Vmin)/(Vmax - Vmin);'
'if Pw >= 0'
'    share = min(max((socS - 0.25)/0.40, 0), 1);'
'    Pth = Pth0*(0.5 + 0.5*min(max((socB - 0.1)/0.2, 0), 1));'
'    Psc = share*Pe + (1 - share)*max(Pe - Pth, 0);'
'    kd = min(max((Vs - Vmin)/(0.05*Vmax), 0), 1);'
'    i_ref = -kd*min(Psc/Vs, Imax);'
'    P_bus = Pe + Paux;'
'else'
'    Pbc = Pbchg*min(max((0.98 - socB)/0.05, 0), 1);'
'    kc = min(max((Vmax - Vs)/(0.04*Vmax), 0), 1);'
'    Pbt = min(Pbsh*Pe, Pbc);'
'    i_ref = min((Pe - Pbt)/Vs, Imax*kc);'
'    Pbat = min(Pe - Vs*i_ref, Pbc);'
'    Pfin = Vs*i_ref + Pbat;'
'    s = 0;'
'    if Pe > 0'
'        s = Pfin/Pe;'
'    end'
'    P_fric = -Pw - s*Psh/eta;'
'    P_bus = -Pfin + Paux;'
'end'
'end'
''
'function [Pw, Psh, Pe] = motorpoint(v, a, Vb, p)'
'm = p(1); kaero = p(2); Crr = p(3); g = p(4); rw = p(5); G = p(6); eta = p(7);'
'Tpk = p(8); Ppk = p(9); vlo = p(23); vhi = p(24);'
'F = m*a + kaero*v^2 + Crr*m*g*double(v > 0.05);'
'Pw = F*v;'
'wm = max(G*v/rw, 1e-3);'
'Psh = 0; Pe = 0;'
'if Pw >= 0'
'    Tm = F*rw/(G*eta);'
'    Psh = Tm*wm;'
'    if Psh > Ppk'
'        Psh = Ppk; Tm = Psh/wm;'
'    end'
'    Tm = min(Tm, Tpk);'
'    Pe = Psh + motorloss(Tm, wm, Vb, p);'
'else'
'    blend = min(max((v - vlo)/(vhi - vlo), 0), 1);'
'    Tm = -F*rw/G*eta*blend;'
'    Psh = Tm*wm;'
'    if Psh > Ppk'
'        Psh = Ppk; Tm = Psh/wm;'
'    end'
'    if Tm > Tpk'
'        Tm = Tpk; Psh = Tm*wm;'
'    end'
'    Pe = max(Psh - motorloss(Tm, wm, Vb, p), 0);'
'end'
'end'
''
'function Ploss = motorloss(T, wm, Vb, p)'
'Iq = abs(T)/p(10);'
'wnl = Vb/p(11);'
'Id = p(17)*min(max((wm/wnl - 0.9)/0.5, 0), 1);'
'I2 = Iq^2 + Id^2;'
'Ploss = 1.5*p(12)*I2 + p(13)*wm^1.5 + p(14)*wm^2 + 1.5*p(15)*I2 + p(16)*sqrt(I2);'
'end'}, newline);
addMLF(mdl,'EMS',[220 40 380 280], emsCode);

% ---- PID current controller
addB(mdl,'simulink/Math Operations/Sum','ErrSum',[450 125 480 155],'Inputs','+-');
pid = addB(mdl,'simulink/Continuous/PID Controller','PID_Current',[520 110 610 170], ...
     'Controller','PID','P','Kp','I','Ki','D','Kd','N','1000', ...
     'LimitOutput','on','UpperSaturationLimit','0.5','LowerSaturationLimit','-0.5');
try, set_param(pid,'UseFilter','on'); catch, end
try, set_param(pid,'AntiWindupMode','clamping'); catch, end

addMLF(mdl,'Duty',[660 105 780 185], strjoin({
'function D = duty(u, Vc, Vb)'
'% Voltage feed-forward (actual DC-link voltage) + PID correction, limited to [0.02, 0.98]'
'D = min(max(Vc/Vb + u, 0.02), 0.98);'
'end'}, newline));

% ---- plant: bidirectional DC-DC (averaged) + supercapacitor
addMLF(mdl,'PlantL',[830 105 950 245], strjoin({
'function diL = plantL(D, Vc, iL, Vb, p)'
'diL = (D*Vb - Vc - (p(26) + p(27))*iL)/p(29);'
'end'}, newline));
addB(mdl,'simulink/Continuous/Integrator','iL_int',[1000 150 1030 180],'InitialCondition','0');
addMLF(mdl,'PlantC',[830 290 950 370], strjoin({
'function dVc = plantC(iL, Vc, p)'
'dVc = (iL - Vc/p(30))/p(28);'
'end'}, newline));
addB(mdl,'simulink/Continuous/Integrator','Vc_int',[1000 315 1030 345],'InitialCondition','Vc0');

% ---- battery: OCV(SOC) + R0 + RC, SOC by coulomb counting, DC-link voltage filter
addMLF(mdl,'Battery',[660 400 840 600], strjoin({
'function [dsoc, dv1, dVbf, Pbatt, Pchem, Ib, Vt] = battery(P_bus, D, iL, Vbf, soc, v1, p)'
'Pt = P_bus + Vbf*D*iL;'
'Voc = ocvpack(soc, p(36));'
'E = Voc - v1;'
'disc = max(E^2 - 4*p(31)*Pt, 0);'
'Ib = (E - sqrt(disc))/(2*p(31));'
'Vt = E - p(31)*Ib;'
'dsoc = -Ib/(3600*p(34));'
'dv1 = Ib/p(33) - v1/(p(33)*p(32));'
'dVbf = (Vt - Vbf)/p(35);'
'Pbatt = Pt;'
'Pchem = Voc*Ib;'
'end'
''
'function V = ocvpack(soc, Ns)'
'xs = [0 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00];'
'ys = [3.00 3.35 3.50 3.60 3.66 3.72 3.78 3.85 3.93 4.03 4.12 4.20];'
's = min(max(soc, 0), 1);'
'V = Ns*ys(12);'
'for k = 1:11'
'    if s <= xs(k+1)'
'        V = Ns*(ys(k) + (ys(k+1) - ys(k))*(s - xs(k))/(xs(k+1) - xs(k)));'
'        break;'
'    end'
'end'
'end'}, newline));
addB(mdl,'simulink/Continuous/Integrator','SOC_int',[900 400 930 430],'InitialCondition','soc0');
addB(mdl,'simulink/Continuous/Integrator','V1_int',[900 450 930 480],'InitialCondition','0');
addB(mdl,'simulink/Continuous/Integrator','Vbf_int',[900 500 930 530],'InitialCondition','Vbf0');

% ---- SC losses and power
addMLF(mdl,'PowerCalc',[660 650 820 730], strjoin({
'function [Ploss, Psc] = pcalc(iL, Vc, p)'
'Ploss = (p(26) + p(27))*iL^2 + Vc^2/p(30);'
'Psc   = Vc*iL;'
'end'}, newline));

% ---- energy counters
y0 = 400;  src = {'Battery/5', 'EMS/3', 'PowerCalc/1'};  nm = {'Eb','Efric','Eloss'};
for k = 1:3
    addB(mdl,'simulink/Math Operations/Gain',['Gain_' nm{k}],[980 y0+(k-1)*60 1020 y0+(k-1)*60+30],'Gain','1/3600');
    addB(mdl,'simulink/Continuous/Integrator',['Int_' nm{k}],[1050 y0+(k-1)*60 1090 y0+(k-1)*60+30],'InitialCondition','0');
    L(src{k}, ['Gain_' nm{k} '/1']);  L(['Gain_' nm{k} '/1'], ['Int_' nm{k} '/1']);
end

% ---- scopes
addB(mdl,'simulink/Sinks/Scope','Scope_Current',[1100 120 1140 170],'NumInputPorts','2');
addB(mdl,'simulink/Sinks/Scope','Scope_Vc',[1100 290 1140 330]);
addB(mdl,'simulink/Sinks/Scope','Scope_Power',[1180 420 1220 490],'NumInputPorts','3');
addB(mdl,'simulink/Sinks/Scope','Scope_SOC',[1180 520 1220 560]);

% ---- wiring
L('v_ref/1','EMS/1');   L('a_ref/1','EMS/2');   L('Vc_int/1','EMS/3');
L('Vbf_int/1','EMS/4'); L('SOC_int/1','EMS/5'); L('pvec/1','EMS/6');
L('EMS/2','ErrSum/1');  L('iL_int/1','ErrSum/2');
L('ErrSum/1','PID_Current/1');
L('PID_Current/1','Duty/1');  L('Vc_int/1','Duty/2');  L('Vbf_int/1','Duty/3');
L('Duty/1','PlantL/1');  L('Vc_int/1','PlantL/2');  L('iL_int/1','PlantL/3');
L('Vbf_int/1','PlantL/4');  L('pvec/1','PlantL/5');
L('PlantL/1','iL_int/1');
L('iL_int/1','PlantC/1');  L('Vc_int/1','PlantC/2');  L('pvec/1','PlantC/3');
L('PlantC/1','Vc_int/1');
L('EMS/1','Battery/1');  L('Duty/1','Battery/2');  L('iL_int/1','Battery/3');
L('Vbf_int/1','Battery/4');  L('SOC_int/1','Battery/5');  L('V1_int/1','Battery/6');
L('pvec/1','Battery/7');
L('Battery/1','SOC_int/1');  L('Battery/2','V1_int/1');  L('Battery/3','Vbf_int/1');
L('iL_int/1','PowerCalc/1');  L('Vc_int/1','PowerCalc/2');  L('pvec/1','PowerCalc/3');
L('EMS/2','Scope_Current/1');  L('iL_int/1','Scope_Current/2');
L('Vc_int/1','Scope_Vc/1');
L('Battery/4','Scope_Power/1');  L('PowerCalc/2','Scope_Power/2');  L('EMS/3','Scope_Power/3');
L('SOC_int/1','Scope_SOC/1');

% ---- logging (To Workspace blocks, timeseries format)
logs = {'out_v','v_ref',1; 'out_Vc','Vc_int',1; 'out_iL','iL_int',1; 'out_iref','EMS',2;
        'out_D','Duty',1; 'out_Pbus','EMS',1; 'out_Pfric','EMS',3;
        'out_Pbatt','Battery',4; 'out_Psc','PowerCalc',2;
        'out_Eb','Int_Eb',1; 'out_Efric','Int_Efric',1; 'out_Eloss','Int_Eloss',1;
        'out_soc','SOC_int',1; 'out_Ib','Battery',6; 'out_Vbus','Vbf_int',1};
for k = 1:size(logs,1)
    addB(mdl,'simulink/Sinks/To Workspace',logs{k,1},[1300 20+36*(k-1) 1380 46+36*(k-1)], ...
         'VariableName',logs{k,1},'SaveFormat','Timeseries','MaxDataPoints','inf');
    L(sprintf('%s/%d', logs{k,2}, logs{k,3}), [logs{k,1} '/1']);
end
save_system(mdl);
end

function blk = addB(mdl, lib, name, pos, varargin)
blk = [mdl '/' name];
add_block(lib, blk, 'Position', pos, varargin{:});
end

function addMLF(mdl, name, pos, code)
blk = [mdl '/' name];
add_block('simulink/User-Defined Functions/MATLAB Function', blk, 'Position', pos);
ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', blk);
ch.Script = code;
end

% ----------------------------------------------------------------------
%  Shared helpers (same motor / battery models as the Simulink blocks)
function B = baselineRun(P, cyc, regenOn)
% Reference run WITHOUT supercapacitor (same motor + battery models).
% regenOn = false : friction-only braking, true : battery-only regen.
Ne = numel(cyc.t);  dt = P.dt_e;
soc = P.soc0; v1 = 0; Vb = ocvPack(soc, P); Eb = 0;
Pb = zeros(Ne,1); Ib_l = Pb; Sc = Pb;
for k = 1:Ne
    Voc = ocvPack(soc, P);
    [Pw, ~, Pe] = motorPoint(cyc.v(k), cyc.a(k), P.m0, Vb, P);
    if Pw >= 0
        Pbus = Pe + P.Paux;
    else
        Pbc  = P.Pb_chg*min(max((0.98 - soc)/0.05, 0), 1);
        if regenOn, Pfin = min(Pe, Pbc); else, Pfin = 0; end
        Pbus = -Pfin + P.Paux;
    end
    E  = Voc - v1;
    Ib = (E - sqrt(max(E^2 - 4*P.R0*Pbus, 0)))/(2*P.R0);
    Vb = E - P.R0*Ib;
    soc = soc - Ib*dt/(3600*P.Q);
    v1  = v1 + dt*(Ib/P.C1 - v1/P.tau1);
    Eb  = Eb + Voc*Ib*dt;
    Pb(k) = Pbus; Ib_l(k) = Ib; Sc(k) = soc;
end
B.WhKm = Eb/3600/(trapz(cyc.t, cyc.v)/1000);
B.Pbatt = Pb; B.Ib = Ib_l; B.soc = Sc;
end

function [Pw, Psh, Pe] = motorPoint(v, a, m, Vb, P)
F  = m*a + P.k_aero*v^2 + P.Crr*m*P.g*(v > 0.05);
Pw = F*v;
wm = max(P.G*v/P.rw, 1e-3);
if Pw >= 0
    Tm  = F*P.rw/(P.G*P.eta_g);
    Psh = Tm*wm;
    if Psh > P.P_pk, Psh = P.P_pk; Tm = Psh/wm; end
    Tm  = min(Tm, P.T_pk);
    Pe  = Psh + motorLoss(Tm, wm, Vb, P);
else
    blend = min(max((v - P.v_lo)/(P.v_hi - P.v_lo), 0), 1);
    Tm  = -F*P.rw/P.G*P.eta_g*blend;
    Psh = Tm*wm;
    if Psh > P.P_pk, Psh = P.P_pk; Tm = Psh/wm; end
    if Tm > P.T_pk,  Tm = P.T_pk;  Psh = Tm*wm; end
    Pe  = max(Psh - motorLoss(Tm, wm, Vb, P), 0);
end
end

function Ploss = motorLoss(T, wm, Vb, P)
Iq  = abs(T)/P.Kt;
wnl = Vb/P.ke_LL;
Id  = P.Ifw*min(max((wm/wnl - 0.9)/0.5, 0), 1);
I2  = Iq^2 + Id^2;
Ploss = 1.5*P.Rs*I2 + P.kfe*wm^1.5 + P.kw*wm^2 ...
      + 1.5*P.R_inv*I2 + P.k_sw*sqrt(I2);
end

function V = ocvPack(soc, P)
V = P.Ns*interp1(P.soc_tab, P.ocv_tab, min(max(soc,0),1));
end

function Pw = wheelPower(v, a, m, P)
F  = m.*a + P.k_aero*v.^2 + P.Crr*m*P.g.*(v > 0.05);
Pw = F.*v;
end

function m = scMass(C, P),  m = 0.5*C*P.Vmax^2/(P.Whkg*3600); end
function c = scCost(C, P),  c = P.cost_kJ*0.5*C*P.Vmax^2/1e3; end

function g = tunePID(C, P)
Rt  = P.ESRC/C + P.RL;
Kp0 = 2*pi*150*P.L/P.Vbus;  Ki0 = Kp0*Rt/P.L;  Kd0 = Kp0*1e-4;
cost = @(th) pidCost(Kp0*exp(th(1)), Ki0*exp(th(2)), Kd0*exp(th(3)), Rt, P);
th = fminsearch(cost, [0 0 0], optimset('Display','off','MaxIter',150,'TolX',1e-3,'TolFun',1e-3));
g = [Kp0*exp(th(1)), Ki0*exp(th(2)), Kd0*exp(th(3))];
end

function J = pidCost(Kp, Ki, Kd, Rt, P)
[t, i] = stepSim([Kp Ki Kd], Rt, P);  ref = 100;
if any(~isfinite(i)) || max(abs(i)) > 1e3, J = 1e9; return; end
itae = sum(t.*abs((ref - i)/ref))*P.dt_c;
os   = max(max(i) - ref, 0)/ref;
J    = itae/1e-7 + 50*os^2;
end

function [t, i] = stepSim(g, Rt, P)
dt = P.dt_c; N = round(0.08/dt); Vc = P.Vc0; Vb = P.Vbus; ref = 100;
Kp = g(1); Ki = g(2); Kd = g(3); alpha = dt/(1e-3 + dt);
iL = 0; ui = 0; ud = 0; ip = 0; i = zeros(N,1); t = (1:N)'*dt;
for k = 1:N
    e = ref - iL;
    ud = ud + alpha*(-Kd*(iL - ip)/dt - ud);  ip = iL;
    uiT = ui + Ki*e*dt;  Dt = Vc/Vb + Kp*e + uiT + ud;
    if Dt > P.Dmin && Dt < P.Dmax, ui = uiT; D = Dt;
    else, D = min(max(Dt, P.Dmin), P.Dmax); end
    iL = iL + dt/P.L*(D*Vb - Vc - Rt*iL);
    i(k) = iL;
    if abs(iL) > 1e5, i(k:end) = NaN; break; end
end
end

function shadeBraking(t, Pw, ax)
% Light-green bands where wheel power < 0 (braking), drawn behind the curves.
isB = Pw(:) < 0;  d = diff([0; isB; 0]);
s = find(d == 1);  e = find(d == -1) - 1;
yl = ylim(ax);
for q = 1:numel(s)
    h = patch(ax, [t(s(q)) t(e(q)) t(e(q)) t(s(q))], [yl(1) yl(1) yl(2) yl(2)], ...
              [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    uistack(h, 'bottom');
end
ylim(ax, yl);
end
