% misinfo_numerics_matlab_equivalent.m
%
% MATLAB equivalent of the Python script misinfo_numerics.py.
% It generates the enriched 3-by-3 numerical study for the delayed
% multi-platform misinformation model with saturated incidence.
%
% The figure is saved automatically in:
%   ~/Documents/FakeNews_DDE_Matlab/misinfo_enriched.png
%
% No empirical data are used.
% The code uses fixed-step method-of-steps RK4 with linear interpolation
% of the stored history buffer, matching the Python implementation.
%
% NOTE:
% The n=20 random network data are hard-coded below from NumPy seed 7,
% because MATLAB and NumPy do not generate the same random numbers from
% the same seed. This makes the MATLAB script reproduce the same network
% instance used by the Python script.

clear; clc; close all;

%% Output folder
outDir = fullfile(getenv('HOME'), 'Documents', 'FakeNews_DDE_Matlab');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Plot style
set(groot, 'defaultAxesFontSize', 10);
set(groot, 'defaultAxesXGrid', 'on');
set(groot, 'defaultAxesYGrid', 'on');
set(groot, 'defaultAxesGridAlpha', 0.25);
set(groot, 'defaultAxesLayer', 'bottom');
set(groot, 'defaultAxesLineWidth', 0.8);

CB = [ ...
    hex2rgb('#0072B2'); ...
    hex2rgb('#D55E00'); ...
    hex2rgb('#009E73'); ...
    hex2rgb('#CC79A7'); ...
    hex2rgb('#E69F00'); ...
    hex2rgb('#56B4E9'); ...
    hex2rgb('#F0E442'); ...
    hex2rgb('#000000')  ...
];

%% Base 2-node parameters from the paper
Lam2 = [10.0; 8.0];
mu2  = [0.5; 0.4];
nu2  = [1.0; 0.8];
al2  = [0.3; 0.5];

tau2 = [2.0, 3.0;
        2.5, 1.5];

b_sub = [0.004,  0.001;
         0.0015, 0.003];

b_end = [0.05, 0.012;
         0.02, 0.045];

R0_sub_2 = spectral_R0(b_sub, Lam2, mu2, nu2);
R0_end_2 = spectral_R0(b_end, Lam2, mu2, nu2);

fprintf('2-node  R0(sub)=%.4f  R0(end)=%.4f\n', R0_sub_2, R0_end_2);

% Compute the positive endemic equilibrium.  The algebraic system also has
% the boundary equilibrium E0, so we do not use fsolve here; the monotone
% fixed-point iteration below selects the positive endemic fixed point when
% R0 > 1.
Ie = endemic_I_positive_fixed_point(b_end, Lam2, mu2, nu2, al2);
Qe = b_end * f_inc(Ie, al2);
Se = Lam2 ./ (mu2 + Qe);
fprintf('2-node  S* = [%g  %g]  I* = [%g  %g]\n', Se(1), Se(2), Ie(1), Ie(2));

%% Figure
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 13 9.5]);
tl = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Panel A: sub-threshold
nexttile;
[t, S, I] = integrate_dde(b_sub, Lam2, mu2, nu2, al2, tau2, 60, 0.02, [], [2.0; 2.0]);
plot(t, I(:,1), 'Color', CB(1,:), 'LineWidth', 1.1); hold on;
plot(t, I(:,2), 'Color', CB(2,:), 'LineWidth', 1.1);
title('(A) $n=2$, $\mathcal{R}_0\approx0.105\leq1$: $\to E_0$', 'Interpreter', 'latex', 'FontSize', 10);
xlabel('t');
ylabel('spreaders $I_i$', 'Interpreter', 'latex');
legend({'$I_1$', '$I_2$'}, 'Interpreter', 'latex', 'FontSize', 8, 'Location', 'northeast');
grid on;

%% Panel B: endemic
nexttile;
[t, S, I] = integrate_dde(b_end, Lam2, mu2, nu2, al2, tau2, 120, 0.02, [], [2.0; 0.05]);
plot(t, I(:,1), 'Color', CB(1,:), 'LineWidth', 1.1); hold on;
plot(t, I(:,2), 'Color', CB(2,:), 'LineWidth', 1.1);
yline(Ie(1), '--', 'Color', CB(1,:), 'LineWidth', 0.9);
yline(Ie(2), '--', 'Color', CB(2,:), 'LineWidth', 0.9);
title('(B) $n=2$, $\mathcal{R}_0\approx1.415>1$: $\to E^*$', 'Interpreter', 'latex', 'FontSize', 10);
xlabel('t');
ylabel('$I_i$', 'Interpreter', 'latex');
legend({'$I_1$', '$I_2$'}, 'Interpreter', 'latex', 'FontSize', 8, 'Location', 'northeast');
grid on;

%% Panel C: R0 heat-map over beta12, beta21
nexttile;
b_diag = 0.02;
b12 = linspace(0, 0.06, 140);
b21 = linspace(0, 0.08, 140);
R0grid = zeros(length(b21), length(b12));

for ii = 1:length(b21)
    for jj = 1:length(b12)
        B = [b_diag, b12(jj);
             b21(ii), b_diag];
        R0grid(ii,jj) = spectral_R0(B, Lam2, mu2, nu2);
    end
end

imagesc(b12, b21, R0grid);
set(gca, 'YDir', 'normal');
colormap(gca, viridis_map(256));
hold on;
contour(b12, b21, R0grid, [1.0, 1.0], 'w', 'LineWidth', 2.0);
text(0.043, 0.014, '$\mathcal{R}_0=1$', 'Interpreter', 'latex', ...
     'Color', 'w', 'FontSize', 8, 'Rotation', -20);
cb = colorbar;
cb.Label.String = 'R_0';
cb.Label.Interpreter = 'tex';
title('(C) $\mathcal{R}_0$ over cross-couplings', 'Interpreter', 'latex', 'FontSize', 10);
xlabel('$\beta_{12}$', 'Interpreter', 'latex');
ylabel('$\beta_{21}$', 'Interpreter', 'latex');
grid on;

%% Large network n = 20, hard-coded from NumPy seed 7
n = 20;
Lam = [6.4578497362437428, 10.679512753440687, 8.6304553886453608, 10.340791066985648, 11.867937071979616, 9.230975222462602, 9.0067227819596276, 6.4323068001585693, 7.6106338806112266, 8.9992950049533604, 10.075379976725642, 10.822434216626252, 8.2856467988912303, 6.3956180814354306, 7.7288735958479613, 11.457561166317682, 7.2803121214794935, 8.7127437709060978, 11.587236118134131, 6.1493953653020883]';
mu = [0.48016467523923678, 0.58503885012409373, 0.36909086370628941, 0.46454697577080911, 0.57273851246601937, 0.33995083372777501, 0.4570237742021297, 0.52512295773061046, 0.50070397226517416, 0.44032585792349421, 0.36145472708933851, 0.44722976672732112, 0.41171540681551766, 0.44322034645547648, 0.4097671157341779, 0.55137539829277815, 0.53059425195585275, 0.39419840316379862, 0.47178759979318619, 0.38281471449920851]';
nu = [0.92642146627320021, 0.87648918297219724, 1.0286997313898791, 0.88517554149401756, 0.92954648894571568, 1.0596620612545207, 0.90649591455691725, 1.1532116345821692, 0.79022580960134081, 1.070559436456632, 0.91118702182157008, 0.91322678634247112, 1.0171899343169193, 0.96145310051417265, 0.90744298921972155, 0.70071344028137905, 0.74613117292337594, 1.0546971968626062, 0.9621727983825985, 1.0480802317584845]';
al = [0.58218732920117033, 0.47316554175016701, 0.22125147626918304, 0.32354107394551884, 0.4370378749293492, 0.29404816290298574, 0.58598839981445072, 0.57801928951711745, 0.5393603523350754, 0.38892959851536091, 0.53659068595935888, 0.25244425693916556, 0.32349346291913417, 0.38519855766176681, 0.4967388802733867, 0.39433009148355869, 0.25475044751898035, 0.33741461188174332, 0.32977046786897718, 0.32016756172721589]';
A = [
    0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1;
    1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0;
    0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1;
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1;
    1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0;
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0;
    0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0
];
tau = [
    0, 0, 0, 2.2558510488333594, 0, 0, 0, 2.4321606161342983, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.617638102550174, 0, 2.3998851856947776;
    3.6970029736879542, 0, 2.8318501130060159, 0, 0, 0, 0, 0, 0, 0, 0, 3.5736873231582598, 0, 2.479949576917539, 0, 0, 0, 0, 0, 0;
    0, 1.2180951983688857, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.1903409504058575, 0, 0, 0, 1.2928709915434815, 0, 3.9631360015094854, 0, 0;
    0, 0, 2.0012583272487454, 0, 2.0025643231413399, 0, 0, 0, 0, 0, 0, 0, 1.2600538651912365, 0, 0, 0, 2.511913461038211, 0, 0, 0;
    0, 0, 0, 3.2178826971135339, 0, 3.638450240140346, 3.9945631770053063, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.6393357661515267, 0, 0;
    0, 0, 0, 0, 3.3169567843306265, 0, 1.9731169249567806, 0, 0, 0, 2.8560387989201863, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1.8004149727619758, 0, 0, 0, 1.4082681422282386, 0, 0, 0, 0, 0, 0, 2.8654847760999269, 0, 0, 1.4445001535271045;
    0, 0, 0, 0, 0, 3.1900684987362653, 3.3340494132976501, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.897789539216574, 1.5507285114105398;
    3.3602327068298981, 0, 0, 2.9549570304811512, 0, 0, 0, 2.1630537585361913, 0, 0, 0, 0, 0, 0, 3.7922987892271594, 0, 0, 0, 0, 0;
    0, 0, 0, 2.9165635910112107, 0, 0, 0, 0, 1.0000195713446485, 0, 0, 0, 2.2080021195696911, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 2.681565918829345, 3.9804174971577324, 0, 0, 0, 3.2997136691794591, 0, 0, 1.7904696339569584, 0, 0, 0;
    3.4754391154794897, 2.4904916000724229, 0, 0, 0, 0, 0, 0, 0, 0, 2.9948646184048409, 0, 0, 0, 0, 0, 0, 1.3387183891618586, 0, 0;
    0, 0, 0, 2.3747741074584106, 0, 0, 2.3632964953830422, 0, 0, 0, 0, 2.0393680700179551, 0, 3.4970215623000698, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.7905716630175426, 0, 1.2926498093687533, 0, 2.0735903683614416, 3.6024596655366201, 0, 0, 0, 0;
    0, 1.7692553001963354, 0, 0, 0, 0, 0, 2.7739379558441173, 0, 0, 0, 0, 0, 1.4617869910952539, 0, 0, 1.6777988074154215, 0, 0, 0;
    0, 1.809042857597269, 0, 0, 0, 0, 2.9519971547205461, 0, 0, 0, 0, 0, 0, 0, 3.2199614571416695, 0, 0, 0, 0, 1.9223489133021139;
    0, 0, 0, 0, 0, 2.053365485583138, 0, 0, 0, 0, 0, 0, 0, 0, 1.1741174438185658, 1.3586166468396859, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 2.9671053310457678, 0, 0, 1.1027866838182074, 3.1009657863799882, 0, 0, 0, 0, 0, 0, 3.4304780836910886, 0, 0, 0;
    1.6984717114622794, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.7765218799013849, 1.7943899680223951, 0, 0, 0, 0, 0, 2.9011798809777942, 0, 0;
    0, 3.9030227658829739, 3.7740108937663379, 0, 0, 0, 0, 0, 0, 1.5580540090700667, 0, 0, 0, 0, 0, 0, 0, 0, 1.0025620999360378, 0
];
b_base = [
    0, 0, 0, 0.015022414373617555, 0, 0, 0, 0.015934530439424416, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.017034889902616972, 0, 0.013396654952870065;
    0.007962758517279022, 0, 0.0050387682786651108, 0, 0, 0, 0, 0, 0, 0, 0, 0.0039377641449612512, 0, 0.018511801606783951, 0, 0, 0, 0, 0, 0;
    0, 0.017877858448971204, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.012950354130923894, 0, 0, 0, 0.0070457515372543446, 0, 0.01570542049067088, 0, 0;
    0, 0, 0.0057012523199033967, 0, 0.0073868975925606636, 0, 0, 0, 0, 0, 0, 0, 0.007462779852718442, 0, 0, 0, 0.014638741168660565, 0, 0, 0;
    0, 0, 0, 0.015079176832643008, 0, 0.0092138924528625508, 0.016429967057223868, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.018255237063706437, 0, 0;
    0, 0, 0, 0, 0.016340770532387312, 0, 0.0073889511755181282, 0, 0, 0, 0.017696876387570462, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0.0066828265743270679, 0, 0, 0, 0.0086758322215206157, 0, 0, 0, 0, 0, 0, 0.014622618488383029, 0, 0, 0.01328383546564791;
    0, 0, 0, 0, 0, 0.0083244432466210577, 0.013986651466502834, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0083574475930036866, 0.012564099938598004;
    0.0086309622544055269, 0, 0, 0.014999162077392584, 0, 0, 0, 0.01609793635283064, 0, 0, 0, 0, 0, 0, 0.019713203167988864, 0, 0, 0, 0, 0;
    0, 0, 0, 0.0058704099592041519, 0, 0, 0, 0, 0.011782080576444868, 0, 0, 0, 0.0093725262646891742, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0.013524131629240238, 0.0083870686780910773, 0, 0, 0, 0.019827859143537734, 0, 0, 0.0051677191094348599, 0, 0, 0;
    0.011964288869625949, 0.010626785743388786, 0, 0, 0, 0, 0, 0, 0, 0, 0.018064959558474314, 0, 0, 0, 0, 0, 0, 0.006961260317894344, 0, 0;
    0, 0, 0, 0.016376753651777647, 0, 0, 0.017578913218135386, 0, 0, 0, 0, 0.018617190092451136, 0, 0.009379846751887631, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.019889207518631394, 0, 0.013325525277193617, 0, 0.019536295832395384, 0.00426947370955107, 0, 0, 0, 0;
    0, 0.016368668595160524, 0, 0, 0, 0, 0, 0.0029952733853831313, 0, 0, 0, 0, 0, 0.0045281439477409518, 0, 0, 0.019564601919532507, 0, 0, 0;
    0, 0.014565615094193428, 0, 0, 0, 0, 0.0075034504804033231, 0, 0, 0, 0, 0, 0, 0, 0.014023391415530211, 0, 0, 0, 0, 0.011487780268775458;
    0, 0, 0, 0, 0, 0.01180517984349372, 0, 0, 0, 0, 0, 0, 0, 0, 0.016349630906982382, 0.019498705574206479, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0.018420546293570754, 0, 0, 0.018854284148730742, 0.0059289622274240816, 0, 0, 0, 0, 0, 0, 0.0079665422427009271, 0, 0, 0;
    0.012239809257748795, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0080107263392258776, 0.0084798532035866691, 0, 0, 0, 0, 0, 0.0089439189083262011, 0, 0;
    0, 0.016901674034575178, 0.0056786989861845501, 0, 0, 0, 0, 0, 0, 0.007547626403741718, 0, 0, 0, 0, 0, 0, 0, 0, 0.014097955912988408, 0
];
I0D = [2.6105119315829084, 1.6368274380645489, 1.0941798632841593, 1.9325889024687695, 1.66597486403814, 1.4946179921935738, 3.6833150312858045, 2.0384768030875833, 2.5138794382512746, 1.84395126564899, 2.3538226256387182, 3.2329229752150574, 3.6993688075576441, 3.5618474837033745, 3.3442736436284473, 1.7256046624077126, 2.0930093323769112, 3.3762432135978448, 3.2902901834046392, 1.888243732721107]';
I0E = [1.1054932202916994, 2.861375211373852, 1.584277440161878, 2.0300511619686654, 2.8760840446437617, 2.8051577979308524, 0.64357976435275555, 1.2377527949003415, 2.8451042950832233, 2.1352779184704076, 2.3751429376676469, 2.493886826157691, 1.112002393129998, 1.692930214735628, 2.9087869648715619, 2.2602837231856689, 1.3716301758603817, 0.8787594138690682, 0.4119583398199162, 1.4532922406586235]';

R0_base = spectral_R0(b_base, Lam, mu, nu);
s_sub = 0.15 / max(R0_base, 1e-9);
s_sup = 2.2  / max(R0_base, 1e-9);

b_net_sub = s_sub * b_base;
b_net_sup = s_sup * b_base;

R0sub = spectral_R0(b_net_sub, Lam, mu, nu);
R0sup = spectral_R0(b_net_sup, Lam, mu, nu);

fprintf('n=20  R0(sub)=%.3f  R0(sup)=%.3f\n', R0sub, R0sup);

%% Panel D: n=20 subcritical
nexttile;
[t, S, I] = integrate_dde(b_net_sub, Lam, mu, nu, al, tau, 80, 0.03, [], I0D);
for i = 1:n
    cidx = mod(i-1, 7) + 1;
    plot(t, I(:,i), 'Color', CB(cidx,:), 'LineWidth', 0.8); hold on;
end
title(sprintf('(D) $n=20$ network, $\\mathcal{R}_0\\approx%.2f\\leq1$', R0sub), ...
      'Interpreter', 'latex', 'FontSize', 10);
xlabel('t');
ylabel('$I_i$ (all platforms)', 'Interpreter', 'latex');
grid on;

%% Panel E: n=20 supercritical
nexttile;
[t, S, I] = integrate_dde(b_net_sup, Lam, mu, nu, al, tau, 140, 0.03, [], I0E);
for i = 1:n
    cidx = mod(i-1, 7) + 1;
    plot(t, I(:,i), 'Color', CB(cidx,:), 'LineWidth', 0.8); hold on;
end
title(sprintf('(E) $n=20$ network, $\\mathcal{R}_0\\approx%.2f>1$: persistence', R0sup), ...
      'Interpreter', 'latex', 'FontSize', 10);
xlabel('t');
ylabel('$I_i$', 'Interpreter', 'latex');
grid on;

%% Panel F: alpha sweep at fixed beta
nexttile;
alphas = [0.05, 0.3, 0.8, 1.5];
for c = 1:length(alphas)
    alv = [alphas(c); alphas(c)];
    [t, S, I] = integrate_dde(b_end, Lam2, mu2, nu2, alv, tau2, 140, 0.02, [], [2.0; 0.05]);
    plot(t, I(:,1), 'Color', CB(c,:), 'LineWidth', 1.1); hold on;
end
title('(F) Vary $\alpha$ at fixed $\mathcal{R}_0=1.415$: peak/level shift, threshold fixed', ...
      'Interpreter', 'latex', 'FontSize', 9);
xlabel('t');
ylabel('$I_1$', 'Interpreter', 'latex');
legend({'$\alpha=0.05$', '$\alpha=0.30$', '$\alpha=0.80$', '$\alpha=1.50$'}, ...
       'Interpreter', 'latex', 'FontSize', 7, 'Location', 'northeast');
grid on;

%% Panel G: delay sweep at fixed beta
nexttile;
delay_scales = [0.0, 1.0, 2.0, 3.5];
for c = 1:length(delay_scales)
    scale = delay_scales(c);
    [t, S, I] = integrate_dde(b_end, Lam2, mu2, nu2, al2, scale*tau2, 160, 0.02, [], [3.0; 0.02]);
    plot(t, I(:,1), 'Color', CB(c,:), 'LineWidth', 1.1); hold on;
end
yline(Ie(1), '--', 'Color', 'k', 'LineWidth', 0.8);
title('(G) Vary delays at fixed $\beta$: same $E^*$, slower/overshoot', ...
      'Interpreter', 'latex', 'FontSize', 9);
xlabel('t');
ylabel('$I_1$', 'Interpreter', 'latex');
legend({'$\tau\times0.0$', '$\tau\times1.0$', '$\tau\times2.0$', '$\tau\times3.5$'}, ...
       'Interpreter', 'latex', 'FontSize', 7, 'Location', 'northeast');
grid on;

%% Panel H: sensitivity tornado of R0
nexttile;
R0b = spectral_R0(b_end, Lam2, mu2, nu2);

hi_beta = spectral_R0(1.2*b_end, Lam2, mu2, nu2) - R0b;
lo_beta = spectral_R0(0.8*b_end, Lam2, mu2, nu2) - R0b;

hi_Lam = spectral_R0(b_end, 1.2*Lam2, mu2, nu2) - R0b;
lo_Lam = spectral_R0(b_end, 0.8*Lam2, mu2, nu2) - R0b;

hi_mu = spectral_R0(b_end, Lam2, 1.2*mu2, nu2) - R0b;
lo_mu = spectral_R0(b_end, Lam2, 0.8*mu2, nu2) - R0b;

hi_nu = spectral_R0(b_end, Lam2, mu2, 1.2*nu2) - R0b;
lo_nu = spectral_R0(b_end, Lam2, mu2, 0.8*nu2) - R0b;

his = [hi_beta; hi_Lam; hi_mu; hi_nu];
los = [lo_beta; lo_Lam; lo_mu; lo_nu];

y = 1:4;
barh(y, his, 'FaceColor', CB(1,:), 'FaceAlpha', 0.8); hold on;
barh(y, los, 'FaceColor', CB(2,:), 'FaceAlpha', 0.8);
xline(0, 'k', 'LineWidth', 0.8);
set(gca, 'YTick', y, 'YTickLabel', {'$\beta$', '$\Lambda$', '$\mu$', '$\nu$'}, ...
         'TickLabelInterpreter', 'latex');
title('(H) Sensitivity of $\mathcal{R}_0$ ($\pm20\%$), base $=1.415$', ...
      'Interpreter', 'latex', 'FontSize', 9);
xlabel('$\Delta\mathcal{R}_0$', 'Interpreter', 'latex');
legend({'+20\%', '-20\%'}, 'Interpreter', 'latex', 'FontSize', 7, 'Location', 'northeast');
grid on;

%% Panel I: moderation intervention on platform 1 only
% Increasing only nu_1 lowers R0 and I_1^*, but it cannot eradicate
% misinformation in this parameter set because K_22 = 1.125 > 1.
nexttile;
factors = linspace(1.0, 5.0, 35);
R0s = zeros(size(factors));
I1s = zeros(size(factors));
I2s = zeros(size(factors));
K22 = b_end(2,2) * Lam2(2) / (mu2(2) * nu2(2));

for k = 1:length(factors)
    fmul = factors(k);
    nu_new = nu2 .* [fmul; 1.0];

    R0s(k) = spectral_R0(b_end, Lam2, mu2, nu_new);

    % Positive endemic equilibrium, not the boundary equilibrium E0.
    Ix = endemic_I_positive_fixed_point(b_end, Lam2, mu2, nu_new, al2);
    I1s(k) = Ix(1);
    I2s(k) = Ix(2);
end

yyaxis left;
hR0 = plot(factors, R0s, 'Color', CB(3,:), 'LineWidth', 1.6); hold on;
hOne = yline(1, '--', 'Color', 'k', 'LineWidth', 0.8);
hK22 = yline(K22, ':', 'Color', CB(3,:), 'LineWidth', 0.9);
ylabel('$\mathcal{R}_0$', 'Interpreter', 'latex');
ylim([1.0, max(R0s)+0.08]);

yyaxis right;
hI1 = plot(factors, I1s, 'Color', CB(1,:), 'LineWidth', 1.2); hold on;
hI2 = plot(factors, I2s, 'Color', CB(2,:), 'LineWidth', 1.2);
ylabel('endemic $I_i^*$', 'Interpreter', 'latex');
ylim([0, max([I1s(:); I2s(:)])*1.12]);

xlabel('factor multiplying $\nu_1$', 'Interpreter', 'latex');
title('(I) Raising $\nu_1$: $\mathcal{R}_0$ decreases, no eradication ($\mathcal{K}_{22}>1$)', ...
      'Interpreter', 'latex', 'FontSize', 8.5);
legend([hR0, hOne, hK22, hI1, hI2], ...
       {'$\mathcal{R}_0$', '$\mathcal{R}_0=1$', '$\mathcal{K}_{22}=1.125$', '$I_1^*$', '$I_2^*$'}, ...
       'Interpreter', 'latex', 'FontSize', 6.5, 'Location', 'east');
grid on;

%% Overall title and save
title(tl, 'Enriched numerical study -- delayed multi-platform misinformation model (all delays $\tau_{ij}$ active)', ...
      'Interpreter', 'latex', 'FontSize', 12);

outFile = fullfile(outDir, 'misinfo_enriched.png');
try
    exportgraphics(fig, outFile, 'Resolution', 200);
catch
    print(fig, outFile, '-dpng', '-r200');
end

fprintf('saved %s\n', outFile);

%% ============================================================
% Local functions
% ============================================================

function value = f_inc(u, alpha)
    value = u ./ (1.0 + alpha .* u);
end

function R0 = spectral_R0(beta, Lam, mu, nu)
    rowScale = Lam ./ (mu .* nu);
    K = bsxfun(@times, beta, rowScale);
    R0 = max(abs(eig(K)));
end

function [Sstar, Istar, ier] = endemic_equilibrium(beta, Lam, mu, nu, alpha, guess)
    %#ok<INUSD>
    % Robust equilibrium selector.
    % If R0 <= 1, the only biologically relevant equilibrium is E0.
    % If R0 > 1, the algebraic equations also admit E0, so a generic solver
    % may converge to the boundary.  We therefore compute the positive
    % endemic equilibrium by monotone fixed-point iteration in I.
    R0 = spectral_R0(beta, Lam, mu, nu);

    if R0 <= 1
        Istar = zeros(length(Lam), 1);
        Sstar = Lam ./ mu;
        ier = true;
        return;
    end

    Istar = endemic_I_positive_fixed_point(beta, Lam, mu, nu, alpha);
    Q = beta * f_inc(Istar, alpha);
    Sstar = Lam ./ (mu + Q);
    ier = all(Istar > 0) && all(isfinite(Istar)) && all(isfinite(Sstar));
end

function Istar = endemic_I_positive_fixed_point(beta, Lam, mu, nu, alpha)
    % Computes the positive endemic fixed point of
    %   I_i = Lambda_i Q_i(I)/(nu_i (mu_i + Q_i(I))),
    % where Q_i(I)=sum_j beta_ij f_j(I_j).
    % Starting from the positive upper vector Lam./nu avoids convergence to
    % the boundary equilibrium I=0 in the supercritical case R0>1.
    n = length(Lam);
    Iold = Lam ./ nu;
    tol = 1e-12;
    maxIter = 100000;

    for it = 1:maxIter
        Q = beta * f_inc(Iold, alpha);
        Inew = Lam .* Q ./ (nu .* (mu + Q));

        if max(abs(Inew - Iold)) < tol
            Istar = max(Inew, 0);
            return;
        end

        % Damping improves robustness for near-threshold cases.
        Iold = 0.5 * Iold + 0.5 * Inew;
    end

    Istar = max(Iold, 0);
end

function F = endemic_residual(v, beta, Lam, mu, nu, alpha)
    n = length(Lam);
    S = v(1:n);
    I = v(n+1:2*n);

    inc = beta * f_inc(I, alpha);
    F = [Lam - mu .* S - S .* inc;
         S .* inc - nu .* I];
end

function [t, S, I] = integrate_dde(beta, Lam, mu, nu, alpha, tau, T, dt, S0, I0)
    n = length(Lam);
    nT = round(T / dt);

    if nargin < 9 || isempty(S0)
        S0 = Lam ./ mu;
    end
    if nargin < 10 || isempty(I0)
        I0 = 0.15 * ones(n,1);
    end

    S0 = S0(:);
    I0 = I0(:);

    taumax = max(tau(:));
    hist = ceil(taumax / dt) + 2;

    S = zeros(nT + 1, n);
    I = zeros(nT + 1, n);

    Ibuf = repmat(I0.', hist, 1);
    Ihist = [Ibuf; zeros(nT + 1, n)];

    S(1,:) = S0.';
    I(1,:) = I0.';
    Ihist(hist + 1,:) = I0.';

    t = linspace(0, T, nT + 1).';

    for kk = 1:nT
        k = kk - 1;

        Sk = S(kk,:).';
        Ik = I(kk,:).';

        [dS1, dI1] = rhs_dde(k,       Sk,                  Ik,                  k, beta, Lam, mu, nu, alpha, tau, dt, hist, nT, Ihist);
        [dS2, dI2] = rhs_dde(k + 0.5, Sk + 0.5*dt*dS1,     Ik + 0.5*dt*dI1,     k, beta, Lam, mu, nu, alpha, tau, dt, hist, nT, Ihist);
        [dS3, dI3] = rhs_dde(k + 0.5, Sk + 0.5*dt*dS2,     Ik + 0.5*dt*dI2,     k, beta, Lam, mu, nu, alpha, tau, dt, hist, nT, Ihist);
        [dS4, dI4] = rhs_dde(k + 1.0, Sk + dt*dS3,         Ik + dt*dI3,         k, beta, Lam, mu, nu, alpha, tau, dt, hist, nT, Ihist);

        S_next = Sk + dt/6 * (dS1 + 2*dS2 + 2*dS3 + dS4);
        I_next = Ik + dt/6 * (dI1 + 2*dI2 + 2*dI3 + dI4);

        I_next = max(I_next, 0.0);

        S(kk+1,:) = S_next.';
        I(kk+1,:) = I_next.';
        Ihist(hist + kk + 1,:) = I_next.';
    end
end

function [dS, dI] = rhs_dde(k_cont, Sv, Iv, kfill, beta, Lam, mu, nu, alpha, tau, dt, hist, nT, Ihist)
    n = length(Lam);
    kdelay = k_cont - tau / dt;
    inc = zeros(n,1);

    for i = 1:n
        Ijdel = delayed_I_vec(kdelay(i,:), kfill, hist, nT, Ihist);
        inc(i) = sum(beta(i,:).' .* f_inc(Ijdel, alpha));
    end

    dS = Lam - mu .* Sv - Sv .* inc;
    dI = Sv .* inc - nu .* Iv;
end

function Ij = delayed_I_vec(kfloats, kfill, hist, nT, Ihist)
    n = length(kfloats);
    Ij = zeros(n,1);

    for j = 1:n
        idx = min(hist + kfloats(j), hist + kfill);
        i0 = floor(idx);
        frac = idx - i0;

        i0 = max(0, min(i0, hist + nT - 1));
        i1 = max(0, min(i0 + 1, hist + nT - 1));

        row0 = i0 + 1;
        row1 = i1 + 1;

        Ij(j) = (1 - frac) * Ihist(row0, j) + frac * Ihist(row1, j);
    end
end

function rgb = hex2rgb(hex)
    if hex(1) == '#'
        hex = hex(2:end);
    end
    rgb = reshape(sscanf(hex, '%2x'), 1, 3) / 255;
end

function cmap = viridis_map(m)
    if nargin < 1
        m = 256;
    end

    base = [
         68,   1,  84
         72,  35, 116
         64,  67, 135
         52,  94, 141
         41, 120, 142
         32, 144, 140
         34, 167, 132
         68, 190, 112
        121, 209,  81
        189, 223,  38
        253, 231,  37
    ] / 255;

    x = linspace(0, 1, size(base,1));
    xi = linspace(0, 1, m);

    cmap = zeros(m,3);
    for k = 1:3
        cmap(:,k) = interp1(x, base(:,k), xi, 'pchip');
    end

    cmap = max(min(cmap, 1), 0);
end
