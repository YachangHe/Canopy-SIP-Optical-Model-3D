% main_SIP_optical
% User-facing entry for the SIP optical reflectance model.
%
% Coordinate convention:
%   SZA and VZA are zenith angles in degrees from the canopy normal.
%   0 deg is vertical. Relative azimuth is abs(wrapTo180(VAA - SAA)).

clear
clc

cfg = struct();

%% 1. Common inputs for both 1D and 3D scenes
cfg.scene_type = '1D';             % '1D' or '3D'

% Sun-view geometry
cfg.geometry.SZA = 20;             % Solar zenith angle [deg]
cfg.geometry.SAA = 0;              % Solar azimuth angle [deg]
cfg.geometry.view_angles = [];     % Optional [VZA, VAA] grid [deg]; [] uses bundled grid

% Leaf optical properties for PROSPECT-D
cfg.leaf.N = 1.5;                  % Leaf structure parameter [-]
cfg.leaf.Cab = 30.0;               % Chlorophyll a+b [ug cm-2]
cfg.leaf.Car = 8.0;                % Carotenoids [ug cm-2]
cfg.leaf.Ant = 0.0;                % Anthocyanins [ug cm-2]
cfg.leaf.Brown = 0.0;              % Brown pigment [-]
cfg.leaf.Cw = 0.015;               % Equivalent water thickness [cm]
cfg.leaf.Cm = 0.012;               % Dry matter [g cm-2]

% Leaf angle distribution
cfg.lad.TypeLidf = 2;              % 1 = Verhoef, 2 = Campbell
cfg.lad.LIDFa = 60;                % Campbell mean leaf angle [deg]
cfg.lad.LIDFb = -0.15;             % Used only for TypeLidf = 1

% Soil and output controls
cfg.soil.scale = 1.0;              % 1 = measured soil, 0 = black soil
cfg.output.save_output = true;

%% 2. Inputs used only by the 1D homogeneous/turbid-medium scene
cfg.turbid.LAI = 3.0;              % Scene LAI [m2 m-2]
cfg.turbid.Omega = 1.0;            % Macroscopic clumping factor [-]
cfg.turbid.hotspot = 0.05;         % Hotspot size parameter [-], must be > 0

%% 3. Inputs used only by the 3D discrete crown scene
cfg.canopy.LAI_crown = 5.0;        % Single-crown LAI [m2 m-2]
cfg.canopy.Height = 1.0;           % Canopy top height [m]
cfg.canopy.Crowndeepth = 0.66;     % Mean within-crown path depth [m]
cfg.canopy.lmax_nadir = 1.0;       % Nadir maximum crown path [m]
cfg.canopy.leaf_diameter = 0.01;   % Crown-scale hotspot leaf-size parameter [m]
cfg.canopy.use_external_crown_center = true;
cfg.canopy.external_crown_center = 0.5;
cfg.numerics.CI_min = 0.01;        % 3D effective CI lower bound
cfg.numerics.CI_max = 2.50;        % 3D effective CI upper bound

[outputs, diagnostics] = sip_optical_core(cfg);
