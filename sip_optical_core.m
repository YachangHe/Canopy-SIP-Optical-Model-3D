function [outputs, diagnostics] = sip_optical_core(user_cfg)
%SIP_OPTICAL_CORE SIP canopy optical model using GO and spectral invariants.
%   Combines geometric-optical canopy structure with spectral-invariant multiple scattering.
%
% Coordinate convention / 坐标约定:
%   All angular inputs are in degrees. Zenith angles are measured from the
%   upward canopy normal; 0 deg is vertical. Relative azimuth is computed as
%   abs(wrapTo180(VAA - SAA)).
%   所有角度均为 degree；天顶角从冠层法线起算，0 度表示垂直方向；
%   相对方位角在模型内部定义为 abs(wrapTo180(VAA - SAA))。

if nargin < 1
    user_cfg = struct();
end

cfg = merge_config(default_config(), user_cfg);
validate_config(cfg);

model_dir = fileparts(mfilename('fullpath'));
addpath(model_dir);

gap_dir = fullfile(model_dir, 'data', 'CI_HET10');
soil_file = fullfile(model_dir, 'data', 'soilnew_1_to_LESS.txt');
scene_type = normalize_scene_type(cfg.scene_type);

gap_tot = load_gap(fullfile(gap_dir, 'A_gap_tot_HET10.mat'), 'gap_tot');
gap_within = load_gap(fullfile(gap_dir, 'A_gap_whithin_HET10.mat'), 'gap_within');
gap_betw = load_gap(fullfile(gap_dir, 'A_gap_betw_HET10.mat'), 'gap_betw');
[va, gap_within_cond, gap_encoding] = prepare_gap_tables(gap_tot, gap_within, gap_betw);

SZA = cfg.geometry.SZA;
SAA = cfg.geometry.SAA;
if strcmp(scene_type, '1D')
    LAI_crown = cfg.turbid.LAI;
else
    LAI_crown = cfg.canopy.LAI_crown;
end
Height = cfg.canopy.Height;
Crowndeepth = cfg.canopy.Crowndeepth;
lmax_nadir = cfg.canopy.lmax_nadir;
leaf_diameter = cfg.canopy.leaf_diameter;

sun_idx = find_gap_row(gap_tot, SZA, SAA);
nadir_idx = find_gap_row(gap_tot, 0, 0);
Pw_s = gap_betw(sun_idx, 3);
Ptot_s = gap_tot(sun_idx, 3);
Pi_s = gap_within_cond(sun_idx);
Pw_h = gap_betw(nadir_idx, 3);

[Themi_gap_gauss, Themi_gap_trapz] = compute_hemispherical_gap(gap_tot);
iD = 1 - Themi_gap_gauss;
assert_probability(iD, 'hemispherical interceptance iD');
crown_cover_nadir = 1 - Pw_h;

% Scene LAI for p-theory / 用于 p-theory 场景 LAI
crown_cover_nadir = 1 - Pw_h;
LAI_scene = LAI_crown * crown_cover_nadir;
require_scalar_range(LAI_scene, eps, Inf, 'Derived scene LAI');

[wave_p, rho_all, tau_all, w_all] = prepare_leaf_optics(cfg.leaf);
nwl = numel(wave_p);
rg_all = prepare_soil_reflectance(soil_file, nwl, cfg.soil.scale);
lidf = prepare_lidf(cfg.lad);

if strcmp(scene_type, '3D')
    [Height_c, Height_c_path, continuity_c, HotSpotPar, go_par] = prepare_path_structure( ...
        cfg, Pw_s, Pw_h);
else
    va = prepare_view_angles(cfg.geometry.view_angles, va);
    [gap_tot, gap_betw, gap_within_cond, Ptot_s, Pi_s, iD, Themi_gap_gauss, Themi_gap_trapz] = ...
        prepare_1d_gap_terms(va, SZA, SAA, lidf, LAI_crown, cfg.turbid.Omega);
    gap_encoding = '1d_homogeneous';
    Pw_s = 0;
    Pw_h = 0;
    crown_cover_nadir = 1;
    LAI_scene = LAI_crown;
    Height_c = NaN;
    Height_c_path = NaN;
    continuity_c = NaN;
    HotSpotPar = cfg.turbid.hotspot;
    go_par = NaN;
end

p = 1 - iD / LAI_scene;
require_scalar_range(p, 0, 1 - eps, 'Recollision probability p');
if any(1 - p .* w_all <= 0)
    error('Energy denominator 1 - p*w must be positive for every wavelength.');
end

nang = size(va, 1);
terms = allocate_outputs(nang, nwl);

for t = 1:nang
    view = compute_view_geometry( ...
        t, va, SZA, SAA, gap_tot, gap_betw, gap_within_cond, ...
        Pw_s, Pi_s, LAI_crown, lidf, Height, Height_c, Crowndeepth, ...
        continuity_c, HotSpotPar, go_par, cfg.numerics.CI_min, cfg.numerics.CI_max, ...
        scene_type, cfg.turbid.Omega);

    [terms, row] = solve_optical_row( ...
        terms, t, view, rho_all, tau_all, w_all, rg_all, ...
        Ptot_s, iD, LAI_scene, p, scene_type);

    terms.K_components(t, :) = row.K_components;
    terms.CI_terms(t, :) = row.CI_terms;
    terms.Pi_terms(t, :) = row.Pi_terms;
    terms.leaf_terms(t, :) = row.leaf_terms;
    terms.phase_terms(t, :) = row.phase_terms;
end

outputs = struct();
outputs.wavelength = wave_p;
outputs.view_angles = va;
outputs.BRF_total = terms.BRF_total;
outputs.BRF_single = terms.BRF_single;
outputs.BRF_veg_single = terms.BRF_veg_single;
outputs.BRF_vegetation = terms.BRF_vegetation;
outputs.BRF_soil = terms.BRF_soil;
outputs.BRF_veg_soil_interaction = terms.BRF_veg_soil_interaction;
outputs.components = struct( ...
    'BRF_leaf_C', terms.BRF_leaf_C, ...
    'BRF_leaf_T', terms.BRF_leaf_T, ...
    'BRF_soil_crown', terms.BRF_soil_crown, ...
    'BRF_soil_open', terms.BRF_soil_open, ...
    'BRF_soil_turbid', terms.BRF_soil_turbid, ...
    'BRF_CM', terms.BRF_CM, ...
    'BRF_GCM', terms.BRF_GCM);

diagnostics = struct();
diagnostics.config = cfg;
diagnostics.scene_type = scene_type;
diagnostics.turbid_Omega = cfg.turbid.Omega;
diagnostics.gap_encoding = gap_encoding;
diagnostics.Pw_s = Pw_s;
diagnostics.Ptot_s = Ptot_s;
diagnostics.Pi_s = Pi_s;
diagnostics.Pw_h = Pw_h;
diagnostics.Themi_gap_gauss = Themi_gap_gauss;
diagnostics.Themi_gap_trapz = Themi_gap_trapz;
diagnostics.crown_cover_nadir = crown_cover_nadir;
diagnostics.LAI_scene = LAI_scene;
diagnostics.iD = iD;
diagnostics.p = p;
diagnostics.Height_c = Height_c;
diagnostics.Height_c_path = Height_c_path;
diagnostics.continuity_c = continuity_c;
diagnostics.HotSpotPar = HotSpotPar;
diagnostics.go_par = go_par;
diagnostics.K_components = terms.K_components;
diagnostics.CI_terms = terms.CI_terms;
diagnostics.Pi_terms = terms.Pi_terms;
diagnostics.leaf_terms = terms.leaf_terms;
diagnostics.phase_terms = terms.phase_terms;

if cfg.output.save_output
    save(cfg.output.output_file, 'outputs', 'diagnostics');
end
end

function cfg = default_config()
model_dir = fileparts(mfilename('fullpath'));
cfg = struct();
cfg.scene_type = '3D';
cfg.geometry = struct('SZA', 20, 'SAA', 0, 'view_angles', []);
cfg.canopy = struct( ...
    'LAI_crown', 5.0, ...
    'Height', 1.0, ...
    'Crowndeepth', 0.66, ...
    'lmax_nadir', 1.0, ...
    'leaf_diameter', 0.01, ...
    'use_external_crown_center', true, ...
    'external_crown_center', 0.5);
cfg.leaf = struct('N', 1.5, 'Cab', 30, 'Car', 8, 'Ant', 0, ...
    'Brown', 0, 'Cw', 0.015, 'Cm', 0.012);
cfg.lad = struct('TypeLidf', 2, 'LIDFa', 60, 'LIDFb', -0.15);
cfg.soil = struct('scale', 1.0);
cfg.turbid = struct('LAI', 3.0, 'Omega', 1.0, 'hotspot', 0.05);
cfg.numerics = struct('CI_min', 0.01, 'CI_max', 2.50);
cfg.output = struct( ...
    'save_output', true, ...
    'output_file', fullfile(model_dir, 'SIP_optical_output.mat'));
end

function cfg = merge_config(cfg, user_cfg)
if isempty(user_cfg)
    return
end
if ~isstruct(user_cfg)
    error('Configuration input must be a struct.');
end
names = fieldnames(user_cfg);
for i = 1:numel(names)
    name = names{i};
    if isfield(cfg, name) && isstruct(cfg.(name)) && isstruct(user_cfg.(name))
        cfg.(name) = merge_config(cfg.(name), user_cfg.(name));
    else
        cfg.(name) = user_cfg.(name);
    end
end
end

function validate_config(cfg)
scene_type = normalize_scene_type(cfg.scene_type);
require_scalar_range(cfg.geometry.SZA, 0, 89.999, 'SZA');
require_scalar_range(cfg.geometry.SAA, -360, 360, 'SAA');
if strcmp(scene_type, '3D')
    require_scalar_range(cfg.canopy.LAI_crown, eps, Inf, 'LAI_crown');
    require_scalar_range(cfg.canopy.Height, eps, Inf, 'Height');
    require_scalar_range(cfg.canopy.Crowndeepth, eps, cfg.canopy.Height, 'Crowndeepth');
    require_scalar_range(cfg.canopy.lmax_nadir, eps, Inf, 'lmax_nadir');
    require_scalar_range(cfg.canopy.leaf_diameter, eps, Inf, 'leaf_diameter');
    require_scalar_range(cfg.canopy.external_crown_center, 0, cfg.canopy.Height, 'external_crown_center');
else
    require_scalar_range(cfg.turbid.LAI, eps, Inf, 'turbid.LAI');
    require_scalar_range(cfg.turbid.Omega, eps, 5, 'turbid.Omega');
    require_scalar_range(cfg.turbid.hotspot, eps, Inf, 'turbid.hotspot');
end
require_scalar_range(cfg.leaf.N, 1, Inf, 'leaf.N');
require_nonnegative_fields(cfg.leaf, {'Cab', 'Car', 'Ant', 'Brown', 'Cw', 'Cm'});
require_scalar_range(cfg.soil.scale, 0, Inf, 'soil.scale');
require_scalar_range(cfg.numerics.CI_min, 0, Inf, 'CI_min');
require_scalar_range(cfg.numerics.CI_max, cfg.numerics.CI_min + eps, Inf, 'CI_max');
end

function require_nonnegative_fields(s, names)
for i = 1:numel(names)
    require_scalar_range(s.(names{i}), 0, Inf, names{i});
end
end

function gap = load_gap(file_name, var_name)
if ~exist(file_name, 'file')
    error('Required gap file is missing: %s', file_name);
end
s = load(file_name, var_name);
if ~isfield(s, var_name)
    error('File %s does not contain variable %s.', file_name, var_name);
end
gap = double(s.(var_name));
if size(gap, 2) ~= 3
    error('%s must have columns [zenith_deg, azimuth_deg, probability].', var_name);
end
end

function [va, Pi, gap_encoding] = prepare_gap_tables(gap_tot, gap_within, gap_betw)
if ~isequal(gap_tot(:, 1:2), gap_within(:, 1:2))
    error('gap_tot and gap_within use different angular grids.');
end
if ~isequal(gap_tot(:, 1:2), gap_betw(:, 1:2))
    error('gap_tot and gap_betw use different angular grids.');
end
assert_probability(gap_tot(:, 3), 'gap_tot');
assert_probability(gap_within(:, 3), 'gap_within');
assert_probability(gap_betw(:, 3), 'gap_betw');

va = gap_tot(:, 1:2);
gap_sum_error = max(abs(gap_tot(:, 3) - gap_betw(:, 3) - gap_within(:, 3)));
gap_cond_error = max(abs(gap_tot(:, 3) - gap_betw(:, 3) ...
    - (1 - gap_betw(:, 3)) .* gap_within(:, 3)));

if gap_sum_error < 1e-8
    Pi = gap_within(:, 3) ./ max(1 - gap_betw(:, 3), eps);
    gap_encoding = 'nonconditional_additive';
elseif gap_cond_error < 1e-8
    Pi = gap_within(:, 3);
    gap_encoding = 'conditional';
else
    Pi = (gap_tot(:, 3) - gap_betw(:, 3)) ./ max(1 - gap_betw(:, 3), eps);
    gap_encoding = 'derived_from_total';
end
Pi = min(max(Pi, eps), 1);
assert_probability(Pi, 'conditional within-crown gap Pi');
end

function [wave_p, rho_all, tau_all, w_all] = prepare_leaf_optics(leaf)
LRT = prospect_DB(leaf.N, leaf.Cab, leaf.Car, leaf.Ant, leaf.Brown, leaf.Cw, leaf.Cm);
if size(LRT, 2) < 3 || size(LRT, 1) < 2001
    error('PROSPECT-D output must contain at least 2001 rows and 3 columns.');
end
wave_p = LRT(1:2001, 1);
rho_all = LRT(1:2001, 2)';
tau_all = LRT(1:2001, 3)';
w_all = rho_all + tau_all;
if any(~isfinite(w_all)) || any(w_all < -1e-10) || any(w_all > 1 + 1e-10)
    error('Leaf albedo w=rho+tau must be finite and within [0,1].');
end
rho_all = min(max(rho_all, 0), 1);
tau_all = min(max(tau_all, 0), 1);
w_all = min(max(w_all, 0), 1);
end

function rg = prepare_soil_reflectance(soil_file, nwl, scale)
if ~exist(soil_file, 'file')
    error('Required soil spectrum is missing: %s', soil_file);
end
soil_data = load(soil_file);
if size(soil_data, 1) < nwl || size(soil_data, 2) < 2
    error('Soil spectrum must contain at least %d rows and two columns.', nwl);
end
rg = soil_data(1:nwl, 2)' .* scale;
assert_probability(rg, 'soil reflectance');
rg = min(max(rg, 0), 1);
end

function lidf = prepare_lidf(lad)
if lad.TypeLidf == 1
    [lidf, ~] = dladgen(lad.LIDFa, lad.LIDFb);
elseif lad.TypeLidf == 2
    [lidf, ~] = campbell(lad.LIDFa);
else
    error('Unsupported TypeLidf. Use 1 for Verhoef or 2 for Campbell.');
end
lidf = lidf(:);
if any(~isfinite(lidf)) || any(lidf < -1e-12)
    error('Leaf inclination distribution must be finite and non-negative.');
end
end

function [Height_c, Height_c_path, c, HotSpotPar, go_par] = prepare_path_structure(cfg, Pw_s, Pw_h)
Height = cfg.canopy.Height;
Crowndeepth = cfg.canopy.Crowndeepth;
lmax_nadir = cfg.canopy.lmax_nadir;
c = (1 - Pw_h) * Crowndeepth / lmax_nadir;
require_scalar_range(c, 0, 1, 'PATH continuity c');

Height_c_path = Height - lmax_nadir + 0.5 * Crowndeepth;
if cfg.canopy.use_external_crown_center
    Height_c = cfg.canopy.external_crown_center;
else
    Height_c = Height_c_path;
end
require_scalar_range(Height_c, 0, Height, 'Crown center height');

HotSpotPar = cfg.canopy.leaf_diameter / Height;
dthr = estimate_dthr_from_gap(Pw_s, Pw_h, cfg.geometry.SZA);
if ~isfinite(dthr)
    dthr = 1.0;
end
go_par = dthr * Crowndeepth;
require_scalar_range(go_par, eps, Inf, 'GO hotspot width go_par');
end

function terms = allocate_outputs(nang, nwl)
z = zeros(nang, nwl);
terms = struct();
terms.BRF_total = z;
terms.BRF_single = z;
terms.BRF_veg_single = z;
terms.BRF_vegetation = z;
terms.BRF_soil = z;
terms.BRF_veg_soil_interaction = z;
terms.BRF_leaf_C = z;
terms.BRF_leaf_T = z;
terms.BRF_soil_crown = z;
terms.BRF_soil_open = z;
terms.BRF_soil_turbid = z;
terms.BRF_CM = z;
terms.BRF_GCM = z;
terms.K_components = zeros(nang, 4);
terms.CI_terms = zeros(nang, 4);
terms.Pi_terms = zeros(nang, 2);
terms.leaf_terms = zeros(nang, 4);
terms.phase_terms = zeros(nang, 6);
end

function view = compute_view_geometry(t, va, SZA, SAA, gap_tot, gap_betw, Pi_table, ...
    Pw_s, Pi_s, LAI_crown, lidf, Height, Height_c, Crowndeepth, c, HotSpotPar, ...
    go_par, CI_min, CI_max, scene_type, Omega)
tto = va(t, 1);
VAA = va(t, 2);
raa = abs(wrap180(VAA - SAA));

Pw_v = gap_betw(t, 3);
Ptot_v = gap_tot(t, 3);
Pi_v = Pi_table(t);

[Gs, Go, k, K, sob, sof] = PHASE(SZA, tto, raa, lidf);
if K <= 0
    error('Projection term K must be positive.');
end

if strcmp(scene_type, '1D')
    % Homogeneous canopy: no between-crown component; all directional
    % sunlit/soil probabilities are handled by turbid-medium sunshade terms.
    Kc = 1;
    Kt = 0;
    Kg = 0;
    Kz = 0;
    CIs = Omega;
    CIo = Omega;
else
    Kg_raw = Pw_s * Pw_v + get_HSF_go(go_par, SZA, SAA, tto, VAA, Pw_s, Pw_v, Height_c);
    Kg = enforce_bounds(Kg_raw, 0, min(Pw_s, Pw_v), 'Kg');
    Kz = Pw_v - Kg;
    Kct = 1 - Pw_v;

    cos_alpha = cosd(SZA) * cosd(tto) + sind(SZA) * sind(tto) * cosd(VAA - SAA);
    alpha = acosd(min(max(cos_alpha, -1), 1));
    beta = alpha * (1 - sin(pi * c / 2));

    same_azimuth = abs(wrap180(VAA - SAA)) < 1e-10;
    if ((Height - Crowndeepth) < Crowndeepth) && (tto > SZA) && same_azimuth
        Kc = Kct;
    else
        Kc = 0.5 * (1 + cosd(beta)) * Kct;
    end
    Kc = enforce_bounds(Kc, 0, Kct, 'Kc');
    Kt = Kct - Kc;

    % Effective clumping constrained by PATH within-crown gaps.
    % 由 PATH 冠内间隙率约束有效聚集指数。
    CIs_raw = -log(max(Pi_s, eps)) * cosd(SZA) / max(Gs * LAI_crown, eps);
    CIo_raw = -log(max(Pi_v, eps)) * cosd(tto) / max(Go * LAI_crown, eps);
    CIs = min(max(CIs_raw, CI_min), CI_max);
    CIo = min(max(CIo_raw, CI_min), CI_max);
end

Pi_s_model = exp(-Gs * CIs * LAI_crown / max(cosd(SZA), eps));
Pi_v_model = exp(-Go * CIo * LAI_crown / max(cosd(tto), eps));

[kc, kg] = sunshade_H(SZA, tto, raa, Gs, Go, CIs, CIo, LAI_crown, HotSpotPar);
[kc_kt, kg_kt] = sunshade_Kt_He(SZA, tto, raa, Gs, Go, CIs, CIo, LAI_crown);

view = struct();
view.Ptot_v = Ptot_v;
view.Pi_v = Pi_v;
view.Kc = Kc;
view.Kt = Kt;
view.Kg = Kg;
view.Kz = Kz;
view.Gs = Gs;
view.Go = Go;
view.k = k;
view.K = K;
view.sob = sob;
view.sof = sof;
view.CIs = CIs;
view.CIo = CIo;
view.Pi_s = Pi_s;
view.Pi_s_model = Pi_s_model;
view.Pi_v_model = Pi_v_model;
view.kc = kc;
view.kc_kt = kc_kt;
view.kg = kg;
view.kg_kt = kg_kt;
end

function [terms, row] = solve_optical_row(terms, t, view, rho, tau, w, rg, Ptot_s, iD, LAI_scene, p, scene_type)
wso = view.sob * rho + view.sof * tau;
BRF_leaf_C = view.Kc .* (wso .* view.kc ./ view.K);
BRF_leaf_T = view.Kt .* (sqrt(view.Pi_s) .* wso .* view.kc_kt ./ view.K);
BRF_soil_canopy_path = view.Kc .* (view.kg .* rg) + view.Kt .* (view.kg_kt .* rg);
BRF_soil_open_path = (view.Kg + view.Kz .* view.Pi_s) .* rg;

if strcmp(scene_type, '1D')
    BRF_soil_crown = zeros(size(rg));
    BRF_soil_open = zeros(size(rg));
    BRF_soil_turbid = BRF_soil_canopy_path + BRF_soil_open_path;
else
    BRF_soil_crown = BRF_soil_canopy_path;
    BRF_soil_open = BRF_soil_open_path;
    BRF_soil_turbid = zeros(size(rg));
end

BRF_veg_single = BRF_leaf_C + BRF_leaf_T;
BRF_soil = BRF_soil_crown + BRF_soil_open + BRF_soil_turbid;
BRF_single = BRF_veg_single + BRF_soil;

i0 = 1 - Ptot_s;
iv = 1 - view.Ptot_v;
t0 = 1 - i0;
tv = 1 - iv;
rho_o = iv / (2 * LAI_scene);
rho_hemi = iD / (2 * LAI_scene);

denom = 1 - p .* w;
if any(denom <= 0)
    error('Energy denominator 1 - p*w must be positive.');
end

Tdn = t0 + i0 .* w .* rho_hemi ./ denom;
Tup = tv + iD .* w .* rho_o ./ denom;
Rdn = iD .* w .* rho_hemi ./ denom;
ground_denom = 1 - rg .* Rdn;
if any(ground_denom <= 0)
    error('Energy denominator 1 - soil_reflectance*Rdn must be positive.');
end

BRF_CM = i0 .* (w .^ 2) .* p .* rho_o ./ denom;
BRF_GCM = rg .* Tdn .* Tup ./ ground_denom - t0 .* rg .* tv;

BRF_vegetation = BRF_veg_single + BRF_CM;
BRF_interaction = BRF_GCM;
BRF_total = BRF_vegetation + BRF_soil + BRF_interaction;

if any(~isfinite(BRF_total))
    error('Non-finite BRF values produced at angular row %d.', t);
end

terms.BRF_leaf_C(t, :) = BRF_leaf_C;
terms.BRF_leaf_T(t, :) = BRF_leaf_T;
terms.BRF_soil_crown(t, :) = BRF_soil_crown;
terms.BRF_soil_open(t, :) = BRF_soil_open;
terms.BRF_soil_turbid(t, :) = BRF_soil_turbid;
terms.BRF_CM(t, :) = BRF_CM;
terms.BRF_GCM(t, :) = BRF_GCM;
terms.BRF_veg_single(t, :) = BRF_veg_single;
terms.BRF_soil(t, :) = BRF_soil;
terms.BRF_single(t, :) = BRF_single;
terms.BRF_vegetation(t, :) = BRF_vegetation;
terms.BRF_veg_soil_interaction(t, :) = BRF_interaction;
terms.BRF_total(t, :) = BRF_total;

row = struct();
row.K_components = [view.Kc, view.Kt, view.Kg, view.Kz];
row.CI_terms = [view.CIs, view.CIo, view.Pi_s_model, view.Pi_v_model];
row.Pi_terms = [view.Pi_s, view.Pi_v];
row.leaf_terms = [view.kc, view.kc_kt, view.kg, view.kg_kt];
row.phase_terms = [view.Gs, view.Go, view.k, view.K, view.sob, view.sof];
end

function idx = find_gap_row(gap_table, zenith, azimuth)
tol = 1e-6;
idx = find(abs(gap_table(:, 1) - zenith) < tol ...
    & abs(wrap180(gap_table(:, 2) - azimuth)) < tol, 1);
if isempty(idx)
    error('No gap row found for zenith %.3f deg and azimuth %.3f deg.', zenith, azimuth);
end
end

function scene_type = normalize_scene_type(scene_type)
if isstring(scene_type)
    scene_type = char(scene_type);
end
if ~ischar(scene_type)
    error('scene_type must be ''3D'' or ''1D''.');
end
switch lower(strtrim(scene_type))
    case {'3d', 'discrete', 'discrete3d'}
        scene_type = '3D';
    case {'1d', 'homogeneous', 'turbid', 'turbid_medium'}
        scene_type = '1D';
    otherwise
        error('scene_type must be ''3D'' or ''1D''.');
end
end

function va = prepare_view_angles(user_va, fallback_va)
if isempty(user_va)
    va = fallback_va;
    return
end
if ~isnumeric(user_va) || size(user_va, 2) ~= 2 || isempty(user_va)
    error('geometry.view_angles must be an N-by-2 matrix: [VZA_deg, VAA_deg].');
end
if any(~isfinite(user_va(:))) || any(user_va(:, 1) < 0) || any(user_va(:, 1) >= 90)
    error('View zenith angles must be finite and within [0, 90) degrees.');
end
va = user_va;
end

function [gap_tot, gap_betw, Pi_table, Ptot_s, Pi_s, iD, T_gauss, T_trapz] = ...
    prepare_1d_gap_terms(va, SZA, SAA, lidf, LAI, Omega)
P_v = zeros(size(va, 1), 1);
for i = 1:size(va, 1)
    P_v(i) = directional_gap_1d(va(i, 1), lidf, LAI, Omega);
end
Ptot_s = directional_gap_1d(SZA, lidf, LAI, Omega);
Pi_s = Ptot_s;
gap_tot = [va, P_v];
gap_betw = [va, zeros(size(P_v))];
Pi_table = min(max(P_v, eps), 1);
T_gauss = hemispherical_gap_1d(lidf, LAI, Omega, 'gauss');
T_trapz = hemispherical_gap_1d(lidf, LAI, Omega, 'trapz');
iD = 1 - T_gauss;
assert_probability(gap_tot(:, 3), '1D total gap');
assert_probability(Pi_table, '1D directional transmittance');
end

function P = directional_gap_1d(theta, lidf, LAI, Omega)
if theta >= 89.999
    P = 0;
    return
end
[Gtheta, ~, ~, ~, ~, ~] = PHASE(theta, 0, 0, lidf);
P = exp(-Gtheta * Omega * LAI / max(cosd(theta), eps));
P = min(max(P, 0), 1);
end

function T = hemispherical_gap_1d(lidf, LAI, Omega, method)
if strcmp(method, 'gauss')
    x = [0; 0.5384693101; -0.5384693101; 0.9061798459; -0.9061798459];
    w = [0.5688888889; 0.4786286705; 0.4786286705; 0.2369268850; 0.2369268850];
    theta = (pi / 4) .* x + pi / 4;
    weight = (pi / 4) .* w;
    P = arrayfun(@(th) directional_gap_1d(rad2deg(th), lidf, LAI, Omega), theta);
    T = sum(weight .* 2 .* P .* sin(theta) .* cos(theta));
else
    theta = linspace(0, pi / 2, 181)';
    P = arrayfun(@(th) directional_gap_1d(rad2deg(th), lidf, LAI, Omega), theta);
    T = trapz(theta, 2 .* P .* sin(theta) .* cos(theta));
end
T = min(max(T, 0), 1);
end

function a = wrap180(a)
a = mod(a + 180, 360) - 180;
end

function dthr = estimate_dthr_from_gap(Pw_s, Pw_h, SZA)
if SZA == 0 || Pw_s <= 0 || Pw_h <= 0 || Pw_s >= 1 || Pw_h >= 1
    dthr = NaN;
    return
end
ratio = log(Pw_s) / log(Pw_h);
term = (ratio ^ 2 - 1) / (tand(SZA) ^ 2);
if term <= 0
    dthr = NaN;
else
    dthr = 1 / sqrt(term);
end
end

function require_scalar_range(x, lo, hi, name)
if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x) || x < lo || x > hi
    error('%s must be a finite scalar within [%g, %g].', name, lo, hi);
end
end

function assert_probability(x, name)
tol = 1e-10;
if any(~isfinite(x(:))) || any(x(:) < -tol) || any(x(:) > 1 + tol)
    error('%s must be finite and within [0,1].', name);
end
end

function x = enforce_bounds(x, lo, hi, name)
tol = 1e-8;
if x < lo - tol || x > hi + tol
    error('%s=%.12g is outside physical bounds [%.12g, %.12g].', ...
        name, x, lo, hi);
end
x = min(max(x, lo), hi);
end
