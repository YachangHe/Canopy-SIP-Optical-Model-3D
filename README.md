# SIP Optical Canopy Reflectance Model

SIP optical canopy reflectance model combining geometric-optical and spectral invariants theories.

## Model Overview

This package provides the optical reflectance branch of the SIP canopy model. It combines:

- Geometric-optical (GO) canopy structure and sun-view geometry.
- Direction-dependent hotspot and mutual-shadowing terms for sun-view geometry.
- Leaf optical properties from PROSPECT-D.
- Spectral-invariant multiple scattering and vegetation-soil interaction terms.

For the 3D discrete scene, the four canopy components are derived from gap probabilities using PATH model. The radiative transfer closure follows the SIP framework that combines geometric-optical and spectral invariants theories.

The user-facing workflow is intentionally small: edit the parameter block in `main_SIP_optical.m`, run the script, and read `outputs` and `diagnostics`.

## Authorship

- The first SIP model for 1D homogeneous canopies was developed by Yelu Zeng, Min Chen, Dalei Hao, and collaborators.
- The 3D SIP model extension combining geometric-optical and spectral invariants theories was developed primarily by Yachang He, Yelu Zeng, and Dalei Hao.

## Scene Types

Set `cfg.scene_type` in `main_SIP_optical.m`:

```matlab
cfg.scene_type = '3D';  % GO + spectral-invariants discrete crown scene
cfg.scene_type = '1D';  % homogeneous/turbid-medium scene
```

In `1D` mode, directional gap probability is computed as:

```matlab
P_gap(theta) = exp(-G(theta) * Omega * LAI / cos(theta))
```

where `cfg.turbid.LAI` is the homogeneous-scene LAI and `cfg.turbid.Omega` is an optional macroscopic clumping factor. The 3D-only crown geometry is ignored in `1D` mode.

The following parameters are not required in `1D` mode:

- `cfg.canopy.Height`
- `cfg.canopy.Crowndeepth`
- `cfg.canopy.lmax_nadir`
- `cfg.canopy.external_crown_center`
- Crown-center and between-crown gap terms used by the 3D component partition

Use `cfg.turbid.hotspot` instead of the 3D crown-scale hotspot width when running the 1D model.

## Coordinate Convention

All angular inputs are in degrees.

| Symbol | Meaning | Unit | Convention |
|---|---|---:|---|
| `SZA` | Solar zenith angle | degree | `0` is vertical/zenith |
| `VZA` | View zenith angle | degree | `0` is nadir |
| `SAA` | Solar azimuth angle | degree | horizontal-plane azimuth |
| `VAA` | View azimuth angle | degree | read from the angular gap-probability table |
| `RAA` | Relative azimuth angle | degree | `abs(wrapTo180(VAA - SAA))` |

## Prerequisites & Installation

- MATLAB R2018b or newer is recommended.
- No MATLAB toolbox is required for the default demo.
- Keep the directory structure unchanged because `data/CI_HET10` and `data/soilnew_1_to_LESS.txt` are loaded by relative path.

Run from MATLAB:

```matlab
cd SIP_Optical_Model
main_SIP_optical
```

## Quick Start

Edit the top block of `main_SIP_optical.m`:

```matlab
cfg.geometry.SZA = 20;
cfg.geometry.SAA = 0;
cfg.scene_type = '3D';      % or '1D'
cfg.turbid.LAI = 3.0;       % used only in 1D mode
cfg.turbid.Omega = 1.0;     % used only in 1D mode
cfg.canopy.LAI_crown = 5.0; % used only in 3D mode
cfg.leaf.Cab = 30.0;
cfg.soil.scale = 1.0;

[outputs, diagnostics] = sip_optical_core(cfg);
```

The demo saves `SIP_optical_output.mat` by default. To keep results only in the MATLAB workspace:

```matlab
cfg.output.save_output = false;
```

## Inputs

| Field | Unit | Description |
|---|---:|---|
| `cfg.geometry.SZA` | degree | Solar zenith angle |
| `cfg.geometry.SAA` | degree | Solar azimuth angle |
| `cfg.geometry.view_angles` | degree | Optional `[VZA, VAA]` grid for 1D mode |
| `cfg.scene_type` | - | `'3D'` discrete crown scene or `'1D'` homogeneous scene |
| `cfg.turbid.LAI` | m2 m-2 | 1D homogeneous-scene LAI |
| `cfg.turbid.Omega` | - | 1D clumping factor |
| `cfg.turbid.hotspot` | - | 1D hotspot size parameter |
| `cfg.canopy.LAI_crown` | m2 m-2 | 3D manually specified single-crown LAI |
| `cfg.canopy.Height` | m | Canopy top height |
| `cfg.canopy.Crowndeepth` | m | Mean path depth after a crown hit |
| `cfg.canopy.leaf_diameter` | m | Leaf-size parameter controlling hotspot width |
| `cfg.leaf.N` | - | PROSPECT-D leaf structure parameter |
| `cfg.leaf.Cab` | ug cm-2 | Chlorophyll a+b content |
| `cfg.leaf.Car` | ug cm-2 | Carotenoid content |
| `cfg.leaf.Ant` | ug cm-2 | Anthocyanin content |
| `cfg.leaf.Cw` | cm | Equivalent water thickness |
| `cfg.leaf.Cm` | g cm-2 | Dry matter content |
| `cfg.lad.TypeLidf` | - | `1` Verhoef, `2` Campbell |
| `cfg.soil.scale` | - | Soil reflectance multiplier |

## Outputs

| Output | Description |
|---|---|
| `outputs.wavelength` | Wavelength grid in nm |
| `outputs.view_angles` | `[VZA, VAA]` angular grid from the gap-probability table |
| `outputs.BRF_total` | Total canopy bidirectional reflectance factor |
| `outputs.BRF_single` | Single-scattering BRF from leaves plus direct soil contribution |
| `outputs.BRF_veg_single` | Leaf single-scattering BRF before canopy multiple scattering |
| `outputs.BRF_vegetation` | Vegetation contribution, including leaf single scattering and canopy multiple scattering |
| `outputs.BRF_soil` | Direct soil contribution. In 3D it is split into crown/open components; in 1D it is stored as the turbid-medium soil term |
| `outputs.BRF_veg_soil_interaction` | Vegetation-soil multiple-interaction contribution |
| `outputs.components.BRF_leaf_C` | 3D crown-facing leaf single-scattering term; in 1D this stores the homogeneous-canopy leaf single-scattering term |
| `outputs.components.BRF_leaf_T` | 3D mutual-shadowed leaf single-scattering term; zero in 1D mode |
| `outputs.components.BRF_soil_crown` | 3D soil contribution viewed through crown-covered gaps; zero in 1D mode |
| `outputs.components.BRF_soil_open` | 3D soil contribution from open/between-crown gaps; zero in 1D mode |
| `outputs.components.BRF_soil_turbid` | 1D homogeneous-canopy direct soil contribution; zero in 3D mode |
| `outputs.components.BRF_CM` | Canopy multiple-scattering contribution from spectral-invariant recurrence |
| `outputs.components.BRF_GCM` | Vegetation-soil coupled multiple-interaction contribution |
| `diagnostics` | Configuration, gap probabilities, scene LAI, p-theory parameters, clumping terms, phase terms, and component probabilities |

## Physical Guardrails

The core function checks:

- LAI, canopy height, crown depth, hotspot size, and spectral parameters are physically valid.
- Gap probabilities and soil reflectance remain within `[0, 1]`.
- Leaf single-scattering albedo satisfies `0 <= rho + tau <= 1`.
- Energy denominators `1 - p*w` and `1 - soil_reflectance*Rdn` stay positive.
- 3D component probabilities (`Kc`, `Kt`, `Kg`, `Kz`) remain inside their physical bounds.

## Citation

If you use the 1D SIP mode, please cite:

```bibtex
@article{Zeng2018SpectralInvariant,
  title   = {Spectral invariant provides a practical modeling approach for future biophysical variable estimations},
  author  = {Zeng, Yelu and Xu, Baodong and Yin, Gaofei and Wu, Shengbiao and Hu, Guoqing and Yan, Kai and Yang, Bin and Song, Wanjuan and Li, Jing},
  journal = {Remote Sensing},
  year    = {2018},
  volume  = {10},
  number  = {10},
  pages   = {1508}
}
```

If you use the 3D SIP mode or the GO + spectral-invariants framework, please cite:

```bibtex
@article{He2025GOSpectralInvariantsSIF,
  title   = {Combining geometric-optical and spectral invariants theories for modeling canopy fluorescence anisotropy},
  author  = {He, Yachang and Zeng, Yelu and Hao, Dalei and Shabanov, Nikolay V. and Huang, Jianxi and Yin, Gaofei and Biriukova, Khelvi and others},
  journal = {Remote Sensing of Environment},
  year    = {2025},
  volume  = {323},
  pages   = {114716},
  doi     = {10.1016/j.rse.2025.114716}
}
```

## Acknowledgements

The authors thank Dr. Weihua Li for providing the PATH_RT model, and Dr. Jianbo Qi and Bang Sun for the guidance on the LESS model.

## License

Add the final open-source license before public release.
