# DeepExcavation2D

2D analytical assessment of ground movements and building damage due to deep excavation.

## Method

Wall deflection profiles (parametric shapes M0–M5) are converted to equivalent
circular cavities superposed using the analytical solution of Zheng et al. (2023).
Building damage is assessed using the Timoshenko equivalent beam model (SSI) from
ASREp, with a greenfield reference computed via the Boscardin & Cording (1989)
limiting tensile strain method.

## Files

| File | Description |
|---|---|
| `Input.py` | All input parameters (wall, soil, building, grid) |
| `main.py` | Runs greenfield model + SSI Timoshenko, prints summary |
| `results_viz.ipynb` | Visualisation: contour plots, profiles, strain results |

## How to run

```bash
cd Examples/DeepExcavation2D
python main.py


References
Ground movement model:

Zheng, C., Franza, A., & Jimenez, R. (2023). Analytical prediction for ground
movements due to deep excavations in soils. Tunnelling and Underground Space
Technology, 141, 105316.
https://doi.org/10.1016/j.tust.2023.105316

Damage assessment:

Boscardin, M. D., & Cording, E. J. (1989). Building response to excavation-induced
settlement. Journal of Geotechnical Engineering, 115(1), 1–21.