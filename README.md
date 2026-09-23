# ANSG (Applied Nuclear Science Group) Analysis
<!---->
## Projects
<!---->
- **73mGe** — Precision gamma-ray measurements of 73mGe using a state-of-the-art position sensitive CZT detector.
- **78mBr** — Half-life measurement of 78mBr by self-coincidence analysis in a LaBr3 detector.  
  V. Fondement, E. Todd, I. Jovanovic, “Updated ⁷⁸ᵐBr 2⁻ → 1⁺ Lifetime from Isomeric Decay Cascade in Active Target,” *Nuclear Physics A* (2026), [doi:10.1016/j.nuclphysa.2026.123485](https://doi.org/10.1016/j.nuclphysa.2026.123485).
- **YAG-PreProcessing** — Pre-processing of YAG detector data.
- **YAP-PSD** — Machine learning-based pulse shape discrimination in YAP:Ce, comparing Random Forest, Gradient Boosting, XGBoost, and MLP models.  
  E. Todd, V. Fondement, I. Jovanovic, “Machine Learning-Enabled Pulse Shape Discrimination in YAP:Ce,” *Nuclear Instruments and Methods in Physics Research A* (2026), [doi:10.1016/j.nima.2026.171883](https://doi.org/10.1016/j.nima.2026.171883).
<!---->
## Dependencies
<!---->
All projects depend on [Analysis-Utilities](https://github.com/ewtodd/Analysis-Utilities), a shared library providing common ROOT-based analysis tools.
<!---->
## Reproducibility
<!---->
This project uses [Nix](https://nixos.org/) to manage dependencies.
Install it by following the instructions [here](https://nixos.org/download/).
Ensure [flakes are enabled](https://nixos.wiki/wiki/Flakes#Enable_flakes_permanently) in your Nix configuration.
The repository is one flake with a development shell per project:
<!---->
```
nix develop .#73mGe
nix develop .#78mBr
nix develop .#YAG-PreProcessing
nix develop .#YAP-PSD
```
<!---->
Every shell pins the same Analysis-Utilities and the ROOT it was built against, and starts in the project directory (the `macros/` directory for 73mGe and 78mBr). The 73mGe shell decrypts its agenix secrets on entry and needs the matching ssh key.
Contact me for access to data.
