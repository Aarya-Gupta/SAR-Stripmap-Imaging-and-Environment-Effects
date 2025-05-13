# SAR Stripmap Imaging and Environmental Effects

## Project Overview

Under the umbrella of **Synthetic Aperture Radar (SAR) Stripmap Imaging**, this project guides you through setting up a sidelooking SAR problem in MATLAB, generating stripmap images of point targets and realistic terrain, and exploring the impact of environmental conditions and system parameters on image quality.

## Project Tasks

1. **SAR Model Setup:**

   * Use MATLAB’s Radar Toolbox to configure a sidelooking SAR with squint angle = 0°.
   * Platform wavelength (λ), velocity (v), and altitude (h) are specified in your code.
2. **Point-Target Imaging:**

   * Transmit 100% duty-cycle chirp pulses with specified chirp rate (\$\alpha\$) and PRI (T$\_	ext{PRI}\$).
   * Define transmit power (P$\_	ext{TX}\$) and antenna gain (G), and assume free-space propagation.
   * Simulate a single ground point target and generate a stripmap image using back-projection or range-Doppler processing.
3. **Terrain Stripmap:**

   * Import a Digital Elevation Model (DEM) in MATLAB.
   * Repeat imaging to produce a terrain stripmap under the same system settings.
4. **Environmental Effects:**

   * Model propagation losses for **rain, fog, and snowfall** at your chosen radar frequency.
   * Generate stripmap images for the point target and terrain under each condition.
5. **Parameter Sweeps:**

   * Vary one parameter at a time—altitude, elevation angle, gain, λ, v, PRI, and P$\_	ext{TX}\$.
   * Document the effect on image focus, SNR, and ambiguity.
6. **Reporting:**

   * **Team Contributions:** Clearly list each member’s responsibilities.
   * **Code:** Submit all “.m” scripts with inline comments.
   * **Results:** Provide plots with labeled axes for every scenario.
   * **Analysis:** Summarize observations, compare cases, and draw insights on parameter impacts.

## Repository Structure

```
SAR_Stripmap_Project/
│   README.md
│
├───DEM/                     # Digital Elevation Model files
│       ground_dem.mat
│       ground_dem.tif
│
├───scripts/                 # MATLAB code
│   ├── sar_setup.m          # Define platform, waveform, and system parameters
│   ├── point_target_imaging.m
│   ├── terrain_imaging.m
│   ├── environmental_effects.m
│   └── parameter_sweep.m
│
├───results/                 # Generated images and data
│   ├── point_clear.png
│   ├── point_rain.png
│   ├── terrain_fog.png
│   ├── sweep_gain.png
│   └── ...
│
└───reports/                 # Final write–up and supplementary materials
    ├── SAR_Project_Report.pdf
    └── team_contributions.txt
```

## Getting Started

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/your-org/SAR_Stripmap_Project.git
   ```
2. **Prerequisites:**

   * MATLAB R2024b or later with Radar Toolbox.
   * Mapping Toolbox (for DEM import) optional but recommended.
3. **Run SAR Setup:**

   * Open `scripts/sar_setup.m`, customize λ, v, h, α, T$\_	ext{PRI}\$, P$\_	ext{TX}\$, and G.
   * Execute to generate basic SAR parameters and plots.
4. **Image Generation:**

   * Run `point_target_imaging.m` and `terrain_imaging.m` for baseline images.
   * Execute `environmental_effects.m` to see weather impacts.
   * Use `parameter_sweep.m` to explore system parameter variations.

## Submission Guidelines

* **Deadline:** May 5, 2025, 12:00 AM (midnight).
* **Deliverables:**

  * Final report in `reports/SAR_Project_Report.pdf`.
  * All MATLAB scripts in `scripts/` with proper comments.
  * Team contributions file `reports/team_contributions.txt`.
* **Evaluation:**

  * Clarity and correctness of code.
  * Quality of images (focus, SNR, ambiguity).
  * Depth of analysis and insights.

---

Good luck, and happy imaging!