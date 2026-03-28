# Mitochondrial Transplantation Metabolic Calculator

An interactive, web-based metabolic modeling tool for evaluating whether mitochondrial transplantation (MT) can overcome the ~200 μm oxygen diffusion limit in engineered tissue constructs.

## Live Demo

Open `index.html` in a browser, or serve locally:

```bash
python3 -m http.server 8080
# Then navigate to http://localhost:8080
```

## What It Models

### 01 — Krogh Oxygen Diffusion Model
Models oxygen penetration depth as a function of diffusion coefficient, surface O₂ concentration, cell density, and MT efficiency boost using the classical Krogh equation:

**L = √(2 · D · C₀ / Q)** where Q = ρ · q / η

### 02 — ATP Budget Analysis
Compares OXPHOS vs glycolytic ATP supply against cellular demand at any given depth and oxygen availability.

### 03 — Buying Time: Survival vs Vascularization
Overlays cell death curves against vascularization timelines to determine if MT can extend survival long enough for vascular networks to establish. Supports three vascularization strategies:
- Pre-vascularized (inosculation ~2–5 days)
- Standard construct (neovascularization ~10–15 days)
- Growth factor-enhanced (~7–10 days)

### 04 — Depth-Resolved Viability Timeline
Simulates cell viability over 30 days at different construct depths, modeling MT benefit decay (exponential half-life) and optional re-dosing strategies.

### 05 — Mathematical Verdict
Automated synthesis of all models into actionable recommendations.

## Key Parameters (Literature-Derived)

| Parameter | Default | Source |
|-----------|---------|--------|
| O₂ diffusion coefficient | 2.0 × 10⁻⁵ cm²/s | Hydrogel/collagen scaffolds |
| Surface O₂ (C₀) | 0.21 mol/m³ | Normoxia at 37°C |
| Cell density | 6.0 × 10⁸ cells/mL | FRESH bioprint specs |
| OCR per cell | 25 amol O₂/cell/s | hiPSC-CM range |
| MT efficiency boost | 1.0–2.0× | McCully/Emani Seahorse XF data |
| MT half-life | 5 days | Peak at ~2 days, exponential decay |

## Cell Type Presets

| Cell Type | OCR (amol O₂/cell/s) |
|-----------|----------------------|
| Cardiomyocyte | 20 |
| Hepatocyte | 40 |
| Renal tubular | 15 |
| MSC | 5 |

## Key Mathematical Insight

Penetration depth scales as **√η** (square root of efficiency boost). A 2× efficiency improvement yields only ~1.4× depth improvement. This fundamental constraint means MT is best positioned as a metabolic resilience layer within a multi-pronged approach (MT + pre-vascularization + oxygen-releasing biomaterials).

## Literature References

- Malda et al., 2008 — O₂ depletion onset in static 3D culture
- Das et al., 2009; Potier et al., 2007 — MSC hypoxia tolerance
- Nazeer et al., 2021 — Neovascularization timelines
- Tremblay et al., 2005 — Inosculation kinetics
- Ramirez-Barbieri et al., 2019; Patel et al., 2023 — MT efficiency data
- McCully/Emani groups — Seahorse XF validation of MT OXPHOS improvement

## Tech Stack

- Vanilla HTML/CSS/JavaScript
- [Chart.js](https://www.chartjs.org/) v4.4.1 (via CDN)
- No build step required
