/* ═══════════════════════════════════════════════════════════
   MITOCHONDRIAL TRANSPLANT METABOLIC CALCULATOR
   Core Computation Engine & Chart Rendering
   
   Integrates models from:
   - metabolic_demand_calculator.html (Krogh + cell presets)
   - buying_time_survival_model.html (survival vs vascularization)
   ═══════════════════════════════════════════════════════════ */

// ── Chart.js Global Config ────────────────────────────────
Chart.defaults.color = '#606060';
Chart.defaults.font.family = "'JetBrains Mono', monospace";
Chart.defaults.font.size = 10;
Chart.defaults.plugins.legend.labels.padding = 12;
Chart.defaults.plugins.legend.labels.usePointStyle = true;
Chart.defaults.plugins.legend.labels.pointStyleWidth = 10;

// ── Cell Type Presets (matching metabolic_demand_calculator.html) ──
const cellTypeData = {
  hepatocyte:     { ocr: 40 },
  cardiomyocyte:  { ocr: 20 },
  renal:          { ocr: 15 },
  msc:            { ocr: 5 },
  custom:         { ocr: 25 }   // default custom value
};

// ── Vascularization Data (matching buying_time_survival_model.html) ──
const vascData = {
  prevascularized:  { onset: 1, complete: 4,  label: 'Inosculation' },
  standard:         { onset: 5, complete: 15, label: 'Neovascularization' },
  growth_factors:   { onset: 3, complete: 10, label: 'GF-enhanced angiogenesis' }
};

// ── State ─────────────────────────────────────────────────
let kroghChart = null;
let atpChart = null;
let buyingTimeChart = null;
let viabilityChart = null;

// ── DOM References ────────────────────────────────────────
const sliders = {
  cellType:        document.getElementById('cell-type'),
  diffCoeff:       document.getElementById('diffusion-coeff'),
  surfaceO2:       document.getElementById('surface-o2'),
  cellDensity:     document.getElementById('cell-density'),
  ocrPerCell:      document.getElementById('ocr-per-cell'),
  mtEfficiency:    document.getElementById('mt-efficiency'),
  atpDemand:       document.getElementById('atp-demand'),
  oxphosEfficiency:document.getElementById('oxphos-efficiency'),
  glycolysisRate:  document.getElementById('glycolysis-rate'),
  localO2Pct:      document.getElementById('local-o2-pct'),
  vascStrategy:    document.getElementById('vasc-strategy'),
  hypoxiaTolerance:document.getElementById('hypoxia-tolerance'),
  o2DepletionDays: document.getElementById('o2-depletion-days'),
  constructThickness: document.getElementById('construct-thickness'),
  viabilityThreshold: document.getElementById('viability-threshold'),
  mtDecayHalflife: document.getElementById('mt-decay-halflife'),
  reDoseInterval:  document.getElementById('re-dose-interval'),
};

const displays = {
  diffCoeff:       document.getElementById('diffusion-coeff-val'),
  surfaceO2:       document.getElementById('surface-o2-val'),
  cellDensity:     document.getElementById('cell-density-val'),
  ocrPerCell:      document.getElementById('ocr-per-cell-val'),
  mtEfficiency:    document.getElementById('mt-efficiency-val'),
  atpDemand:       document.getElementById('atp-demand-val'),
  oxphosEfficiency:document.getElementById('oxphos-efficiency-val'),
  glycolysisRate:  document.getElementById('glycolysis-rate-val'),
  localO2Pct:      document.getElementById('local-o2-pct-val'),
  hypoxiaTolerance:document.getElementById('hypoxia-tolerance-val'),
  o2DepletionDays: document.getElementById('o2-depletion-days-val'),
  constructThickness: document.getElementById('construct-thickness-val'),
  viabilityThreshold: document.getElementById('viability-threshold-val'),
  mtDecayHalflife: document.getElementById('mt-decay-halflife-val'),
  reDoseInterval:  document.getElementById('re-dose-interval-val'),
};

// ── Cell Type Dropdown Handler ────────────────────────────
function handleCellTypeChange() {
  const type = sliders.cellType.value;
  if (type !== 'custom') {
    const ocr = cellTypeData[type].ocr;
    sliders.ocrPerCell.value = ocr;
    displays.ocrPerCell.textContent = ocr.toFixed(0);
  }
  updateAll();
}

// ── Physics Functions ─────────────────────────────────────

/**
 * Krogh penetration depth (1D planar slab)
 * L = sqrt(2 * D * C0 / Q)
 * 
 * Matches metabolic_demand_calculator.html:
 *   D = 2e-9 m²/s, C0 = 0.2 mol/m³
 *   Q = OCR * density (volumetric consumption)
 * 
 * Our version is generalized to accept adjustable D and C0.
 */
function kroghPenetration(D, C0, rho, q, eta) {
  // D: slider value × 10⁻⁵ cm²/s → m²/s
  const D_m2s = D * 1e-9;

  // C0: mol/m³ (direct)

  // rho: × 10⁸ cells/mL → cells/m³
  const rho_m3 = rho * 1e14;

  // q: amol/cell/s = 10⁻¹⁸ mol/cell/s
  const q_mol = q * 1e-18;

  // Q (mol/m³/s) = rho * q / eta
  const Q = rho_m3 * q_mol / eta;

  if (Q <= 0) return Infinity;

  const L_m = Math.sqrt(2 * D_m2s * C0 / Q);
  return L_m * 1e6; // μm
}

/**
 * Oxygen concentration profile C(x)
 */
function oxygenProfile(x_um, L_um, C0) {
  if (x_um >= L_um) return 0;
  return C0 * (1 - (x_um * x_um) / (L_um * L_um));
}

/**
 * ATP supply rate at a given depth
 * 
 * MT effect: improves OXPHOS efficiency (more ATP per O₂ consumed).
 * Validated by McCully/Emani groups showing 1.2–2× basal respiration
 * improvement via Seahorse XF assays.
 * 
 * Note: At depths where O₂ = 0, only glycolysis contributes.
 * MT does not affect glycolytic ATP production in this model —
 * this is consistent with current literature where MT benefits
 * are primarily mitochondrial (ETC function, membrane potential).
 */
function atpSupplyAtDepth(localO2Fraction, mtEfficiency, poRatio, glycolyticFraction, demand) {
  const oxphosBaseline = (1 - glycolyticFraction) * demand;
  const oxphosActual = oxphosBaseline * localO2Fraction * (poRatio / 2.5) * mtEfficiency;
  const glycolyticActual = glycolyticFraction * demand;
  return oxphosActual + glycolyticActual;
}

/**
 * MT efficiency over time with decay and optional re-dosing
 */
function mtEfficiencyOverTime(t_days, eta0, halflife, reDoseInterval) {
  let totalBoost = 0;
  const boost0 = (eta0 - 1) * Math.exp(-0.693 * t_days / halflife);
  totalBoost += boost0;

  if (reDoseInterval > 0) {
    let doseTime = reDoseInterval;
    while (doseTime <= t_days) {
      const timeSinceDose = t_days - doseTime;
      const boostN = (eta0 - 1) * Math.exp(-0.693 * timeSinceDose / halflife);
      totalBoost += boostN;
      doseTime += reDoseInterval;
    }
    totalBoost = Math.min(totalBoost, (eta0 - 1) * 3);
  }

  return 1 + totalBoost;
}

/**
 * Sigmoid function (matching buying_time_survival_model.html exactly)
 */
function sigmoid(x, midpoint, steepness) {
  return 1 / (1 + Math.exp(steepness * (x - midpoint)));
}

/**
 * Vascularization progress curve (matching buying_time_survival_model.html)
 * Linear ramp from onset to complete
 */
function vascCurve(day, onset, complete) {
  if (day < onset) return 0;
  if (day >= complete) return 100;
  return 100 * (day - onset) / (complete - onset);
}

// ── Chart Creation ────────────────────────────────────────

function createKroghChart() {
  const ctx = document.getElementById('krogh-chart').getContext('2d');
  kroghChart = new Chart(ctx, {
    type: 'line',
    data: {
      labels: [],
      datasets: [
        {
          label: 'Baseline (no MT)',
          data: [],
          borderColor: 'rgba(255, 255, 255, 0.25)',
          backgroundColor: 'rgba(255, 255, 255, 0.02)',
          fill: true,
          tension: 0.3,
          pointRadius: 0,
          borderWidth: 1.5,
          borderDash: [4, 2],
        },
        {
          label: 'With Mito-Transplant',
          data: [],
          borderColor: 'rgba(255, 255, 255, 0.85)',
          backgroundColor: 'rgba(255, 255, 255, 0.05)',
          fill: true,
          tension: 0.3,
          pointRadius: 0,
          borderWidth: 2,
        },
        {
          label: 'Hypoxia Threshold',
          data: [],
          borderColor: 'rgba(255, 255, 255, 0.12)',
          borderDash: [6, 4],
          pointRadius: 0,
          borderWidth: 1,
          fill: false,
        }
      ]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      interaction: { mode: 'index', intersect: false },
      scales: {
        x: {
          title: { display: true, text: 'Depth from surface (μm)', color: '#4e4e5e', font: { size: 10, family: "'JetBrains Mono', monospace" } },
          grid: { color: 'rgba(255, 255, 255, 0.04)' },
          ticks: { maxTicksLimit: 10 },
        },
        y: {
          title: { display: true, text: 'O₂ Concentration (mol/m³)', color: '#4e4e5e', font: { size: 10, family: "'JetBrains Mono', monospace" } },
          grid: { color: 'rgba(255, 255, 255, 0.04)', drawBorder: false },
          min: 0,
        }
      },
      plugins: {
        legend: { position: 'top' },
        tooltip: {
          callbacks: {
            label: ctx => `${ctx.dataset.label}: ${ctx.parsed.y.toFixed(4)} mol/m³`
          }
        }
      }
    }
  });
}

function createAtpChart() {
  const ctx = document.getElementById('atp-chart').getContext('2d');
  atpChart = new Chart(ctx, {
    type: 'bar',
    data: {
      labels: ['OXPHOS', 'Glycolysis', 'Total Supply', 'Demand'],
      datasets: [
        {
          label: 'Baseline',
          data: [0, 0, 0, 0],
          backgroundColor: [
            'rgba(255, 255, 255, 0.12)',
            'rgba(255, 255, 255, 0.08)',
            'rgba(255, 255, 255, 0.18)',
            'rgba(255, 255, 255, 0.05)',
          ],
          borderColor: 'rgba(255, 255, 255, 0.2)',
          borderWidth: 1,
          borderRadius: 2,
        },
        {
          label: 'With MT',
          data: [0, 0, 0, 0],
          backgroundColor: [
            'rgba(255, 255, 255, 0.4)',
            'rgba(255, 255, 255, 0.25)',
            'rgba(255, 255, 255, 0.55)',
            'rgba(255, 255, 255, 0.05)',
          ],
          borderColor: 'rgba(255, 255, 255, 0.6)',
          borderWidth: 1,
          borderRadius: 2,
        }
      ]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        x: {
          grid: { color: 'rgba(255, 255, 255, 0.04)' },
        },
        y: {
          title: { display: true, text: 'μmol ATP/min/10⁶ cells', color: '#4e4e5e', font: { size: 10, family: "'JetBrains Mono', monospace" } },
          grid: { color: 'rgba(255, 255, 255, 0.04)' },
          min: 0,
        }
      },
      plugins: {
        legend: { position: 'top' },
        tooltip: {
          callbacks: {
            label: ctx => `${ctx.dataset.label}: ${ctx.parsed.y.toFixed(3)} μmol/min`
          }
        }
      }
    }
  });
}

function createBuyingTimeChart() {
  const ctx = document.getElementById('buying-time-chart').getContext('2d');
  buyingTimeChart = new Chart(ctx, {
    type: 'line',
    data: {
      labels: [],
      datasets: [
        {
          label: 'Baseline viability',
          data: [],
          borderColor: 'rgba(255, 255, 255, 0.25)',
          backgroundColor: 'rgba(255, 255, 255, 0.02)',
          fill: true,
          tension: 0.4,
          pointRadius: 0,
          borderWidth: 1.5,
          borderDash: [4, 2],
        },
        {
          label: 'With mito transplant',
          data: [],
          borderColor: 'rgba(255, 255, 255, 0.85)',
          backgroundColor: 'rgba(255, 255, 255, 0.05)',
          fill: true,
          tension: 0.4,
          pointRadius: 0,
          borderWidth: 2,
        },
        {
          label: 'Vascularization progress',
          data: [],
          borderColor: 'rgba(255, 255, 255, 0.45)',
          backgroundColor: 'rgba(255, 255, 255, 0.02)',
          fill: true,
          tension: 0.3,
          pointRadius: 0,
          borderWidth: 1.5,
          borderDash: [8, 4],
        }
      ]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      interaction: { mode: 'index', intersect: false },
      scales: {
        x: {
          title: { display: true, text: 'Days post-implantation', color: '#4e4e5e', font: { size: 10, family: "'JetBrains Mono', monospace" } },
          grid: { color: 'rgba(255, 255, 255, 0.04)' },
          ticks: {
            maxTicksLimit: 12,
            callback: function(val, idx) {
              const v = parseFloat(this.getLabelForValue(val));
              return v % 2 === 0 ? v : '';
            }
          },
        },
        y: {
          min: 0, max: 100,
          title: { display: true, text: '% viability / vascularization', color: '#4e4e5e', font: { size: 10, family: "'JetBrains Mono', monospace" } },
          grid: { color: 'rgba(255, 255, 255, 0.04)' },
        }
      },
      plugins: {
        legend: { position: 'top' },
        tooltip: {
          callbacks: {
            label: ctx => `${ctx.dataset.label}: ${Math.round(ctx.parsed.y)}%`
          }
        }
      }
    }
  });
}

function createViabilityChart() {
  const ctx = document.getElementById('viability-chart').getContext('2d');
  viabilityChart = new Chart(ctx, {
    type: 'line',
    data: {
      labels: [],
      datasets: [
        {
          label: 'Baseline @ Surface',
          data: [],
          borderColor: 'rgba(255, 255, 255, 0.15)',
          borderDash: [3, 3],
          pointRadius: 0,
          borderWidth: 1,
          fill: false,
        },
        {
          label: 'Baseline @ Center',
          data: [],
          borderColor: 'rgba(255, 255, 255, 0.3)',
          borderDash: [4, 2],
          pointRadius: 0,
          borderWidth: 1.5,
          fill: false,
        },
        {
          label: 'MT @ Surface',
          data: [],
          borderColor: 'rgba(255, 255, 255, 0.5)',
          borderDash: [3, 3],
          pointRadius: 0,
          borderWidth: 1,
          fill: false,
        },
        {
          label: 'MT @ Center',
          data: [],
          borderColor: 'rgba(255, 255, 255, 0.9)',
          backgroundColor: 'rgba(255, 255, 255, 0.04)',
          pointRadius: 0,
          borderWidth: 2.5,
          fill: true,
        },
        {
          label: 'Viability Threshold',
          data: [],
          borderColor: 'rgba(255, 255, 255, 0.12)',
          borderDash: [6, 4],
          pointRadius: 0,
          borderWidth: 1,
          fill: false,
        },
      ]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      interaction: { mode: 'index', intersect: false },
      scales: {
        x: {
          title: { display: true, text: 'Days', color: '#4e4e5e', font: { size: 10, family: "'JetBrains Mono', monospace" } },
          grid: { color: 'rgba(255, 255, 255, 0.04)' },
        },
        y: {
          title: { display: true, text: 'ATP Supply / Demand (%)', color: '#4e4e5e', font: { size: 10, family: "'JetBrains Mono', monospace" } },
          grid: { color: 'rgba(255, 255, 255, 0.04)' },
          min: 0,
          max: 150,
        }
      },
      plugins: {
        legend: { position: 'top' },
        tooltip: {
          callbacks: {
            label: ctx => `${ctx.dataset.label}: ${ctx.parsed.y.toFixed(1)}%`
          }
        }
      }
    }
  });
}

// ── Update Functions ──────────────────────────────────────

function getParams() {
  return {
    D: parseFloat(sliders.diffCoeff.value),
    C0: parseFloat(sliders.surfaceO2.value),
    rho: parseFloat(sliders.cellDensity.value),
    q: parseFloat(sliders.ocrPerCell.value),
    eta: parseFloat(sliders.mtEfficiency.value),
    atpDemand: parseFloat(sliders.atpDemand.value),
    poRatio: parseFloat(sliders.oxphosEfficiency.value),
    glycFrac: parseFloat(sliders.glycolysisRate.value) / 100,
    localO2Pct: parseFloat(sliders.localO2Pct.value) / 100,
    vascStrategy: sliders.vascStrategy.value,
    hypoxiaTolerance: parseFloat(sliders.hypoxiaTolerance.value),
    o2DepletionDays: parseFloat(sliders.o2DepletionDays.value),
    thickness: parseFloat(sliders.constructThickness.value),
    viabThreshold: parseFloat(sliders.viabilityThreshold.value) / 100,
    halflife: parseFloat(sliders.mtDecayHalflife.value),
    reDose: parseFloat(sliders.reDoseInterval.value),
  };
}

function updateDisplays(p) {
  displays.diffCoeff.textContent = p.D.toFixed(1);
  displays.surfaceO2.textContent = p.C0.toFixed(2);
  displays.cellDensity.textContent = p.rho.toFixed(1);
  displays.ocrPerCell.textContent = p.q.toFixed(0);
  displays.mtEfficiency.textContent = p.eta.toFixed(1) + '×';
  displays.atpDemand.textContent = p.atpDemand.toFixed(1);
  displays.oxphosEfficiency.textContent = p.poRatio.toFixed(1);
  displays.glycolysisRate.textContent = (p.glycFrac * 100).toFixed(0) + '%';
  displays.localO2Pct.textContent = (p.localO2Pct * 100).toFixed(0) + '%';
  displays.hypoxiaTolerance.textContent = p.hypoxiaTolerance.toFixed(0) + 'h';
  displays.o2DepletionDays.textContent = p.o2DepletionDays.toFixed(1);
  displays.constructThickness.textContent = p.thickness.toFixed(0);
  displays.viabilityThreshold.textContent = (p.viabThreshold * 100).toFixed(0) + '%';
  displays.mtDecayHalflife.textContent = p.halflife.toFixed(0);
  displays.reDoseInterval.textContent = p.reDose.toFixed(0);
}

function updateKroghChart(p) {
  const L_baseline = kroghPenetration(p.D, p.C0, p.rho, p.q, 1.0);
  const L_mt = kroghPenetration(p.D, p.C0, p.rho, p.q, p.eta);

  const maxX = Math.min(Math.max(L_mt * 1.3, L_baseline * 1.3, 500), 5000);
  const nPoints = 100;

  const labels = [];
  const baselineData = [];
  const mtData = [];
  const thresholdData = [];

  const hypoxiaThreshold = p.C0 * 0.05;

  for (let i = 0; i <= nPoints; i++) {
    const x = (i / nPoints) * maxX;
    labels.push(x.toFixed(0));
    baselineData.push(oxygenProfile(x, L_baseline, p.C0));
    mtData.push(oxygenProfile(x, L_mt, p.C0));
    thresholdData.push(hypoxiaThreshold);
  }

  kroghChart.data.labels = labels;
  kroghChart.data.datasets[0].data = baselineData;
  kroghChart.data.datasets[1].data = mtData;
  kroghChart.data.datasets[2].data = thresholdData;
  kroghChart.update('none');

  // Update result cards
  const depthEl = document.getElementById('penetration-depth-val');
  const foldEl = document.getElementById('fold-improvement-val');
  const volEl = document.getElementById('viable-volume-val');

  depthEl.textContent = L_mt >= 10000 ? '>10 mm' : (L_mt < 1000 ? L_mt.toFixed(0) + ' μm' : (L_mt / 1000).toFixed(2) + ' mm');
  document.getElementById('penetration-depth-sub').textContent = `Baseline: ${L_baseline.toFixed(0)} μm`;

  const fold = L_mt / L_baseline;
  foldEl.textContent = fold.toFixed(1) + '×';
  document.getElementById('fold-improvement-sub').textContent = `vs. baseline: ${L_baseline.toFixed(0)} μm`;

  const halfThickness = p.thickness / 2;
  const viableFracBaseline = Math.min(L_baseline / halfThickness, 1.0);
  const viableFracMt = Math.min(L_mt / halfThickness, 1.0);
  volEl.textContent = (viableFracMt * 100).toFixed(1) + '%';
  document.getElementById('viable-volume-sub').textContent = `Baseline: ${(viableFracBaseline * 100).toFixed(1)}%`;

  document.getElementById('penetration-result').className = 'result-card ' + (L_mt > 500 ? 'positive' : L_mt > 200 ? 'neutral' : 'negative');
  document.getElementById('fold-improvement-result').className = 'result-card ' + (fold >= 3 ? 'positive' : fold >= 1.5 ? 'neutral' : 'negative');
  document.getElementById('viable-volume-result').className = 'result-card ' + (viableFracMt >= 0.8 ? 'positive' : viableFracMt >= 0.4 ? 'neutral' : 'negative');

  return { L_baseline, L_mt, fold, viableFracBaseline, viableFracMt };
}

function updateAtpChart(p) {
  const demand = p.atpDemand;
  const localO2 = p.localO2Pct;

  const oxphosBase = (1 - p.glycFrac) * demand * localO2 * (p.poRatio / 2.5);
  const glycBase = p.glycFrac * demand;
  const totalBase = oxphosBase + glycBase;

  const oxphosMt = (1 - p.glycFrac) * demand * localO2 * (p.poRatio / 2.5) * p.eta;
  const glycMt = p.glycFrac * demand;
  const totalMt = oxphosMt + glycMt;

  atpChart.data.datasets[0].data = [oxphosBase, glycBase, totalBase, demand];
  atpChart.data.datasets[1].data = [oxphosMt, glycMt, totalMt, demand];
  atpChart.update('none');

  const supplyEl = document.getElementById('atp-supply-val');
  const deficitEl = document.getElementById('atp-deficit-val');
  const viabEl = document.getElementById('viability-val');

  supplyEl.textContent = totalMt.toFixed(2);
  document.getElementById('atp-supply-sub').textContent = `Baseline: ${totalBase.toFixed(2)} μmol/min`;

  const balance = totalMt - demand;
  const balancePct = (totalMt / demand * 100).toFixed(0);
  deficitEl.textContent = (balance >= 0 ? '+' : '') + balance.toFixed(2);
  document.getElementById('atp-deficit-sub').textContent = `${balancePct}% of demand met`;

  let viabStatus, viabColor;
  if (totalMt / demand >= 0.8) { viabStatus = 'Viable'; viabColor = 'positive'; }
  else if (totalMt / demand >= 0.5) { viabStatus = 'Stressed'; viabColor = 'neutral'; }
  else if (totalMt / demand >= 0.2) { viabStatus = 'Critical'; viabColor = 'negative'; }
  else { viabStatus = 'Non-viable'; viabColor = 'negative'; }

  viabEl.textContent = viabStatus;
  document.getElementById('viability-sub').textContent = `${balancePct}% demand covered`;

  document.getElementById('atp-supply-result').className = 'result-card ' + (totalMt >= totalBase * 1.5 ? 'positive' : 'neutral');
  document.getElementById('atp-deficit-result').className = 'result-card ' + (balance >= 0 ? 'positive' : 'negative');
  document.getElementById('viability-result').className = 'result-card ' + viabColor;

  return { totalBase, totalMt, demand, balance };
}

/**
 * "Buying Time" model — matches buying_time_survival_model.html logic
 * 
 * Core idea: cells deplete local O2 by day ~5 (o2DepletionDays).
 * After that, they survive on hypoxia tolerance for hypTolDays.
 * MT boost extends BOTH phases linearly (same as original model).
 * Vascularization races to complete before cells die.
 */
function updateBuyingTimeChart(p) {
  const vasc = vascData[p.vascStrategy];
  const boost = p.eta;
  const hypTolDays = p.hypoxiaTolerance / 24;
  const o2Depletion = p.o2DepletionDays;

  // Critical day = when cells die (O2 depletion + hypoxia tolerance period)
  const baseCriticalDay = o2Depletion + hypTolDays;
  // MT extends both the pre-depletion phase and hypoxia tolerance
  const mitoCriticalDay = (o2Depletion * boost) + (hypTolDays * boost);

  const maxDay = Math.max(20, Math.ceil(mitoCriticalDay) + 2, vasc.complete + 3);
  const labels = [];
  const baseViab = [];
  const mitoViab = [];
  const vascProg = [];

  for (let d = 0; d <= maxDay; d += 0.5) {
    labels.push(d.toFixed(1));

    // Sigmoid cell death curves (matching original model exactly)
    const decayStart = o2Depletion;
    const decayStartMito = o2Depletion * boost;

    const bv = d < decayStart ? 100 : Math.max(0, 100 * sigmoid(d, baseCriticalDay, 1.2));
    const mv = d < decayStartMito ? 100 : Math.max(0, 100 * sigmoid(d, mitoCriticalDay, 1.2 / boost));

    baseViab.push(Math.round(bv * 10) / 10);
    mitoViab.push(Math.round(mv * 10) / 10);
    vascProg.push(Math.round(vascCurve(d, vasc.onset, vasc.complete) * 10) / 10);
  }

  buyingTimeChart.data.labels = labels;
  buyingTimeChart.data.datasets[0].data = baseViab;
  buyingTimeChart.data.datasets[1].data = mitoViab;
  buyingTimeChart.data.datasets[2].data = vascProg;
  buyingTimeChart.update('none');

  // Result cards
  const baseDays = baseCriticalDay.toFixed(1);
  const mitoDays = mitoCriticalDay.toFixed(1);
  const gained = (mitoCriticalDay - baseCriticalDay).toFixed(1);

  document.getElementById('base-survival-val').textContent = baseDays + ' d';
  document.getElementById('base-survival-sub').textContent = 'O₂ depletion + hypoxia window';

  document.getElementById('mt-survival-val').textContent = mitoDays + ' d';
  document.getElementById('mt-survival-sub').textContent = '+' + gained + ' days gained';

  const baseBridges = baseCriticalDay >= vasc.complete;
  const mitoBridges = mitoCriticalDay >= vasc.complete;

  const bridgeEl = document.getElementById('bridge-gap-val');
  const bridgeSubEl = document.getElementById('bridge-gap-sub');

  if (mitoBridges && !baseBridges) {
    bridgeEl.textContent = '✓ Yes!';
    bridgeSubEl.textContent = `MT bridges the gap to ${vasc.label}`;
    document.getElementById('bridge-gap-result').className = 'result-card positive';
  } else if (mitoBridges && baseBridges) {
    bridgeEl.textContent = '✓ Safety';
    bridgeSubEl.textContent = `+${gained}d margin added`;
    document.getElementById('bridge-gap-result').className = 'result-card positive';
  } else {
    const gap = (vasc.complete - mitoCriticalDay).toFixed(1);
    bridgeEl.textContent = '✗ No';
    bridgeSubEl.textContent = `${gap}d gap remains`;
    document.getElementById('bridge-gap-result').className = 'result-card negative';
  }

  document.getElementById('base-survival-result').className = 'result-card ' + (baseBridges ? 'neutral' : 'negative');
  document.getElementById('mt-survival-result').className = 'result-card ' + (mitoBridges ? 'positive' : 'neutral');

  return { baseCriticalDay, mitoCriticalDay, baseBridges, mitoBridges, vasc, gained };
}

function updateViabilityChart(p) {
  const days = 30;
  const nPoints = 60;
  const labels = [];

  const dsBaseSurf = [];
  const dsBaseCenter = [];
  const dsMtSurf = [];
  const dsMtCenter = [];
  const dsThreshold = [];

  const halfThickness = p.thickness / 2;

  for (let i = 0; i <= nPoints; i++) {
    const t = (i / nPoints) * days;
    labels.push(t.toFixed(1));

    const eta_t = mtEfficiencyOverTime(t, p.eta, p.halflife, p.reDose);

    const L_base = kroghPenetration(p.D, p.C0, p.rho, p.q, 1.0);
    const L_mt_t = kroghPenetration(p.D, p.C0, p.rho, p.q, eta_t);

    const o2FracCenterBase = halfThickness >= L_base ? 0 : (1 - (halfThickness * halfThickness) / (L_base * L_base));
    const o2FracCenterMt = halfThickness >= L_mt_t ? 0 : (1 - (halfThickness * halfThickness) / (L_mt_t * L_mt_t));

    const supplySurfBase = atpSupplyAtDepth(1.0, 1.0, p.poRatio, p.glycFrac, p.atpDemand);
    const supplyCenterBase = atpSupplyAtDepth(o2FracCenterBase, 1.0, p.poRatio, p.glycFrac, p.atpDemand);
    const supplySurfMt = atpSupplyAtDepth(1.0, eta_t, p.poRatio, p.glycFrac, p.atpDemand);
    const supplyCenterMt = atpSupplyAtDepth(o2FracCenterMt, eta_t, p.poRatio, p.glycFrac, p.atpDemand);

    dsBaseSurf.push(Math.min((supplySurfBase / p.atpDemand) * 100, 150));
    dsBaseCenter.push(Math.min((supplyCenterBase / p.atpDemand) * 100, 150));
    dsMtSurf.push(Math.min((supplySurfMt / p.atpDemand) * 100, 150));
    dsMtCenter.push(Math.min((supplyCenterMt / p.atpDemand) * 100, 150));
    dsThreshold.push(p.viabThreshold * 100);
  }

  viabilityChart.data.labels = labels;
  viabilityChart.data.datasets[0].data = dsBaseSurf;
  viabilityChart.data.datasets[1].data = dsBaseCenter;
  viabilityChart.data.datasets[2].data = dsMtSurf;
  viabilityChart.data.datasets[3].data = dsMtCenter;
  viabilityChart.data.datasets[4].data = dsThreshold;
  viabilityChart.update('none');

  let daysViableMt = 0;
  let daysViableBase = 0;
  for (let i = 0; i <= nPoints; i++) {
    const t = (i / nPoints) * days;
    if (dsMtCenter[i] >= p.viabThreshold * 100) daysViableMt = t;
    if (dsBaseCenter[i] >= p.viabThreshold * 100) daysViableBase = t;
  }

  const daysEl = document.getElementById('days-viable-val');
  daysEl.textContent = daysViableMt >= days ? '>30' : daysViableMt.toFixed(1);
  document.getElementById('days-viable-sub').textContent = `Baseline: ${daysViableBase.toFixed(1)} days`;

  const targetEl = document.getElementById('target-met-val');
  const meetsTarget = daysViableMt >= 15;
  targetEl.textContent = meetsTarget ? '✓ Yes' : '✗ No';
  document.getElementById('target-met-sub').textContent = meetsTarget
    ? `${(daysViableMt - 15).toFixed(1)} days margin`
    : `${(15 - daysViableMt).toFixed(1)} days short`;

  document.getElementById('days-viable-result').className = 'result-card ' + (daysViableMt >= 15 ? 'positive' : daysViableMt >= 7 ? 'neutral' : 'negative');
  document.getElementById('target-met-result').className = 'result-card ' + (meetsTarget ? 'positive' : 'negative');

  return { daysViableMt, daysViableBase, meetsTarget };
}

function updateVerdict(kroghResults, atpResults, buyingTimeResults, viabilityResults, p) {
  const vh = document.getElementById('verdict-header');
  const vb = document.getElementById('verdict-body');
  const vr = document.getElementById('verdict-recommendations');

  const { fold, L_mt, viableFracMt } = kroghResults;
  const { totalMt, demand, balance } = atpResults;
  const { mitoBridges, baseBridges, mitoCriticalDay, baseCriticalDay, vasc, gained } = buyingTimeResults;
  const { daysViableMt, meetsTarget } = viabilityResults;

  // Determine overall verdict
  let verdictClass, verdictTitle, verdictDesc;

  if (fold >= 3 && meetsTarget && balance >= 0) {
    verdictClass = 'positive';
    verdictTitle = '✅ The Math Supports Your Hypothesis';
    verdictDesc = 'With these parameters, mitochondrial transplantation meaningfully shifts the metabolic demand and extends cell survival in the construct.';
  } else if (mitoBridges && !baseBridges) {
    verdictClass = 'positive';
    verdictTitle = '✅ MT Bridges the Vascularization Gap';
    verdictDesc = `Mitochondrial transplantation extends survival from ${baseCriticalDay.toFixed(1)} to ${mitoCriticalDay.toFixed(1)} days — enough to bridge the gap until ${vasc.label.toLowerCase()} completes at day ${vasc.complete}. This supports the "buying time" argument.`;
  } else if (fold >= 1.5 || (meetsTarget && balance >= -0.5) || mitoBridges) {
    verdictClass = 'partial';
    verdictTitle = '⚠️ Partially Promising — Multi-pronged Approach Needed';
    verdictDesc = 'Mitochondrial transplantation offers meaningful improvement, but would benefit from combination with vascularization strategies for optimal results.';
  } else {
    verdictClass = 'negative';
    verdictTitle = '❌ The Math Is Challenging';
    verdictDesc = 'Under these parameters, mitochondrial transplantation alone is unlikely to produce sufficient metabolic shift. Consider it as one component of a multi-pronged strategy.';
  }

  vh.className = 'verdict-header ' + verdictClass;
  vh.innerHTML = `<h3>${verdictTitle}</h3><p>${verdictDesc}</p>`;

  // Metrics
  vb.innerHTML = `
    <div class="verdict-metric">
      <span class="verdict-metric-label">Penetration depth improvement</span>
      <span class="verdict-metric-value ${fold >= 3 ? 'good' : fold >= 1.5 ? 'warn' : 'bad'}">${fold.toFixed(1)}×</span>
    </div>
    <div class="verdict-metric">
      <span class="verdict-metric-label">√η scaling relationship</span>
      <span class="verdict-metric-value warn">L ∝ √η — depth improvement is square root of efficiency gain</span>
    </div>
    <div class="verdict-metric">
      <span class="verdict-metric-label">Depth with current MT settings</span>
      <span class="verdict-metric-value ${L_mt > 500 ? 'good' : L_mt > 200 ? 'warn' : 'bad'}">${L_mt.toFixed(0)} μm</span>
    </div>
    <div class="verdict-metric">
      <span class="verdict-metric-label">ATP balance at set depth</span>
      <span class="verdict-metric-value ${balance >= 0 ? 'good' : 'bad'}">${balance >= 0 ? '+' : ''}${balance.toFixed(2)} μmol/min</span>
    </div>
    <div class="verdict-metric">
      <span class="verdict-metric-label">Buying time: baseline → MT survival</span>
      <span class="verdict-metric-value ${mitoBridges ? 'good' : 'warn'}">${baseCriticalDay.toFixed(1)}d → ${mitoCriticalDay.toFixed(1)}d (+${gained}d)</span>
    </div>
    <div class="verdict-metric">
      <span class="verdict-metric-label">Bridges ${vasc.label} gap? (completes day ${vasc.complete})</span>
      <span class="verdict-metric-value ${mitoBridges ? 'good' : 'bad'}">${mitoBridges ? 'Yes ✓' : 'No ✗'}</span>
    </div>
    <div class="verdict-metric">
      <span class="verdict-metric-label">Days viable at center (depth-resolved)</span>
      <span class="verdict-metric-value ${meetsTarget ? 'good' : daysViableMt >= 7 ? 'warn' : 'bad'}">${daysViableMt >= 30 ? '>30' : daysViableMt.toFixed(1)} days (target: 15)</span>
    </div>
    <div class="verdict-metric">
      <span class="verdict-metric-label">Viable volume fraction</span>
      <span class="verdict-metric-value ${viableFracMt >= 0.8 ? 'good' : viableFracMt >= 0.4 ? 'warn' : 'bad'}">${(viableFracMt * 100).toFixed(1)}%</span>
    </div>
  `;

  // Recommendations
  let recs = '';

  recs += `
    <div class="recommendation-item">
      <span class="rec-icon caution">!</span>
      <span><strong>Key constraint:</strong> Penetration depth scales as √η, meaning a 2× efficiency boost yields only ~1.4× depth improvement. Literature suggests MT provides 1.2–2× improvement in metabolic efficiency.</span>
    </div>
  `;

  if (mitoBridges && !baseBridges) {
    recs += `
      <div class="recommendation-item">
        <span class="rec-icon strategy">→</span>
        <span><strong>The "buying time" argument works!</strong> With ${vasc.label.toLowerCase()}, baseline cells die at day ~${baseCriticalDay.toFixed(1)} but MT-treated cells survive to day ~${mitoCriticalDay.toFixed(1)}, successfully bridging until vascularization at day ${vasc.complete}.</span>
      </div>
    `;
  } else if (!mitoBridges) {
    recs += `
      <div class="recommendation-item">
        <span class="rec-icon strategy">→</span>
        <span><strong>Frame as adjuvant therapy:</strong> Position MT not as a standalone solution but as one leg of a multi-pronged strategy: (1) MT for metabolic resilience, (2) pre-vascularization channels for O₂ delivery, (3) oxygen-releasing biomaterials for bridging the avascular period. Try switching to "pre-vascularized" strategy above to see if the gap closes.</span>
      </div>
    `;
  }

  recs += `
    <div class="recommendation-item">
      <span class="rec-icon experiment">⬡</span>
      <span><strong>Seahorse XF experiment:</strong> Measure OCR (pmol/min) in hiPSC-CMs before and after MT under normoxic and hypoxic conditions. This gives you the real η value. Literature ranges 1.2–2.0× for basal respiration improvement.</span>
    </div>
  `;

  recs += `
    <div class="recommendation-item">
      <span class="rec-icon experiment">⬡</span>
      <span><strong>Hypoxia tolerance assay:</strong> Culture cells with and without MT at different O₂ levels (1%, 5%, 21%) and measure viability at 3, 7, 10, 15 days. This tests whether MT extends the survival window even without solving the diffusion limit per se.</span>
    </div>
  `;

  recs += `
    <div class="recommendation-item">
      <span class="rec-icon strategy">→</span>
      <span><strong>Proposed framing:</strong> "The Krogh model shows penetration depth scales as √η. A ${p.eta.toFixed(1)}× efficiency boost extends the viable zone by ~${((fold - 1) * 100).toFixed(0)}%, and MT extends cell survival from ${baseCriticalDay.toFixed(1)} to ${mitoCriticalDay.toFixed(1)} days — ${mitoBridges ? 'successfully bridging' : 'narrowing the gap to'} the ${vasc.label.toLowerCase()} window. We propose MT as a metabolic resilience layer within a vascularized construct."</span>
    </div>
  `;

  vr.innerHTML = `<h4>Recommendations & Next Steps</h4>${recs}`;
}

function updateAll() {
  const p = getParams();
  updateDisplays(p);
  const kroghResults = updateKroghChart(p);
  const atpResults = updateAtpChart(p);
  const buyingTimeResults = updateBuyingTimeChart(p);
  const viabilityResults = updateViabilityChart(p);
  updateVerdict(kroghResults, atpResults, buyingTimeResults, viabilityResults, p);
}

// ── Initialize ────────────────────────────────────────────
function init() {
  createKroghChart();
  createAtpChart();
  createBuyingTimeChart();
  createViabilityChart();
  updateAll();

  // Scroll reveal for sections
  const sectionObserver = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        entry.target.classList.add('visible');
        sectionObserver.unobserve(entry.target);
      }
    });
  }, { threshold: 0.08, rootMargin: '0px 0px -40px 0px' });

  document.querySelectorAll('.calc-section').forEach(el => sectionObserver.observe(el));

  // Cell type dropdown
  sliders.cellType.addEventListener('change', handleCellTypeChange);

  // All sliders and selects
  Object.values(sliders).forEach(el => {
    if (el.tagName === 'SELECT') {
      el.addEventListener('change', updateAll);
    } else {
      el.addEventListener('input', updateAll);
    }
  });
}

document.addEventListener('DOMContentLoaded', init);
