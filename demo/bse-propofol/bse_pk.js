/* BSE Demo — Propofol (Marsh 3-compartment PK)
   RK4 integration of:
   dA1/dt = Dose(t) + k21*A2 + k31*A3 - (k10+k12+k13)*A1
   dA2/dt = k12*A1 - k21*A2
   dA3/dt = k13*A1 - k31*A3
   Plasma concentration: C1 = A1 / (V1_per_kg * weight)
*/

document.addEventListener("DOMContentLoaded", () => {
  const ids = [
    "weight", "bolus_mg",
    "infusion_mg_min", "infusion_duration_min",
    "t_end_min", "dt_min"
  ];

  const inputs = {};
  const displays = {};

  ids.forEach(id => {
    inputs[id] = document.getElementById(id);
    displays[id] = document.getElementById("val_" + id);
    inputs[id].addEventListener("input", update);
  });

  // --- Marsh constants (min^-1) ---
  const k10 = 0.119;
  const k12 = 0.112;
  const k21 = 0.055;
  const k13 = 0.042;
  const k31 = 0.0033;

  // Volumes (L/kg)
  const V1_per_kg = 0.228;
  const V2_per_kg = 0.463;
  const V3_per_kg = 2.893;

  function readParams() {
    const p = {};
    ids.forEach(id => {
      p[id] = parseFloat(inputs[id].value);
      displays[id].innerText = (Number.isFinite(p[id]) ? p[id] : 0);
    });

    // Guardrails
    if (!Number.isFinite(p.weight) || p.weight <= 0) p.weight = 70;
    if (!Number.isFinite(p.bolus_mg) || p.bolus_mg < 0) p.bolus_mg = 0;

    if (!Number.isFinite(p.infusion_mg_min) || p.infusion_mg_min < 0) p.infusion_mg_min = 0;
    if (!Number.isFinite(p.infusion_duration_min) || p.infusion_duration_min < 0) p.infusion_duration_min = 0;

    if (!Number.isFinite(p.t_end_min) || p.t_end_min <= 0) p.t_end_min = 240;
    if (!Number.isFinite(p.dt_min) || p.dt_min <= 0) p.dt_min = 0.2;

    // Cap max steps to avoid browser freeze
    const MAX_STEPS = 20000;
    const steps = p.t_end_min / p.dt_min;
    if (steps > MAX_STEPS) {
      p.dt_min = p.t_end_min / MAX_STEPS;
      // reflect adjusted dt in UI (without triggering loops)
      inputs.dt_min.value = p.dt_min.toFixed(4);
      displays.dt_min.innerText = inputs.dt_min.value;
    }

    return p;
  }

  // Infusion (mg/min) applied from t=0 to t=infusion_duration
  function doseRate(t, p) {
    if (p.infusion_mg_min <= 0 || p.infusion_duration_min <= 0) return 0;
    return (t <= p.infusion_duration_min) ? p.infusion_mg_min : 0;
  }

  // derivatives
  function deriv(t, A1, A2, A3, p) {
    const input = doseRate(t, p);
    const dA1 = input + k21*A2 + k31*A3 - (k10 + k12 + k13)*A1;
    const dA2 = k12*A1 - k21*A2;
    const dA3 = k13*A1 - k31*A3;
    return [dA1, dA2, dA3];
  }

  // RK4 solver
  function solveRK4(p) {
    let t = 0;

    // Initial conditions:
    // instantaneous bolus into central compartment
    let A1 = p.bolus_mg;
    let A2 = 0;
    let A3 = 0;

    const arr_t = [t];
    const arr_A1 = [A1];
    const arr_A2 = [A2];
    const arr_A3 = [A3];

    const nSteps = Math.ceil(p.t_end_min / p.dt_min);

    for (let i = 0; i < nSteps; i++) {
      if (t >= p.t_end_min) break;
      const h = Math.min(p.dt_min, p.t_end_min - t);

      const k1 = deriv(t, A1, A2, A3, p);
      const k2 = deriv(t + 0.5*h, A1 + 0.5*h*k1[0], A2 + 0.5*h*k1[1], A3 + 0.5*h*k1[2], p);
      const k3 = deriv(t + 0.5*h, A1 + 0.5*h*k2[0], A2 + 0.5*h*k2[1], A3 + 0.5*h*k2[2], p);
      const k4 = deriv(t + h,     A1 + h*k3[0],     A2 + h*k3[1],     A3 + h*k3[2],     p);

      A1 += (h/6) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
      A2 += (h/6) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
      A3 += (h/6) * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);

      // Avoid tiny negatives from numerical error
      A1 = Math.max(0, A1);
      A2 = Math.max(0, A2);
      A3 = Math.max(0, A3);

      t += h;

      arr_t.push(t);
      arr_A1.push(A1);
      arr_A2.push(A2);
      arr_A3.push(A3);
    }

    // compute C1
    const V1 = V1_per_kg * p.weight; // L
    const arr_C1 = arr_A1.map(a => a / V1); // mg/L

    return { t: arr_t, A1: arr_A1, A2: arr_A2, A3: arr_A3, C1: arr_C1 };
  }

  function computeMetrics(sol) {
    // Metrics on C1
    let cmax = -Infinity;
    let tmax = 0;
    let auc = 0;

    for (let i = 0; i < sol.C1.length; i++) {
      const c = sol.C1[i];
      if (c > cmax) {
        cmax = c;
        tmax = sol.t[i];
      }
      if (i > 0) {
        const dt = sol.t[i] - sol.t[i - 1];
        auc += 0.5 * (sol.C1[i] + sol.C1[i - 1]) * dt;
      }
    }

    document.getElementById("out_cmax").innerText = (cmax >= 0 ? cmax.toFixed(3) : "—");
    document.getElementById("out_tmax").innerText = tmax.toFixed(2);
    document.getElementById("out_auc").innerText = auc.toFixed(2);
  }

  function renderPlots(sol) {
    // Plot 1: C1
    const traceC1 = {
      x: sol.t,
      y: sol.C1,
      mode: "lines",
      name: "Plasma C1 (mg/L)"
    };

    const layoutC1 = {
      title: "Central (Plasma) Concentration — C1",
      xaxis: { title: "Time (min)", rangemode: "tozero" },
      yaxis: { title: "C1 (mg/L)", rangemode: "tozero" },
      margin: { l: 60, r: 20, t: 50, b: 50 },
      hovermode: "x unified",
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: "#e5e7eb" }
    };

    Plotly.react("plot_c1", [traceC1], layoutC1, { responsive: true, displaylogo: false });

    // Plot 2: Amounts
    const trA1 = { x: sol.t, y: sol.A1, mode: "lines", name: "A1 (Central, mg)" };
    const trA2 = { x: sol.t, y: sol.A2, mode: "lines", name: "A2 (Fast peripheral, mg)" };
    const trA3 = { x: sol.t, y: sol.A3, mode: "lines", name: "A3 (Slow peripheral, mg)" };

    const layoutA = {
      title: "Amounts by Compartment",
      xaxis: { title: "Time (min)", rangemode: "tozero" },
      yaxis: { title: "Amount (mg)", rangemode: "tozero" },
      margin: { l: 60, r: 20, t: 50, b: 50 },
      hovermode: "x unified",
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: "#e5e7eb" }
    };

    Plotly.react("plot_amounts", [trA1, trA2, trA3], layoutA, { responsive: true, displaylogo: false });
  }

  function update() {
    const p = readParams();
    const sol = solveRK4(p);
    computeMetrics(sol);
    renderPlots(sol);
  }

  // First render
  update();
});
