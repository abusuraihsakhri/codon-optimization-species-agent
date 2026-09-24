"use strict";

const $ = (id) => document.getElementById(id);

const SAMPLE = {
  dna: "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCAC",
  protein: "MASKGEELFTGVVPILVELDGDVNGH",
  organism: "e_coli",
};

const PYODIDE_INDEX = "https://cdn.jsdelivr.net/pyodide/v314.0.7/full/";

const BRIDGE_SOURCE = `
import json

def browser_analyze(dna, protein, organism, rare_threshold):
    if organism not in ORGANISM_RSCU:
        raise ValueError("Unknown target organism")

    seq = _normalize_nucleotide_sequence(dna, require_codon_aligned=True)
    if not seq:
        raise ValueError("Enter a coding DNA or RNA sequence")

    rscu = ORGANISM_RSCU[organism]
    translated = "".join(
        CODON_TABLE[seq[i:i + 3]]
        for i in range(0, len(seq), 3)
    )

    supplied_protein = _normalize_protein_sequence(protein) if protein.strip() else translated
    cai = calculate_cai(seq, rscu)
    gc = analyze_gc_content(seq)
    rare = identify_rare_codons(seq, rscu, float(rare_threshold))
    hairpins = detect_hairpins(seq)
    optimized = full_optimization(supplied_protein, seq, organism)

    payload = {
        "organism": organism,
        "normalized_dna": seq,
        "protein": supplied_protein,
        "translation": translated,
        "cai": cai.cai,
        "num_cai_codons": cai.num_codons,
        "gc_content": gc.gc_content,
        "gc_positions": gc.gc_at_positions,
        "rare_threshold": float(rare_threshold),
        "rare_codons": [
            {
                "position": item.position,
                "codon": item.codon,
                "amino_acid": item.amino_acid,
                "rscu": item.rscu,
                "frequency": item.frequency,
            }
            for item in rare
        ],
        "hairpins": [
            {
                "position": item.position,
                "stem_length": item.stem_length,
                "pseudo_energy": item.free_energy,
                "sequence": item.sequence,
            }
            for item in hairpins
        ],
        "optimized_sequence": optimized.optimized_sequence,
        "optimized_cai": optimized.optimized_cai,
        "optimized_gc": optimized.optimized_gc,
        "num_changes": optimized.num_changes,
        "changes": [
            {"position": p, "original": old, "optimized": new}
            for p, old, new in optimized.changes
        ],
    }
    return json.dumps(payload)
`;

let pyodide = null;
let latestResult = null;

function setRuntimeStatus(text, state) {
  const el = $("runtimeStatus");
  el.textContent = text;
  el.className = "status-pill " + state;
}

function setError(message) {
  const box = $("inputError");
  box.textContent = message;
  box.hidden = !message;
}

function normalizeVisibleSequence(value) {
  return value.replace(/\s+/g, "").toUpperCase();
}

function updateDnaMeta() {
  const seq = normalizeVisibleSequence($("dnaSequence").value);
  $("dnaMeta").textContent = seq.length + " nt · " + Math.floor(seq.length / 3) + " codons";
}

function setMobilePanel(panelId) {
  document.querySelectorAll(".mobile-tab").forEach((button) => {
    button.classList.toggle("active", button.dataset.panel === panelId);
  });
  ["inputPanel", "resultsPanel"].forEach((id) => {
    $(id).classList.toggle("mobile-hidden", id !== panelId);
  });
}

function formatFraction(value) {
  return (Number(value) * 100).toFixed(1) + "%";
}

function clearRareTable(message) {
  const tbody = $("rareTableBody");
  tbody.replaceChildren();
  const row = document.createElement("tr");
  const cell = document.createElement("td");
  cell.colSpan = 4;
  cell.className = "empty-cell";
  cell.textContent = message;
  row.appendChild(cell);
  tbody.appendChild(row);
}

function renderRareTable(items) {
  const tbody = $("rareTableBody");
  tbody.replaceChildren();

  if (!items.length) {
    clearRareTable("No codons below the selected threshold.");
    return;
  }

  items.slice(0, 60).forEach((item) => {
    const row = document.createElement("tr");
    [item.codon, item.amino_acid, String(item.position + 1), Number(item.rscu).toFixed(3)].forEach((value) => {
      const cell = document.createElement("td");
      cell.textContent = value;
      row.appendChild(cell);
    });
    tbody.appendChild(row);
  });
}

function renderResult(data) {
  latestResult = data;

  const hostLabels = {
    e_coli: "E. coli",
    human: "Human",
    yeast: "S. cerevisiae",
  };

  $("resultTitle").textContent = hostLabels[data.organism] + " reference analysis";
  $("caiValue").textContent = Number(data.cai).toFixed(3);
  const delta = Number(data.optimized_cai) - Number(data.cai);
  $("caiDelta").textContent = "Optimized Δ " + (delta >= 0 ? "+" : "") + delta.toFixed(3);

  $("gcValue").textContent = formatFraction(data.gc_content);
  $("gcPositions").textContent =
    "GC1 " + formatFraction(data.gc_positions["1"] ?? data.gc_positions[1]) +
    " · GC2 " + formatFraction(data.gc_positions["2"] ?? data.gc_positions[2]) +
    " · GC3 " + formatFraction(data.gc_positions["3"] ?? data.gc_positions[3]);

  $("rareValue").textContent = String(data.rare_codons.length);
  $("rareDetail").textContent = "RSCU < " + Number(data.rare_threshold).toFixed(2);

  if (data.hairpins.length) {
    const hp = data.hairpins[0];
    $("hairpinValue").textContent = "Detected";
    $("hairpinDetail").textContent =
      "Stem " + hp.stem_length + " bp · pseudo-score " + Number(hp.pseudo_energy).toFixed(1);
  } else {
    $("hairpinValue").textContent = "None";
    $("hairpinDetail").textContent = "No qualifying inverted repeat";
  }

  $("optimizedSummary").textContent = data.optimized_sequence.length + " nt · " + data.protein.length + " aa";
  $("optimizedSequence").textContent = data.optimized_sequence;
  $("optimizedCai").textContent = "Optimized CAI " + Number(data.optimized_cai).toFixed(3);
  $("optimizedGc").textContent = "Optimized GC " + formatFraction(data.optimized_gc);
  $("changeCount").textContent = "Changes " + data.num_changes;

  $("rareHeading").textContent =
    data.rare_codons.length + (data.rare_codons.length === 1 ? " codon" : " codons") + " below threshold";
  renderRareTable(data.rare_codons);

  $("copyButton").disabled = false;
  $("downloadButton").disabled = false;

  if (window.matchMedia("(max-width: 980px)").matches) {
    setMobilePanel("resultsPanel");
  }
}

async function analyze() {
  setError("");
  const button = $("analyzeButton");
  const dna = $("dnaSequence").value;
  const protein = $("proteinSequence").value;
  const organism = $("organism").value;
  const threshold = Number($("threshold").value);

  if (!normalizeVisibleSequence(dna)) {
    setError("Enter a coding DNA or RNA sequence.");
    return;
  }
  if (!Number.isFinite(threshold) || threshold < 0) {
    setError("Rare-codon threshold must be a non-negative number.");
    return;
  }

  button.disabled = true;
  button.querySelector("span").textContent = "Analyzing…";

  try {
    pyodide.globals.set("js_dna", dna);
    pyodide.globals.set("js_protein", protein);
    pyodide.globals.set("js_organism", organism);
    pyodide.globals.set("js_threshold", threshold);

    const raw = await pyodide.runPythonAsync(
      "browser_analyze(js_dna, js_protein, js_organism, js_threshold)"
    );
    renderResult(JSON.parse(raw));
  } catch (error) {
    const message = String(error?.message || error).replace(/^PythonError:\s*/, "").split("\n")[0];
    setError(message);
    if (window.matchMedia("(max-width: 980px)").matches) {
      setMobilePanel("inputPanel");
    }
  } finally {
    button.querySelector("span").textContent = "Analyze sequence";
    button.disabled = false;
  }
}

async function copyOptimized() {
  if (!latestResult) return;
  const text = latestResult.optimized_sequence;
  try {
    await navigator.clipboard.writeText(text);
  } catch {
    const temp = document.createElement("textarea");
    temp.value = text;
    temp.style.position = "fixed";
    temp.style.opacity = "0";
    document.body.appendChild(temp);
    temp.select();
    document.execCommand("copy");
    temp.remove();
  }
  const button = $("copyButton");
  const original = button.textContent;
  button.textContent = "Copied";
  window.setTimeout(() => { button.textContent = original; }, 1200);
}

function downloadResult() {
  if (!latestResult) return;
  const blob = new Blob([JSON.stringify(latestResult, null, 2)], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = "codon-analysis.json";
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

function clearWorkspace() {
  $("dnaSequence").value = "";
  $("proteinSequence").value = "";
  $("organism").value = "e_coli";
  $("threshold").value = "0.30";
  latestResult = null;
  setError("");
  updateDnaMeta();

  $("resultTitle").textContent = "Ready for analysis";
  $("caiValue").textContent = "—";
  $("caiDelta").textContent = "Reference-adapted index";
  $("gcValue").textContent = "—";
  $("gcPositions").textContent = "GC1 — · GC2 — · GC3 —";
  $("rareValue").textContent = "—";
  $("rareDetail").textContent = "Below selected RSCU threshold";
  $("hairpinValue").textContent = "—";
  $("hairpinDetail").textContent = "Sequence-only screen";
  $("optimizedSummary").textContent = "No result yet";
  $("optimizedSequence").textContent = "Run an analysis to generate an RSCU-based sequence.";
  $("optimizedCai").textContent = "Optimized CAI —";
  $("optimizedGc").textContent = "Optimized GC —";
  $("changeCount").textContent = "Changes —";
  $("rareHeading").textContent = "No result yet";
  clearRareTable("No analysis yet.");
  $("copyButton").disabled = true;
  $("downloadButton").disabled = true;
}

function loadSample() {
  $("dnaSequence").value = SAMPLE.dna;
  $("proteinSequence").value = SAMPLE.protein;
  $("organism").value = SAMPLE.organism;
  updateDnaMeta();
  setError("");
}

function initTheme() {
  const saved = localStorage.getItem("codon-theme");
  const preferredDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
  const theme = saved || (preferredDark ? "dark" : "light");
  document.documentElement.dataset.theme = theme;
  $("themeToggle").setAttribute("aria-label", theme === "dark" ? "Switch to light theme" : "Switch to dark theme");
}

function toggleTheme() {
  const next = document.documentElement.dataset.theme === "dark" ? "light" : "dark";
  document.documentElement.dataset.theme = next;
  localStorage.setItem("codon-theme", next);
  $("themeToggle").setAttribute("aria-label", next === "dark" ? "Switch to light theme" : "Switch to dark theme");
}

async function initRuntime() {
  try {
    setRuntimeStatus("Loading Python runtime…", "loading");
    pyodide = await loadPyodide({ indexURL: PYODIDE_INDEX });

    const response = await fetch("./codon_optimization/engine.py", { cache: "no-store" });
    if (!response.ok) {
      throw new Error("Unable to load the analysis engine (" + response.status + ")");
    }

    const source = await response.text();
    await pyodide.runPythonAsync(source);
    await pyodide.runPythonAsync(BRIDGE_SOURCE);

    setRuntimeStatus("Python ready", "ready");
    $("analyzeButton").disabled = false;
  } catch (error) {
    console.error(error);
    setRuntimeStatus("Runtime failed", "error");
    setError("Python runtime failed to load. Check the network connection and reload the page.");
  }
}

document.addEventListener("DOMContentLoaded", () => {
  initTheme();
  updateDnaMeta();

  $("themeToggle").addEventListener("click", toggleTheme);
  $("dnaSequence").addEventListener("input", updateDnaMeta);
  $("sampleButton").addEventListener("click", loadSample);
  $("clearButton").addEventListener("click", clearWorkspace);
  $("analyzeButton").addEventListener("click", analyze);
  $("copyButton").addEventListener("click", copyOptimized);
  $("downloadButton").addEventListener("click", downloadResult);

  document.querySelectorAll(".mobile-tab").forEach((button) => {
    button.addEventListener("click", () => setMobilePanel(button.dataset.panel));
  });

  initRuntime();
});
