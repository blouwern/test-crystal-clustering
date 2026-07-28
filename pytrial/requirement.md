# Original user requirement
<!-- user instrucion log -->

## NO.1 <2026-07-27 10:40>
- Read and understand the whole program in `/test-crystal-clustering`, note any `note.md` for quick information.
- Complete the `scoring.py` numerical analyse

## NO.2 <2026-07-27 11:14>
- Check line 167 to line 177 of file `reconstruction.py` refine it to make the default parameter setting work.

## NO.3 <2026-07-28 17:37>
- Fix `scoring.py` ValueError: chisquare sum of observed vs expected mismatch — build full range [0..max_val] + tail bin, merge low-expected bins instead of masking, so probability mass is fully accounted.

## NO.4 <2026-07-28 17:58>
- Parameterize `scoring.py`: add `expected_lambda` param (None → auto-estimate MLE from data with ddof=1, specified → use value with ddof=0). Make inline "expected value" label and "% exactly X clusters" adapt to the lambda.
- Improve `reconstruction.py` file naming: add `output_path` param (None → auto-name, specified → use value). Wire CLI to accept optional output filename and expected_lambda.

## NO.5 <2026-07-28 18:17>
- Add `NTrack` branch to `evt_stitch_child.py` output ROOT file to store ground-truth track count per event.
- Refine `scoring.py`: add `event_by_event` mode (per-event cluster vs truth comparison) and empirical distribution from data file (two-sample KS, chi-squared vs ground-truth). Add `data_path` parameter to read `NTrack` branch.
- Create `sweeps/` directory: `config.py` (parameter grid), `run_sweep.py` (automated sweeping + CSV summary), `analyze.py` (heatmaps + best-param printout), `results/` (ROOT outputs), `analysis_output/` (CSV + plots).
- `Nth` has negligible effect; `Eth=10` is optimal (42.9% exact match, MAE=0.74, r=0.84). `Eth=5` over-clusters, `Eth=20` under-clusters, `Eth=50` misses nearly everything.

## NO.6 <2026-07-28 18:28>
- Fine sweep: Eth ∈ [5,6,7,8,9,10,11,12,13,14,15,20], Nth ∈ [1,2,3] (36 combinations).
- **Best exact match: Eth=11, Nth=2 (43.27%)** — marginally beats Eth=10,Nth=1 (43.17%).
- **Best MAE/RMSE: Eth=10, Nth=1** (MAE=0.7449, RMSE=1.0859).
- **Best mean proximity to truth: Eth=9, Nth=1** (mean=2.95 vs truth 3.01).
- **Overall optimum: Eth=9–10, Nth=1** — balances all metrics. Nth=1 favored for Eth≤9; Nth=2 slightly better for Eth≥11 but difference is <0.5%.
- Correlation peaks at Eth=6 (r=0.856) due to consistent over-prediction (+0.7 bias); MAE degrades at lower Eth due to noise seeding.

## NO.7 <2026-07-28 18:51>
- Add `.gitignore` for Python artifacts, ROOT outputs, plots, sweep results.
- Create `pipeline.md` documenting both `reconstruction.py` (algorithm, API, CLI, output schema, data flow) and `scoring.py` (modes, API, examples); renamed from `scoring.md`.

# Agent suggestion
<!-- agent may give the suggestion of instructions to be given -->
