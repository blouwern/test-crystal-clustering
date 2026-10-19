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

## NO.8 <2026-10-16>
- Recovered the missing Git LFS event data (`data/processed/*.root` were 133-byte pointers; `git-lfs` was not installed). Installed `git-lfs` 3.5.1 into a git-ignored `.tools/`, fetched `data/processed/**` from the public LFS endpoint, verified 1M events per dataset.
- Refined the sweep grid: `Eth_wave` densified to 16 points (0.1/0.2/0.3/0.5/0.75/1/1.5/2/2.5/3/4/5/6/8/10/12 MeV), `Nth` trimmed to 0/1/2, `Eth_seed` unchanged. 85 x 3 x 3 = 765 points, 0 failures.
- Added physics diagnostics to each point: median containment, unclustered energy, missed lit energy (energy above `Eth_wave` never absorbed), all-energy-in-one-cluster fraction.
- Rewrote `sweeps/analyze.py`: per-`Nth` heatmaps (rows = dataset + mean, columns = `Nth`) for 7 metrics, plus annotated zoom (0.1-1.5 MeV) and full-range contour overviews, wave-threshold profiles, trade-off scatter and best-parameter tables.
- Results: `Nth=0, seed=7, wave=0.1` is best (mean exact 0.8677, energy capture 0.9755). The optimum wave threshold moves **down to 0.1-0.2 MeV** from the previously recommended 0.5 MeV; `Eth_wave` is a monotonic over-splitting/energy-loss knob whose transition lies in 0.3-2 MeV (the old coarse grid missed it). `Nth=3` confirmed to be a large loss (-33 pp). Acceptance floor ~0.88 (12.5% of events deposit too little energy to seed).
- New deliverables under `pytrial/results/sweep/26-10-16/`: `sweep_summary.csv` (765 rows, 32 cols), `README.md`/`.html`, `best_params.md`/`.html`, 7 per-`Nth` heatmaps (+ full-range variants and per-`Nth` single figures), 3 wave profiles, trade-off scatter and best-point bars.

## NO.9 <2026-10-17>
- Follow-up on the two open items of the 26-10-16 report (§7 #2 and #3). New campaign stored in a separate folder `pytrial/results/sweep/26-10-16-refined/`.
- Grid: `Eth_seed` 5, 6, 7, 8, 9, 10, 12, 15, 20 (densified 5-10 MeV); `Eth_wave` **0.02, 0.05** + 0.1 ... 12 (16 points); `Nth` 0/1/2. 124 x 3 x 3 = **1116 points, 0 failures**.
- Made the grid env-overridable (`ECAL_SWEEP_SEEDS` / `ECAL_SWEEP_WAVES`), decoupled the deliverable folder (`ECAL_SWEEP_RESULTS_DATE`) from the working-area cache tag (`ECAL_SWEEP_TAG`), and added `sweeps/run_refined.sh` for one-command reproduction.
- Result #2: the sub-0.1 MeV plateau is confirmed — 0.02/0.05/0.1 MeV differ by <0.001 exact-match for `Nth=0/1` (0.0035 for `Nth=2`) and <0.0015 energy capture, with `frac_3plus` unchanged, so no noise-induced over-splitting.
- Result #3: the seed optimum is a broad 6-8 MeV plateau (peak at 7) for `Nth=0/1`; `Nth=2` still prefers the 5 MeV edge. The densification improved `mu_rest` from 0.8597 (seed=5, wave=0.1) to **0.8614 (seed=6, wave=0.02)**, +0.17 pp; mean optima are 0.8678 / 0.8641 / 0.8393 for `Nth=0/1/2`. The previous recommendation (`Nth=0/1, seed~7, wave=0.02-0.1`) is confirmed.
- Added a dedicated 0.02-0.5 MeV "low-threshold zoom" heatmap; fixed the `best_params.md` "mean" rows (they previously reported the single best dataset instead of the dataset average).
- Deliverables: `sweep_summary.csv` (1116 rows, 32 cols), `README.md`/`.html`, `best_params.md`/`.html`, per-`Nth` heatmaps (zoom / low-zoom / full-range) for 7 metrics, 3 wave profiles, trade-off scatter, best-point bars.

## NO.10 <2026-10-18>
- Wrote the consolidated algorithm performance report `pytrial/results/sweep/perf_report_1018.md` (667 lines) covering: data structures, the two-threshold wavefront algorithm and parameters, the cluster flow, all three campaigns (2301 points), metrics, results, failure modes, known issues, and the next-step plan.
- Key synthesis: the single-particle optimum is `Nth=0, Eth_seed=7 MeV, Eth_wave=0.02-0.1 MeV` (exact 0.8678, energy capture 0.9764, containment median 0.990, over-split 0.74%); the 12.48% zero-cluster rate is exactly the energy acceptance floor (`max(E) < Eth_seed`), and the single-particle metric is saturated (0.7 pp from the floor).
- Future plan in the report: (1) build a proper multi-particle stitcher (per-primary truth labels, overlap control, N=1-5, cluster-level efficiency/purity and merge-rate-vs-distance metrics) before judging algorithm changes; (2) barrier-aware propagation (energy-ratio alpha/beta), local-maximum/watershed seeding, split/merge/reassign post-processing, adaptive thresholds, geometry/shape terms, and set-based/C++ performance work.
- Fixed a correctness issue found during the synthesis: `scoring._read_cluster_counts` silently dropped trailing zero-cluster events, biasing low-acceptance points (312/1116 points had `N_events < 50000`; all `Nth=0` and Top-10 points were unaffected). The function now respects the true event count; archived CSVs are flagged in the report for optional re-scoring.

# Agent suggestion
<!-- agent may give the suggestion of instructions to be given -->
