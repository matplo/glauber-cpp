# Glauber C++ for Pb–Pb

A standalone C++17 Monte Carlo Glauber generator for minimum-bias Pb-208 collisions. Build and run it with CMake and the C++ standard library. The implementation follows the published model and the documented Pbpnrw parameters.

[Project documentation, Pb–Pb figures, and validation](https://matplo.github.io/glauber-cpp/)

## Run

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
build/glauber-pbpb pbpb-5020 60000 72913
```

The arguments are a new output directory, the number of impact-parameter proposals, and an unsigned random seed. Optional fourth and fifth arguments set the inelastic nucleon–nucleon cross section in mb and the centrality bin width in percent. The defaults are 67.6 mb at 5.02 TeV and 5%-wide bins. The bin width must divide 100. Repeat runs with the same executable and seed give the same sample.

- `events.csv` records inelastic events ordered by impact parameter, with proposal index, `b`, `N_part`, `N_coll`, event `T_AA`, and centrality class.
- `centrality.csv` records the boundaries and class means. Its `T_AA` standard errors hold the estimated boundaries fixed.
- `metadata.json` records the model settings, seed, trial counts and hadronic Pb–Pb cross section estimate. It is written after the CSV files close successfully.

Proposals sample uniformly in transverse area out to 20 fm. The generator accepts an event if at least one pair meets the inelastic hard-sphere criterion. It assigns percentile classes by ranking accepted events in `b`. Thus the classes represent impact-parameter centrality. `T_AA = N_coll / sigma_NN_inelastic` is in inverse mb. The first and last `b` boundaries are the proposal limits.

The model uses separate proton and neutron two-parameter Fermi profiles, radial reweighting, a 0.4 fm minimum separation, and center-of-mass recentering with a 0.1 fm shift limit. It uses continuous nucleon placement. The published improved benchmark also has a lattice option. [Glauber model reference](https://arxiv.org/abs/1710.07098)

## TGlauberMC comparison

Use ROOT only for this check. A separate Conda environment with `root_base` is sufficient. Put [TGlauberMC v3.2 `tglaubermc.hh` and `tglaubermc.cxx`](https://github.com/matplo/TGlauberMC) in a reference working directory alongside a copy of [the reference macro](validation/reference_fast.C). Run the macro in that directory. It uses Pbpnrw, continuous placement and the same cross section and proposal limits as the native run.

```sh
for i in 1 2 3 4; do
  root -l -b -q "reference_fast.C(15000,$((72920+i)),\"reference-shard-$i.csv\")"
done
```

From the source checkout, with paths adjusted to the reference working directory:

```sh
python3 validation/compare.py pbpb-5020 \
  --reference /path/to/reference/reference-shard-{1,2,3,4}.csv \
  --reference-trials-per-file 15000 --reference-seeds 72921 72922 72923 72924 \
  --reference-source /path/to/reference/tglaubermc.{hh,cxx} \
  --root-version "$(root-config --version)" --output comparison.json
```

The comparison reranks each sample independently and recomputes its percentile boundaries in 20 jackknife replicates. It checks the hadronic cross section and compares the `b` boundaries, mean `N_part` and mean `T_AA` against four combined Monte Carlo standard errors. This is a statistical comparison with independent seeds. Parameter and centrality-model uncertainties require separate studies.

For the [validation run](validation/pbpb-5020.json), each implementation drew 60,000 proposals and accepted 36,598 (C++) and 36,533 (TGlauberMC) inelastic events. All 20 bins passed. The largest absolute pulls were 1.03 for a boundary, 1.11 for mean N_part and 1.30 for mean T_AA. The inelastic cross sections were 7.665 ± 0.025 b and 7.651 ± 0.025 b, respectively. The [native CSV summary](validation/native-run/centrality.csv) and [metadata](validation/native-run/metadata.json) accompany the report.
