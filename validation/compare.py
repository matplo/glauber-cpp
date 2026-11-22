"""Compare independent ROOT-free and TGlauberMC PbPb event samples."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import mean


def read_events(path, trials, shard=0):
    events = []
    seen = set()
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            trial, b, npart, ncoll = int(row["trial"]), float(row["b_fm"]), int(row["npart"]), int(row["ncoll"])
            if not (trial not in seen and 1 <= trial <= trials and
                    0 <= b <= 20 and 2 <= npart <= 416 and 1 <= ncoll <= 208**2 and npart <= 2*ncoll):
                raise ValueError(f"Invalid inelastic event in {path}: {row}")
            seen.add(trial)
            events.append((b, npart, ncoll, (trial + shard) % 20))
    return events


def measurements(events, bins, sigma):
    n = len(events)
    if n < 2*bins:
        raise ValueError("Too few accepted events")
    values = []
    for index in range(bins):
        first, last = n*index//bins, n*(index+1)//bins
        sample = events[first:last]
        boundary = (events[last-1][0] + events[last][0])/2 if last < n else 20
        values.append((boundary, mean(e[1] for e in sample), mean(e[2] for e in sample)/sigma))
    return values


def summarize(events, bins, sigma):
    events = sorted(events)
    result = measurements(events, bins, sigma)
    resampled = [measurements([e for e in events if e[3] != block], bins, sigma) for block in range(20)]
    errors = []
    for samples in zip(*resampled):
        errors.append(tuple(math.sqrt(19/20*sum((v-mean(values))**2 for v in values)) for values in zip(*samples)))
    return result, errors


def check_native_summary(directory, measured, events):
    width = 100//len(measured)
    with (directory/"centrality.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != len(measured):
        raise ValueError("Unexpected number of native bins")
    for index, (row, values) in enumerate(zip(rows, measured)):
        expected = len(events)*(index+1)//len(rows) - len(events)*index//len(rows)
        actual = [float(row[k]) for k in ("b_high_fm", "npart_mean", "taa_mean_mb_inverse")]
        if (int(row["centrality_low_pct"]) != index*width or int(row["centrality_high_pct"]) != (index+1)*width
            or int(row["events"]) != expected
            or not all(math.isclose(a, b, rel_tol=1e-10) for a, b in zip(actual, values))):
            raise ValueError(f"Native CSV is inconsistent with event ranks: {row}")


def comparison(native, reference, native_error, reference_error):
    error = math.hypot(native_error, reference_error)
    pull = (native-reference)/error if error else None
    return dict(native=native, reference=reference, native_se=native_error, reference_se=reference_error,
                pull=pull, passed=abs(pull) <= 4 if pull is not None else native == reference)


def cross_section(count, trials):
    fraction = count/trials
    area = 4*math.pi
    return area*fraction, area*math.sqrt(fraction*(1-fraction)/trials)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("native", type=Path)
    p.add_argument("--reference", type=Path, nargs="+", required=True)
    p.add_argument("--reference-trials-per-file", type=int, required=True)
    p.add_argument("--reference-seeds", type=int, nargs="+", required=True)
    p.add_argument("--reference-source", type=Path, nargs="+", required=True)
    p.add_argument("--root-version", required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if len(a.reference_seeds) != len(a.reference) or len(set(a.reference_seeds)) != len(a.reference):
        p.error("supply one distinct seed per reference file")
    metadata = json.loads((a.native/"metadata.json").read_text())
    sigma = metadata["sigma_nn_inelastic_mb"]
    if sigma != 67.6 or metadata["max_impact_parameter_fm"] != 20:
        p.error("the supplied TGlauberMC macro uses 67.6 mb and bmax=20 fm")
    if metadata["seed"] in a.reference_seeds:
        p.error("use independent reference seeds")
    width = metadata["centrality_bin_width_pct"]
    native = read_events(a.native/"events.csv", metadata["trials"])
    reference = [event for i,path in enumerate(a.reference) for event in read_events(path, a.reference_trials_per_file, i)]
    if len(native) != metadata["inelastic_events"]:
        raise ValueError("Native event count differs from metadata")
    nval, nerr = summarize(native, 100//width, sigma)
    rval, rerr = summarize(reference, 100//width, sigma)
    check_native_summary(a.native, nval, native)
    rows = []
    for i, (n, r, ne, re) in enumerate(zip(nval, rval, nerr, rerr)):
        rows.append(dict(centrality_percent=[i*width,(i+1)*width], **{
            key: comparison(*values)
            for key, values in zip(("b_high_fm","npart_mean","taa_mean_mb_inverse"), zip(n,r,ne,re))
        }))
    reference_trials = len(a.reference)*a.reference_trials_per_file
    ns, nse = cross_section(len(native), metadata["trials"])
    rs, rse = cross_section(len(reference), reference_trials)
    xsec = comparison(ns, rs, nse, rse)
    report = dict(
        native=dict(metadata, events_sha256=hashlib.sha256((a.native/"events.csv").read_bytes()).hexdigest()),
        reference=dict(root_version=a.root_version, trials=reference_trials, inelastic_events=len(reference),
                       seeds=a.reference_seeds,
                       events_sha256={f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in a.reference},
                       source_sha256={f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                                                               for f in a.reference_source}),
        method="Independent inelastic impact-parameter quantiles. Recompute bin boundaries and means "
               "in each of 20 delete-one-block jackknife replicates. Compare with combined Monte Carlo "
               "standard errors and a four-standard-error threshold. The 100% endpoint is fixed by the proposal.",
        sigma_pbpb_inelastic_barn=xsec,
        bins=rows,
        passed=xsec["passed"] and all(value["passed"] for row in rows for key,value in row.items()
                                     if key != "centrality_percent"),
    )
    with a.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"TGlauberMC validation {'passed' if report['passed'] else 'failed'} for {len(rows)} centrality bins")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
