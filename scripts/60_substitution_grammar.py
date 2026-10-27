#!/usr/bin/env python
"""E62 -- what explains how a single substitution changes a peptide's response?

The roadmap asks whether a position-conditioned substitution model beats fixed
matrices (BLOSUM62, Grantham, physicochemical descriptors) at explaining
response change, under peptide-family holdout rather than a random split.

Substrate, chosen after measuring the alternatives. Single mutants *relative to
the wild type* are too few: the Flu-M1 mutagenesis table holds 69 of them
covering 20 of 171 possible (position, residue) cells, with positions 2, 4 and 6
entirely absent, and it is a top-10,000 table so substitutions that killed the
response were never recorded. The roadmap's own fallback applies. What it also
says, and what works, is *peptide pairs within one screen that differ at one
residue* -- not necessarily from the wild type:

    Flu(4TCR)  7,122 pairs   EBV(4TCR) 3,790   1G4(4TCR) 1,587
    Flu-mut    1,221         DMF5-9mer 1,191   DMF5(4TCR)  247
    TCR59         15         TCR63         1   <- B27 screens unusable here

**The design that matters is the 4TCR block**: four receptors screened against
the same NNK 9-mer library. The same substitution at the same position is then
measured on four different TCRs, which separates the two explanations the
roadmap says must not be conflated -- an effect shared across all four is
library or presentation, an effect that differs between them is TCR-specific
recognition. Nothing in a single screen can make that split.

Two things held to deliberately:

  * the response is the **later-round gain** (`r3_r2_ratio`), because E61
    measured that first-round abundance runs *negative* against the functional
    readout while the later gain runs +0.687;
  * holdout is by **connected component of the Hamming-1 graph**, so a pair's
    own neighbours cannot sit in the training fold.

Run: .venv/bin/python scripts/60_substitution_grammar.py
"""

from __future__ import annotations

import glob
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports/35_substitution_grammar.md"
SEED = 0

# Kyte-Doolittle hydropathy, side-chain volume (A^3), charge.
HYDRO = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
         "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
         "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
         "Y": -1.3, "V": 4.2}
VOLUME = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5,
          "Q": 143.8, "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7,
          "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7,
          "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0}
CHARGE = {**{a: 1.0 for a in "KRH"}, **{a: -1.0 for a in "DE"}}


def blosum62():
    from Bio.Align import substitution_matrices
    return substitution_matrices.load("BLOSUM62")


def grantham_table():
    """Grantham 1974, from the composition/polarity/volume formula."""
    c = {"A": 0.0, "R": 0.65, "N": 1.33, "D": 1.38, "C": 2.75, "Q": 0.89,
         "E": 0.92, "G": 0.74, "H": 0.58, "I": 0.0, "L": 0.0, "K": 0.33,
         "M": 0.0, "F": 0.0, "P": 0.39, "S": 1.42, "T": 0.71, "W": 0.13,
         "Y": 0.20, "V": 0.0}
    pol = {"A": 8.1, "R": 10.5, "N": 11.6, "D": 13.0, "C": 5.5, "Q": 10.5,
           "E": 12.3, "G": 9.0, "H": 10.4, "I": 5.2, "L": 4.9, "K": 11.3,
           "M": 5.7, "F": 5.2, "P": 8.0, "S": 9.2, "T": 8.6, "W": 5.4,
           "Y": 6.2, "V": 5.9}
    vol = {"A": 31, "R": 124, "N": 56, "D": 54, "C": 55, "Q": 85, "E": 83,
           "G": 3, "H": 96, "I": 111, "L": 111, "K": 119, "M": 105, "F": 132,
           "P": 32.5, "S": 32, "T": 61, "W": 170, "Y": 136, "V": 84}
    a1, b1, g1 = 1.833, 0.1018, 0.000399
    out = {}
    for x in c:
        for y in c:
            out[(x, y)] = float(np.sqrt(
                a1 * (c[x] - c[y]) ** 2 + b1 * (pol[x] - pol[y]) ** 2
                + g1 * (vol[x] - vol[y]) ** 2))
    return out


def hamming_pairs(sequences, response):
    """All pairs differing at exactly one position, via position-masked keys."""
    length = len(sequences[0])
    rows = []
    for p in range(length):
        buckets = defaultdict(list)
        for s in sequences:
            buckets[s[:p] + "." + s[p + 1:]].append(s)
        for members in buckets.values():
            if len(members) < 2:
                continue
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    x, y = members[i], members[j]
                    rows.append({"a_seq": x, "b_seq": y, "pos": p + 1,
                                 "a": x[p], "b": y[p],
                                 "delta": abs(response[x] - response[y])})
    return pd.DataFrame(rows)


def families(frame, sequences):
    """Connected components of the Hamming-1 graph, so a pair's neighbours
    cannot land in the training fold."""
    index = {s: i for i, s in enumerate(sequences)}
    rows = [index[s] for s in frame["a_seq"]]
    cols = [index[s] for s in frame["b_seq"]]
    graph = sp.csr_matrix((np.ones(len(rows)), (rows, cols)),
                          shape=(len(sequences), len(sequences)))
    _, labels = connected_components(graph, directed=False)
    return frame["a_seq"].map(lambda s: labels[index[s]])


def main() -> int:
    bl = blosum62()
    gr = grantham_table()

    screens = {}
    for path in sorted(glob.glob(
            str(ROOT / "data/external/seeker/4TCR_screen/top10000_*.csv"))):
        name = path.split("4TCR-")[1].split("-SCD")[0]
        screens[name] = pd.read_csv(path)
    print(f"=== the 4TCR block: {len(screens)} receptors on one NNK library ===")
    print("  Response is r3_r2_ratio, the later-round gain, because E61 found")
    print("  first-round abundance runs negative against the flow readout.\n")

    pairs = {}
    for name, frame in screens.items():
        frame = frame[frame["aa"].astype(str).str.isalpha()].copy()
        length = frame["aa"].str.len().mode().iloc[0]
        frame = frame[frame["aa"].str.len() == length]
        # These tables are ranked per NT CONSTRUCT, not per peptide: 10,000
        # rows hold 8,545 distinct 9-mers and one peptide can carry up to 16
        # synonymous encodings. Treating rows as peptides manufactures
        # duplicate "pairs" with a delta of exactly zero -- a first run showed
        # median |delta| 0.000 and 42,810 Flu pairs against the 7,122 a
        # deduplicated count gives. Reads are summed over a peptide's
        # encodings, which is the unit the biology is about.
        counts = (frame.assign(**{c: frame[c].fillna(0.0) for c in
                                  ("count_round1", "count_round2",
                                   "count_round3")})
                  .groupby("aa")[["count_round1", "count_round2",
                                  "count_round3"]].sum())
        totals = counts.sum()
        # A peptide can be unobserved in a round (NaN -> 0 reads), so the ratio
        # needs a declared pseudocount rather than silently becoming NaN.
        alpha = 1.0
        pct2 = (counts["count_round2"] + alpha) / (totals["count_round2"] + alpha)
        pct3 = (counts["count_round3"] + alpha) / (totals["count_round3"] + alpha)
        response = dict(np.log10(pct3 / pct2))
        sequences = list(counts.index)
        block = hamming_pairs(sequences, response)
        block["family"] = families(block, sequences)
        block["blosum"] = [-float(bl[x, y]) for x, y in zip(block["a"], block["b"])]
        block["grantham"] = [gr[(x, y)] for x, y in zip(block["a"], block["b"])]
        block["d_hydro"] = [abs(HYDRO[x] - HYDRO[y])
                            for x, y in zip(block["a"], block["b"])]
        block["d_volume"] = [abs(VOLUME[x] - VOLUME[y])
                             for x, y in zip(block["a"], block["b"])]
        block["d_charge"] = [abs(CHARGE.get(x, 0) - CHARGE.get(y, 0))
                             for x, y in zip(block["a"], block["b"])]
        pairs[name] = block
        print(f"  {name:<6} {len(counts):>6,} peptides, {len(block):>6,} "
              f"pairs, {block['family'].nunique():>4} families, "
              f"median |delta| {block['delta'].median():.3f}")

    print("\n=== Part 1: fixed matrices against a position-conditioned model ===")
    print("  Group-held-out by Hamming-1 connected component, 5 folds.\n")
    results = []
    for name, block in pairs.items():
        y = block["delta"].to_numpy()
        groups = block["family"].to_numpy()
        designs = {
            "BLOSUM62": block[["blosum"]].to_numpy(float),
            "Grantham": block[["grantham"]].to_numpy(float),
            "physicochemical": block[["d_hydro", "d_volume", "d_charge"]]
                               .to_numpy(float),
            "position only": pd.get_dummies(block["pos"]).to_numpy(float),
            "position + Grantham": np.hstack([
                pd.get_dummies(block["pos"]).to_numpy(float),
                block[["grantham"]].to_numpy(float)]),
            "position x chemistry": np.hstack([
                pd.get_dummies(block["pos"]).to_numpy(float),
                np.column_stack([
                    pd.get_dummies(block["pos"]).to_numpy(float)
                    * block[c].to_numpy(float)[:, None]
                    for c in ("d_hydro", "d_volume", "d_charge")])]),
        }
        for label, x in designs.items():
            out = np.empty(len(y))
            splitter = GroupKFold(n_splits=min(5, len(np.unique(groups))))
            for train, test in splitter.split(x, y, groups):
                out[test] = Ridge(alpha=1.0).fit(x[train], y[train]).predict(x[test])
            rho = spearmanr(out, y).statistic
            results.append({"screen": name, "model": label,
                            "held_out_spearman": float(rho), "n_pairs": len(y)})
        best = max((r for r in results if r["screen"] == name),
                   key=lambda r: r["held_out_spearman"])
        print(f"  {name:<6} best: {best['model']:<22} "
              f"rho {best['held_out_spearman']:+.3f}   " + "  ".join(
                  f"{r['model'][:4]} {r['held_out_spearman']:+.3f}"
                  for r in results if r["screen"] == name))
    table = pd.DataFrame(results)

    print("\n=== Part 2: the negative control ===")
    print("  Hold the position fixed, shuffle which chemical substitution it")
    print("  was. Any model that still works was reading position alone.\n")
    rng = np.random.default_rng(SEED)
    control = []
    for name, block in pairs.items():
        shuffled = block.copy()
        for _, idx in block.groupby("pos").groups.items():
            order = rng.permutation(len(idx))
            for column in ("grantham", "d_hydro", "d_volume", "d_charge"):
                shuffled.loc[idx, column] = block.loc[idx, column].to_numpy()[order]
        y = shuffled["delta"].to_numpy()
        groups = shuffled["family"].to_numpy()
        x = np.hstack([pd.get_dummies(shuffled["pos"]).to_numpy(float),
                       shuffled[["grantham"]].to_numpy(float)])
        out = np.empty(len(y))
        for train, test in GroupKFold(
                n_splits=min(5, len(np.unique(groups)))).split(x, y, groups):
            out[test] = Ridge(alpha=1.0).fit(x[train], y[train]).predict(x[test])
        real = table[(table["screen"] == name)
                     & (table["model"] == "position + Grantham")
                     ]["held_out_spearman"].iloc[0]
        control.append({"screen": name, "real": real,
                        "chemistry_shuffled": float(spearmanr(out, y).statistic)})
        print(f"  {name:<6} position+Grantham {real:+.3f}  ->  shuffled "
              f"{control[-1]['chemistry_shuffled']:+.3f}")
    control = pd.DataFrame(control)

    print("\n=== Part 3: shared across receptors, or TCR-specific? ===")
    print("  The same (position, a, b) cell measured on four receptors. A")
    print("  shared effect is library or presentation; a divergent one is")
    print("  recognition. Only the 4TCR block can make this split.\n")
    cells = {}
    for name, block in pairs.items():
        key = block.assign(cell=block["pos"].astype(str) + ":"
                           + np.minimum(block["a"], block["b"])
                           + np.maximum(block["a"], block["b"]))
        cells[name] = key.groupby("cell")["delta"].mean()
    common = pd.DataFrame(cells).dropna()
    print(f"  cells measured on all {len(cells)} receptors: {len(common)}")
    if len(common) >= 20:
        corr = common.corr(method="spearman")
        print(corr.round(3).to_string())
        off = corr.to_numpy()[np.triu_indices(len(corr), 1)]
        print(f"\n  mean cross-receptor rho {off.mean():+.3f}")

        # A low cross-receptor correlation only means "receptor-specific" if a
        # receptor agrees with ITSELF. Split each receptor's pairs for a cell
        # into halves and correlate the halves: that is the noise ceiling, and
        # without it this comparison cannot be read at all.
        print("\n  Noise ceiling -- does a receptor agree with itself? Split")
        print("  each cell's pairs in half, 50 repeats:\n")
        ceiling = {}
        for name, block in pairs.items():
            key = block.assign(cell=block["pos"].astype(str) + ":"
                               + np.minimum(block["a"], block["b"])
                               + np.maximum(block["a"], block["b"]))
            key = key[key["cell"].isin(common.index)]
            draws = []
            for _ in range(50):
                shuffled = key.sample(frac=1.0, random_state=int(
                    rng.integers(1 << 30)))
                half = shuffled.groupby("cell").cumcount() % 2
                a = shuffled[half == 0].groupby("cell")["delta"].mean()
                b = shuffled[half == 1].groupby("cell")["delta"].mean()
                both = pd.concat([a, b], axis=1).dropna()
                if len(both) >= 20:
                    draws.append(spearmanr(both.iloc[:, 0],
                                           both.iloc[:, 1]).statistic)
            # Spearman-Brown: a split-half uses half the pairs per cell.
            half = float(np.mean(draws)) if draws else np.nan
            ceiling[name] = (2 * half / (1 + half)
                             if np.isfinite(half) and half > -1 else np.nan)
            print(f"    {name:<6} split-half rho {ceiling[name]:+.3f} "
                  f"({len(draws)} usable repeats)")
        mean_ceiling = float(np.nanmean(list(ceiling.values())))
        # A cross-receptor correlation is attenuated by BOTH receptors' noise,
        # so comparing the raw number with the ceiling understates sharing.
        # Disattenuate per pair, and DROP any receptor whose own reliability is
        # too low to divide by -- a first version concluded "nothing can be
        # concluded" from the raw comparison alone.
        print("\n  Disattenuated per pair, r_obs / sqrt(rel_a * rel_b):")
        usable = [n for n in ceiling if ceiling[n] >= 0.3]
        dropped = [n for n in ceiling if ceiling[n] < 0.3]
        corrected = []
        for i, a in enumerate(corr.index):
            for b in list(corr.index)[i + 1:]:
                value = corr.loc[a, b] / np.sqrt(ceiling[a] * ceiling[b])
                flag = ("  <- dropped, reliability "
                        f"{min(ceiling[a], ceiling[b]):.2f}"
                        if a in dropped or b in dropped else "")
                print(f"    {a}-{b:<6} observed {corr.loc[a, b]:+.3f}  "
                      f"disattenuated {value:+.3f}{flag}")
                if a in usable and b in usable:
                    corrected.append(value)
        if dropped:
            print(f"    dropped for low self-reliability: {dropped}")
        if corrected:
            shared = float(np.mean(corrected))
            print(f"\n  Among the {len(usable)} usable receptors, the shared "
                  f"component is {shared:.0%};")
            print(f"  about {1 - shared:.0%} of a substitution's effect is "
                  "receptor-specific.")
            print("  Wide uncertainty: dividing by a reliability near 0.5 "
                  "roughly")
            print("  doubles the estimate and its error bars with it.")
        print(f"\n    mean within-receptor reliability {mean_ceiling:+.3f}")
        print(f"    mean cross-receptor agreement      {off.mean():+.3f}")
        if np.isnan(mean_ceiling) or mean_ceiling < 0.3:
            print("\n  The ceiling is as low as the cross-receptor number, so")
            print("  a receptor does not agree with itself either. NOTHING can")
            print("  be concluded about receptor-specificity from this -- the")
            print("  per-cell estimates are noise-dominated.")
        elif off.mean() < 0.5 * mean_ceiling:
            print("\n  Cross-receptor agreement is well below the ceiling, so")
            print("  the effect size IS receptor-specific rather than a")
            print("  property of the library or of presentation.")
        else:
            print("\n  Cross-receptor agreement approaches the ceiling: the")
            print("  effect is largely shared, i.e. library or presentation.")
    else:
        corr = pd.DataFrame()
        print("  too few shared cells -- NOT ESTIMABLE")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E62 - substitution grammar on the 4TCR NNK block\n\n"
        "Response is `r3_r2_ratio` (later-round gain, per E61). Holdout is by "
        "connected component of the Hamming-1 graph.\n\n"
        "## Held-out Spearman by model\n\n"
        + table.pivot(index="model", columns="screen",
                      values="held_out_spearman").round(4).to_markdown()
        + "\n\n## Negative control (chemistry shuffled within position)\n\n"
        + control.round(4).to_markdown(index=False)
        + "\n\n## Cross-receptor agreement of substitution effects\n\n"
        + (corr.round(4).to_markdown() if len(corr) else "Not estimable.")
        + "\n", encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
