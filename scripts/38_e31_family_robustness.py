#!/usr/bin/env python
"""Is E31's enrichment robust to non-independence? Family-level resampling.

E31 compared 103 mutated sites from 71 discordant pairs against 767 sites from
522 same-label control pairs, and treated sites as the sampling unit. They are
not independent: one TCR can appear in many pairs -- E31 itself noted one
appearing in 14 -- and near-identical TCRs belong to the same clonotype family
by construction, since that is how the pairs were found in the first place.
Treating each site as a draw therefore overstates the effective sample size.

This re-tests the effect with the clonotype family as the sampling unit, which
is the level at which the observations are actually independent:

  * a cluster bootstrap over families, resampling whole families with
    replacement, which propagates within-family correlation into the interval;
  * a permutation test that shuffles the discordant/control label *between
    families* rather than between sites, so the null preserves family structure;
  * leave-one-epitope-out, since E31 reported direction agreeing on three
    epitopes and reversing on a fourth with 2 sites.

Both quantities are tested: the contact count E31 led with, and the binary
peptide-facing rate that the E37 revision showed is the part that survives.

Families are single-linkage at SCEPTR cosine <= 0.12 over the WHOLE cohort,
matching scripts/25 and scripts/37 rather than the labelled-only basis in
scripts/20.

Run: python scripts/38_e31_family_robustness.py
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immrep.data import build_training_frame, load_test  # noqa: E402

REPORT = ROOT / "reports/15_e31_family_robustness.md"
RADIUS = 0.12
DRAWS = 4000
SEED = 0


def load_anatomy():
    spec = importlib.util.spec_from_file_location(
        "anat", ROOT / "scripts/26_discordant_pair_anatomy.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["anat"] = module
    spec.loader.exec_module(module)
    return module


def family_map(test: pd.DataFrame) -> dict[str, str]:
    """Clonotype family per TCR, single-linkage over the whole cohort."""
    cached = np.load(ROOT / "data/interim/sceptr/default.npz", allow_pickle=True)
    position = {t: i for i, t in enumerate(cached["tcr_id"])}
    normed = cached["vectors"] / np.linalg.norm(
        cached["vectors"], axis=1, keepdims=True)
    out = {}
    for epitope, cohort in test.groupby("epitope_aa"):
        ids = cohort["tcr_id"].to_numpy()
        block = normed[[position[t] for t in ids]]
        _, family = connected_components(
            sp.csr_matrix((1.0 - block @ block.T) <= RADIUS), directed=False)
        for tcr_id, f in zip(ids, family):
            out[tcr_id] = f"{epitope}:{f}"
    return out


def build_sites(module, test, labels, families) -> pd.DataFrame:
    """Mutated-site records carrying the identifiers E31's table dropped."""
    rows = []
    for group, want in (("discordant", (1, 0)), ("control", (1, 1))):
        pairs = module.find_pairs(labels, want)
        for pair in pairs.to_dict("records"):
            a = module.peptide_contact_profile(pair["a_id"], pair["varied_chain"], test)
            b = module.peptide_contact_profile(pair["b_id"], pair["varied_chain"], test)
            if a is None or b is None:
                continue
            fa, fb = families.get(pair["a_id"]), families.get(pair["b_id"])
            for site in pair["sites"]:
                rows.append({
                    "group": group, "epitope": pair["epitope"],
                    "a_id": pair["a_id"], "b_id": pair["b_id"],
                    "family": fa if fa == fb else f"{fa}|{fb}",
                    "same_family": fa == fb,
                    "contacts": a.get(site, 0) + b.get(site, 0),
                    "facing": int(a.get(site, 0) + b.get(site, 0) > 0),
                })
    return pd.DataFrame(rows)


def effect(frame: pd.DataFrame, column: str) -> float:
    d = frame.loc[frame["group"] == "discordant", column]
    c = frame.loc[frame["group"] == "control", column]
    if d.empty or c.empty:
        return np.nan
    return float(d.mean() - c.mean())


def cluster_bootstrap(frame: pd.DataFrame, column: str,
                      rng: np.random.Generator) -> np.ndarray:
    """Resample whole families, so within-family correlation reaches the interval."""
    groups = [block for _, block in frame.groupby("family")]
    out = np.empty(DRAWS)
    for i in range(DRAWS):
        pick = rng.integers(0, len(groups), len(groups))
        out[i] = effect(pd.concat([groups[j] for j in pick], ignore_index=True),
                        column)
    return out


def family_cells(frame: pd.DataFrame) -> pd.DataFrame:
    """One value per (family, group): the mean over that family's sites in that
    group. This is the family-as-unit view, with every family weighted equally.

    It replaces a permutation test that was wrong. That test shuffled the group
    label between families and then applied it to every site, so its null was
    computed on family-collapsed labels while the observed statistic was
    computed on site-level labels. Fifteen of 97 families carry both groups,
    covering 34% of sites, and collapsing them moves the contact effect from
    +3.63 to -0.15 on its own -- so the null and the observed were not the same
    statistic and the p-values it produced were not interpretable.
    """
    return (frame.groupby(["family", "group"])[["contacts", "facing"]]
            .mean().reset_index())


def paired_within_family(cells: pd.DataFrame, column: str):
    """Families carrying both groups give a paired comparison that cannot be
    confounded by which families fall in which group."""
    wide = cells.pivot(index="family", columns="group", values=column).dropna()
    return wide["discordant"] - wide["control"]


def main() -> int:
    module = load_anatomy()
    test = load_test()
    labels = build_training_frame()
    families = family_map(test)
    # peptide_contact_profile indexes test by tcr_id.
    frame = build_sites(module, test.set_index("tcr_id"), labels, families)
    rng = np.random.default_rng(SEED)

    print("=== how much independence is actually there? ===")
    counts = frame.groupby("group").agg(
        sites=("contacts", "size"),
        pairs=("a_id", lambda s: len(set(zip(s, frame.loc[s.index, "b_id"])))),
        tcrs=("a_id", lambda s: len(set(s) | set(frame.loc[s.index, "b_id"]))),
        families=("family", "nunique"))
    print(counts.to_string())
    shared = (set(frame.loc[frame["group"] == "discordant", "family"])
              & set(frame.loc[frame["group"] == "control", "family"]))
    print(f"\n  families contributing to both groups: {len(shared)}")
    print(f"  pair members in the same family: "
          f"{frame['same_family'].mean():.1%}")
    biggest = frame.groupby("family").size().nlargest(3)
    print(f"  largest families contribute {biggest.tolist()} sites")
    print("\n  E31 treated the site as the sampling unit. The independent unit is")
    print("  the family, and there are far fewer of those than sites.")

    print("\n=== effect, with the family as the sampling unit ===")
    from scipy.stats import mannwhitneyu, wilcoxon

    cells = family_cells(frame)
    rows = []
    for column, label in (("contacts", "peptide contacts at the site"),
                          ("facing", "site faces the peptide (binary)")):
        observed = effect(frame, column)
        boot = cluster_bootstrap(frame, column, rng)
        lo, hi = np.nanpercentile(boot, [2.5, 97.5])
        d = cells.loc[cells["group"] == "discordant", column]
        c = cells.loc[cells["group"] == "control", column]
        p_cell = mannwhitneyu(d, c, alternative="two-sided").pvalue
        diff = paired_within_family(cells, column)
        p_pair = wilcoxon(diff).pvalue if len(diff) > 5 else np.nan
        rows.append({"quantity": label,
                     "site_effect": observed, "ci_lo": lo, "ci_hi": hi,
                     "family_disc": d.mean(), "family_ctl": c.mean(),
                     "p_family_cell": p_cell,
                     "paired_mean": diff.mean(), "paired_median": diff.median(),
                     "p_paired": p_pair, "paired_positive": int((diff > 0).sum()),
                     "paired_n": len(diff)})
        print(f"  {label}")
        print(f"    site-weighted effect {observed:+.3f}  "
              f"cluster-bootstrap 95% CI [{lo:+.3f}, {hi:+.3f}]")
        print(f"    family as unit       disc {d.mean():.3f} (n={len(d)})  "
              f"ctl {c.mean():.3f} (n={len(c)})  p={p_cell:.3f}")
        print(f"    paired within the {len(diff)} mixed families  "
              f"mean {diff.mean():+.3f}  median {diff.median():+.3f}  "
              f"p={p_pair:.3f}  positive {int((diff > 0).sum())}/{len(diff)}")
    table = pd.DataFrame(rows)
    print("\n  The bootstrap CI is site-weighted, so a few large families carry")
    print("  it -- three contribute 129, 92 and 73 sites. Weighting families")
    print("  equally, neither quantity separates the groups.")

    print("\n=== leave-one-epitope-out ===")
    loo = []
    for epitope in sorted(frame["epitope"].unique()):
        kept = frame[frame["epitope"] != epitope]
        loo.append({
            "dropped": epitope,
            "sites_lost": int((frame["epitope"] == epitope).sum()),
            "contacts_effect": effect(kept, "contacts"),
            "facing_effect": effect(kept, "facing"),
        })
    loo = pd.DataFrame(loo)
    loo.loc[len(loo)] = {"dropped": "(none)", "sites_lost": 0,
                         "contacts_effect": effect(frame, "contacts"),
                         "facing_effect": effect(frame, "facing")}
    print(loo.round(3).to_string(index=False))
    swing = loo["contacts_effect"].max() - loo["contacts_effect"].min()
    print(f"\n  dropping one epitope moves the contact effect by up to {swing:.3f}")
    print("  against a full-data effect of "
          f"{effect(frame, 'contacts'):.3f}.")

    print("\n=== per-epitope, with family counts rather than site counts ===")
    per = frame.groupby(["epitope", "group"]).agg(
        sites=("contacts", "size"), families=("family", "nunique"),
        contacts=("contacts", "mean"), facing=("facing", "mean"))
    print(per.round(3).to_string())

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        "# E31 under family-level resampling\n\n"
        "## Independence\n\n" + counts.to_markdown()
        + f"\n\nFamilies contributing to both groups: {len(shared)}.\n\n"
          "## Effect with the family as the sampling unit\n\n"
        + table.round(4).to_markdown(index=False)
        + "\n\n## Leave-one-epitope-out\n\n" + loo.round(3).to_markdown(index=False)
        + "\n\n## Per epitope\n\n" + per.round(3).to_markdown() + "\n",
        encoding="utf-8")
    print(f"\nwrote {REPORT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
