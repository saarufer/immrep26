#!/usr/bin/env python
"""Short MD on discordant pairs: is the interface difference dynamically stable?

Plan card C11, which I twice called infeasible and twice for the wrong reason.
The first diagnosis blamed installing OpenMM with --no-deps; the actual cause is
that the PyPI openmm wheel ships no CUDA plugin at all, so only conda-forge
provides one.

What this can address. E31 found that discordant negatives mutate at positions
carrying more peptide contact than same-label controls, and E37 showed that
survives energy minimisation. Both are static, and both read a *contact count*,
which a static audit here showed is driven by a few large values: on the relaxed
structures the count ranks discordant sites above controls (5.9 vs 3.1, matching
E37's 7.16 vs 3.53) while the minimum site-to-peptide distance ranks them the
other way (6.10 A vs 5.39 A). Those two readings of the same quantity disagree
because 60% of sites have no contact at all at 4.5 A, so the count is zero for
most sites and large for a few. This run therefore reports both, and reports
them per frame, so persistence is measured rather than a single snapshot.

What it cannot address, per the plan's own instruction: this is not an affinity
calculation, trajectory length is not evidence of binding, and nothing here
should be read as k_on or k_off. Replica-to-replica spread is reported alongside
every number precisely so that a difference smaller than the spread cannot be
mistaken for a result. Two nanoseconds measures short-timescale persistence of a
relaxed pose, not converged equilibrium.

Design note, and the reason this script simulates complexes rather than pairs.
The first version resolved each pair's mutated CDR3 offsets to residue numbers
and passed them into the simulation. It did so by taking the *union* of the two
members' residue numbers, which is wrong: the two TCRs in a pair have different
V genes, so the CDR3 begins at a different position in each domain -- 88 in one
member of the first pair against 90 in the other -- and the union probed each
structure at two residues outside its own mutated site. Three of sixteen pairs
were affected and all three were discordant, so the contamination fell entirely
on the group where the effect is claimed. Its smoke test showed a clean 6.2
against 0.0 site contacts in the expected direction; with the correct residues
both members read 0.0, and the whole signal had come from the contaminating
residues. So the mutated site no longer enters the simulation at all. Each
unique complex is simulated once per replica and every CDR3 residue of both
chains is measured; which residues were mutated is applied afterwards, as a
lookup. That also means a TCR appearing in two pairs is simulated once.

Checked and not a problem: OpenMM renames chains on write, M and P becoming C
and D, so the peptide is identified by residue count rather than chain letter,
and the element column is written, so added hydrogens are correctly excluded.

Run: micromamba run -p ./.mmenv-openmm python scripts/34_md_interface_stability.py
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

RELAXED = ROOT / "data/interim/relaxed"
PAIRS = ROOT / "data/interim/md_pairs.json"
OUT = ROOT / "data/interim/md_residue_contacts.csv"
TRACE = ROOT / "data/interim/md_interface_trace.csv"

N_PAIRS = 20                # per group, matching the minimisation draw
REPLICAS = 3
NANOSECONDS = float(os.environ.get("MD_NS", 2.0))
TIMESTEP_FS = 4.0
SAMPLE_EVERY_PS = 20.0
TEMPERATURE = 310.0
CONTACT_CUTOFF = 4.5
SEED = 0

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C",
}


def chain_residues(pdb_path: Path) -> dict[str, list[tuple[int, str]]]:
    """Per-chain ordered (residue number, one-letter code) read from the file the
    simulation will actually load, so the numbering cannot disagree with it."""
    out: dict[str, list[tuple[int, str]]] = {}
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain, number, name = line[21], int(line[22:26]), line[17:20].strip()
        residues = out.setdefault(chain, [])
        if not residues or residues[-1][0] != number:
            residues.append((number, THREE_TO_ONE.get(name, "X")))
    return out


def locate_loop(pdb_path: Path, chain_key: str, loop: str) -> list[int]:
    """Residue numbers of the CDR3, found in this structure's own sequence.

    Locating the loop in the residue sequence rather than computing
    `full_seq.index(loop) + min(resid)` drops the assumption that the relaxed
    file preserved both the domain offset and the numbering origin, and fails
    loudly instead of silently returning the wrong residues."""
    residues = chain_residues(pdb_path)[chain_key]
    sequence = "".join(code for _, code in residues)
    first = sequence.find(loop)
    if first < 0:
        raise ValueError(f"CDR3 {loop} not found in chain {chain_key}")
    if sequence.find(loop, first + 1) >= 0:
        raise ValueError(f"CDR3 {loop} occurs more than once in chain {chain_key}")
    return [residues[first + offset][0] for offset in range(len(loop))]


def build_pairs(test) -> dict:
    """Select pairs and record each member's own mutated-site residue numbers.

    Kept in this script rather than a throwaway heredoc so the selection is
    reproducible and the residue mapping is reviewable next to the analysis. The
    numbers are written for provenance and for the audit in the analysis step;
    the simulation does not read them.
    """
    import importlib.util

    from immrep.data import build_training_frame

    spec = importlib.util.spec_from_file_location(
        "anat", ROOT / "scripts/26_discordant_pair_anatomy.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["anat"] = module
    spec.loader.exec_module(module)

    labels = build_training_frame()
    out: dict[str, list[dict]] = {}
    for group, want in (("discordant", (1, 0)), ("control", (1, 1))):
        # The same draw as the minimisation run, so the relaxed structures exist.
        candidates = module.find_pairs(labels, want).sample(
            n=N_PAIRS, random_state=SEED)
        entries = []
        for row in candidates.itertuples(index=False):
            if not all((RELAXED / f"{i}.pdb").exists() for i in (row.a_id, row.b_id)):
                continue
            key = "A" if row.varied_chain == "alpha" else "B"
            column = f"cdr3_{row.varied_chain}_aa"
            try:
                sites = {}
                for role in ("a", "b"):
                    tcr_id = getattr(row, f"{role}_id")
                    loop = test.loc[tcr_id, column]
                    numbers = locate_loop(RELAXED / f"{tcr_id}.pdb", key, loop)
                    sites[role] = {key: [numbers[o] for o in row.sites]}
            except ValueError as exc:
                print(f"  skipped {row.a_id}|{row.b_id}: {exc}")
                continue
            entries.append({
                "pair": f"{row.a_id}|{row.b_id}", "epitope": row.epitope,
                "a_id": row.a_id, "b_id": row.b_id,
                "a_pos": int(row.a_pos), "b_pos": int(row.b_pos),
                "varied_chain": row.varied_chain, "chain_key": key,
                "cdr3_offsets": list(row.sites),
                "a_site_residues": sites["a"], "b_site_residues": sites["b"],
            })
        out[group] = entries
    PAIRS.parent.mkdir(parents=True, exist_ok=True)
    PAIRS.write_text(json.dumps(out, indent=1))
    return out


def heavy_indices(topology, chain_id, residue_numbers=None):
    out = []
    for atom in topology.atoms():
        if atom.residue.chain.id != chain_id:
            continue
        if atom.element is not None and atom.element.symbol == "H":
            continue
        if residue_numbers is not None and int(atom.residue.id) not in residue_numbers:
            continue
        out.append(atom.index)
    return np.asarray(out, dtype=int)


def peptide_chain_id(topology) -> str | None:
    """Identify the peptide by residue count: OpenMM renamed the chains on write,
    so keying on the letter P silently finds nothing."""
    for chain in topology.chains():
        if 8 <= len(list(chain.residues())) <= 11:
            return chain.id
    return None


def run_one(path: Path, loops: dict[str, list[int]], seed: int) -> dict:
    """One trajectory. Returns per-frame contact counts and minimum distances for
    every CDR3 residue, plus the whole TCR-peptide interface contact count."""
    from openmm import LangevinMiddleIntegrator, Platform, app, unit
    from scipy.spatial import cKDTree

    pdb = app.PDBFile(str(path))
    field = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = field.createSystem(
        pdb.topology, nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=2.0 * unit.nanometer, constraints=app.HBonds,
        hydrogenMass=1.5 * unit.amu)
    integrator = LangevinMiddleIntegrator(
        TEMPERATURE * unit.kelvin, 1 / unit.picosecond,
        TIMESTEP_FS * unit.femtoseconds)
    integrator.setRandomNumberSeed(seed)
    simulation = app.Simulation(pdb.topology, system, integrator,
                                Platform.getPlatformByName("CUDA"),
                                {"Precision": "mixed"})
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy(maxIterations=200)
    simulation.context.setVelocitiesToTemperature(TEMPERATURE * unit.kelvin, seed)

    peptide = peptide_chain_id(pdb.topology)
    if peptide is None:
        return {}
    peptide_atoms = heavy_indices(pdb.topology, peptide)
    residue_atoms = {
        (key, number): heavy_indices(pdb.topology, key, {number})
        for key, numbers in loops.items() for number in numbers
    }
    residue_atoms = {k: v for k, v in residue_atoms.items() if len(v)}
    tcr_atoms = np.concatenate(
        [heavy_indices(pdb.topology, c) for c in ("A", "B")])

    steps_total = int(round(NANOSECONDS * 1e6 / TIMESTEP_FS))
    steps_between = int(round(SAMPLE_EVERY_PS * 1000 / TIMESTEP_FS))
    counts = {k: [] for k in residue_atoms}
    distances = {k: [] for k in residue_atoms}
    interface = []
    for _ in range(steps_total // steps_between):
        simulation.step(steps_between)
        state = simulation.context.getState(getPositions=True)
        xyz = np.asarray(
            state.getPositions(asNumpy=True).value_in_unit(unit.angstrom))
        tree = cKDTree(xyz[peptide_atoms])
        for key, atoms in residue_atoms.items():
            hits = tree.query_ball_point(xyz[atoms], CONTACT_CUTOFF)
            counts[key].append(sum(len(h) for h in hits))
            distances[key].append(float(tree.query(xyz[atoms])[0].min()))
        interface.append(int(sum(
            len(h) for h in tree.query_ball_point(xyz[tcr_atoms], CONTACT_CUTOFF))))
    return {"counts": counts, "distances": distances, "interface": interface}


def main() -> int:
    import pandas as pd

    from immrep.data import load_test

    test = load_test().set_index("tcr_id")
    print("resolving pairs and CDR3 numbering ...")
    pairs = build_pairs(test)
    for group, entries in pairs.items():
        print(f"  {group}: {len(entries)} pairs")

    complexes = sorted({entry[f"{role}_id"]
                        for entries in pairs.values() for entry in entries
                        for role in ("a", "b")})

    # Every CDR3 residue of both chains, so the mutated site is a later lookup.
    loops, skipped = {}, []
    for tcr_id in complexes:
        path = RELAXED / f"{tcr_id}.pdb"
        found = {}
        for key, column in (("A", "cdr3_alpha_aa"), ("B", "cdr3_beta_aa")):
            try:
                found[key] = locate_loop(path, key, test.loc[tcr_id, column])
            except (ValueError, KeyError) as exc:
                skipped.append((tcr_id, key, str(exc)))
        if found:
            loops[tcr_id] = found
    for tcr_id, key, message in skipped:
        print(f"  {tcr_id} chain {key}: {message}")

    total = len(loops) * REPLICAS
    print(f"\n{len(loops)} complexes x {REPLICAS} replicas x {NANOSECONDS} ns")
    print(f"estimated {total * NANOSECONDS * 1.385 / 60:.1f} h of integration\n")

    done = set()
    if OUT.exists():
        existing = pd.read_csv(OUT)
        done = set(zip(existing["tcr_id"], existing["replica"]))
        print(f"resuming: {len(done)} trajectories already recorded\n")

    start, finished = time.time(), 0
    for index, (tcr_id, found) in enumerate(sorted(loops.items()), 1):
        for replica in range(REPLICAS):
            if (tcr_id, replica) in done:
                continue
            try:
                result = run_one(RELAXED / f"{tcr_id}.pdb", found, seed=1000 + replica)
            except Exception as exc:
                print(f"  {tcr_id} r{replica}: {type(exc).__name__}: {str(exc)[:70]}")
                continue
            if not result:
                print(f"  {tcr_id} r{replica}: no peptide chain")
                continue

            offsets = {(key, number): offset
                       for key, numbers in found.items()
                       for offset, number in enumerate(numbers)}
            rows = []
            for key, series in result["counts"].items():
                count = np.asarray(series, dtype=float)
                distance = np.asarray(result["distances"][key], dtype=float)
                rows.append({
                    "tcr_id": tcr_id, "replica": replica, "chain": key[0],
                    "resid": key[1], "cdr3_offset": offsets[key],
                    "n_frames": len(count),
                    "count_mean": count.mean(), "count_sd": count.std(),
                    "dist_mean": distance.mean(), "dist_sd": distance.std(),
                    "dist_min": distance.min(),
                    "frac_in_contact": float((distance <= CONTACT_CUTOFF).mean()),
                })
            frame = pd.DataFrame(rows)
            OUT.parent.mkdir(parents=True, exist_ok=True)
            frame.to_csv(OUT, mode="a", header=not OUT.exists(), index=False)

            interface = np.asarray(result["interface"], dtype=float)
            pd.DataFrame({
                "tcr_id": tcr_id, "replica": replica,
                "frame": np.arange(len(interface)),
                "ps": np.arange(len(interface)) * SAMPLE_EVERY_PS,
                "interface_contacts": interface,
            }).to_csv(TRACE, mode="a", header=not TRACE.exists(), index=False)
            finished += 1

        elapsed = (time.time() - start) / 60
        rate = elapsed / finished if finished else 0
        print(f"  {index}/{len(loops)} complexes  {elapsed:.0f} min"
              f"  eta {rate * (total - len(done) - finished) / 60:.1f} h", flush=True)

    print(f"\nwrote {OUT.relative_to(ROOT)} and {TRACE.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
