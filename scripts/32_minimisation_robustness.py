#!/usr/bin/env python
"""Does the E31 interface result survive energy minimisation? (plan cards C10, C04)

E31 concluded that discordant negatives mutate at positions carrying more
peptide contact than same-label controls do, which is the project's one positive
mechanistic result. It was measured on raw Boltz output. Predicted coordinates
carry clashes and unrelaxed side chains, so a contact count read off them could
be an artefact of the pose rather than a property of the interface.

The plan's own warning is the standard to meet: if the difference is unstable to
pose or parameters, no conclusion is warranted. So this relaxes each structure
under one fixed protocol -- PDBFixer for missing terminal atoms and hydrogens,
Amber14 with GBn2 implicit solvent, L-BFGS minimisation -- and re-measures at
two minimisation depths to see whether the E31 difference moves.

What this deliberately does **not** compute is an interface binding energy.
Doing it properly needs separated-state topologies, whose termini and atom
counts differ from the complex, so the subtraction would not be comparable; and
a single-point interaction energy on an unrelaxed separated state is not a
binding free energy whatever it is called. Buried surface area, hydrogen bonds
and salt bridges are reported instead, which are measurable and interpretable.

Run: python scripts/32_minimisation_robustness.py
"""

from __future__ import annotations

import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# Leave cores for the GPU job's dataloader workers.
os.environ.setdefault("OPENMM_CPU_THREADS", "6")
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from immrep.data import build_training_frame, load_test  # noqa: E402

PREDICTIONS = ROOT / "data/interim/boltz/out/boltz_results_inputs/predictions"
OUT = ROOT / "data/interim/minimisation_robustness.csv"
REPORT = ROOT / "reports/10_minimisation_robustness.md"
RELAXED = ROOT / "data/interim/relaxed"

CONTACT_CUTOFF = 4.5
HBOND_CUTOFF = 3.5      # donor-acceptor heavy-atom distance
SALT_CUTOFF = 4.0
N_PAIRS = 20            # discordant pairs, plus the same number of controls
STAGES = (200, 1000)    # cumulative minimisation iterations
SEED = 0

CHARGED_POS = {("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"),
               ("LYS", "NZ"), ("HIS", "ND1"), ("HIS", "NE2")}
CHARGED_NEG = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
POLAR = set("NO")


def interface_descriptors(positions, topology) -> dict[str, float]:
    """Contacts, hydrogen-bond-like and salt-bridge-like pairs, and buried area
    between the TCR chains (A, B) and the peptide (P)."""
    import freesasa
    from scipy.spatial import cKDTree

    coords, meta = [], []
    for atom in topology.atoms():
        if atom.element is not None and atom.element.symbol == "H":
            continue
        chain = atom.residue.chain.id
        coords.append(positions[atom.index])
        meta.append((chain, atom.residue.name, atom.name, atom.element.symbol
                     if atom.element else "X"))
    xyz = np.asarray([[c.x * 10, c.y * 10, c.z * 10] for c in coords])  # nm -> A
    chains = np.array([m[0] for m in meta])

    tcr = np.isin(chains, ["A", "B"])
    peptide = chains == "P"
    if not tcr.any() or not peptide.any():
        return {}

    tree = cKDTree(xyz[tcr])
    pairs = tree.query_ball_point(xyz[peptide], CONTACT_CUTOFF)
    contacts = int(sum(len(p) for p in pairs))

    tcr_meta = [m for m, keep in zip(meta, tcr) if keep]
    tcr_xyz = xyz[tcr]
    peptide_indices = np.where(peptide)[0]
    hbonds = salts = 0
    for slot, (atom_index, hits) in enumerate(zip(peptide_indices, pairs)):
        _, pres, pname, pelem = meta[atom_index]
        here = xyz[atom_index]
        for hit in hits:
            _, res, name, elem = tcr_meta[hit]
            distance = float(np.linalg.norm(here - tcr_xyz[hit]))
            if elem in POLAR and pelem in POLAR and distance <= HBOND_CUTOFF:
                hbonds += 1
            if distance <= SALT_CUTOFF and (
                ((res, name) in CHARGED_POS and (pres, pname) in CHARGED_NEG)
                or ((res, name) in CHARGED_NEG and (pres, pname) in CHARGED_POS)
            ):
                salts += 1

    # Buried area: SASA of the parts minus SASA of the whole, over the peptide.
    def sasa(mask) -> float:
        sub = xyz[mask]
        radii = [1.8 if m[3] == "C" else 1.5 for m, keep in zip(meta, mask) if keep]
        result = freesasa.calcCoord(sub.flatten().tolist(), radii)
        return float(result.totalArea())

    buried = sasa(peptide) + sasa(tcr) - sasa(peptide | tcr)
    return {"contacts": contacts, "hbonds": hbonds, "salt_bridges": salts,
            "buried_area": buried}


def minimise(path: Path) -> list[dict]:
    from openmm import LangevinMiddleIntegrator, Platform, app, unit
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(path))
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # never invent unresolved residues
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)      # one fixed protonation assumption

    field = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = field.createSystem(
        fixer.topology, nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds)
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    simulation = app.Simulation(fixer.topology, system, integrator,
                                Platform.getPlatformByName("CPU"))
    simulation.context.setPositions(fixer.positions)

    out = []
    done = 0
    for stage in (0,) + STAGES:
        if stage:
            simulation.minimizeEnergy(maxIterations=stage - done)
            done = stage
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        record = {"stage": stage,
                  "energy": state.getPotentialEnergy().value_in_unit(
                      unit.kilojoule_per_mole)}
        record.update(interface_descriptors(state.getPositions(), fixer.topology))
        out.append(record)
        # Keep the relaxed coordinates. The first version of this script stored
        # only whole-interface descriptors, which cannot test E31 at all: E31
        # measures contact at the *mutated CDR3 position*, and both members of a
        # pair are near-identical TCRs whose whole interfaces are bound to look
        # alike. Saving the structures makes the position-specific test possible
        # without relaxing everything again.
        if stage == STAGES[-1]:
            RELAXED.mkdir(parents=True, exist_ok=True)
            with open(RELAXED / f"{path.parent.name}.pdb", "w") as handle:
                app.PDBFile.writeFile(fixer.topology, state.getPositions(), handle)
    return out


def main() -> int:
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "anat", ROOT / "scripts/26_discordant_pair_anatomy.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["anat"] = module
    spec.loader.exec_module(module)

    labels = build_training_frame()
    discordant = module.find_pairs(labels, (1, 0)).sample(
        n=N_PAIRS, random_state=SEED)
    control = module.find_pairs(labels, (1, 1)).sample(n=N_PAIRS, random_state=SEED)

    jobs = []
    for group, frame in (("discordant", discordant), ("control", control)):
        for row in frame.itertuples(index=False):
            for role, tcr_id in (("a", row.a_id), ("b", row.b_id)):
                jobs.append({"group": group, "pair": f"{row.a_id}|{row.b_id}",
                             "role": role, "tcr_id": tcr_id,
                             "epitope": row.epitope,
                             "pos": row.a_pos if role == "a" else row.b_pos})
    jobs = pd.DataFrame(jobs).drop_duplicates("tcr_id").to_dict("records")
    print(f"{len(jobs)} complexes to relax "
          f"({sum(j['group'] == 'discordant' for j in jobs)} discordant members)")
    print(f"estimated {len(jobs) * 50 / 60:.0f} minutes on CPU\n")

    rows, start = [], time.time()
    for index, job in enumerate(jobs, 1):
        path = PREDICTIONS / job["tcr_id"] / f"{job['tcr_id']}_model_0.pdb"
        if not path.exists():
            continue
        try:
            for record in minimise(path):
                rows.append({**job, **record})
        except Exception as exc:                      # structure-specific failures
            print(f"  {job['tcr_id']}: {type(exc).__name__}: {str(exc)[:80]}")
        if index % 10 == 0:
            print(f"  {index}/{len(jobs)}  ({(time.time() - start) / 60:.0f} min elapsed)")

    frame = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT, index=False)
    print(f"\nwrote {OUT.relative_to(ROOT)}  {frame.tcr_id.nunique()} complexes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
