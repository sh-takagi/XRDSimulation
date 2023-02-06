import numpy as np
import pandas as pd
import xrayutilities as xu
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from tempfile import NamedTemporaryFile
import matplotlib.pyplot as plt
import optuna
import sys

# cifファイルから原子と位置情報を取得。
def getAtomsPos(path):
    data = []
    with open(path, mode="r", encoding="utf-8") as f:
        data = f.readlines()

    is_atom = False
    atoms = []
    pos = []
    for d in data:
        if d == "loop_\n" and is_atom: break
        if is_atom:
            d_split = d.split()
            atoms.append(d_split[1])
            pos.append(d_split[2:5])
        if d == "_atom_site_occupancy\n": is_atom = True

    atoms = np.array(atoms)
    pos = np.array(pos).astype(np.float64)
    return atoms, pos

# 実験データからtargetを生成。
def csv2target(path):
    df = pd.read_csv(path)
    return np.array(df[" yobs"] - df[" bkg"])

# 格子、原子、位置情報からXRDパターンをシミュレーションする。
def simXRD(lattice, atoms, pos):
    latt = Lattice.from_parameters(
            lattice["a"],
            lattice["b"],
            lattice["c"],
            lattice["alpha"],
            lattice["beta"],
            lattice["gamma"]
    )

    structure = Structure(latt, atoms, pos)

    tmp_cif = NamedTemporaryFile(delete=False)
    structure.to("cif", tmp_cif.name)
    xu_cif = xu.materials.cif.CIFFile(tmp_cif.name)
    xu_crystal = xu.materials.material.Crystal(name="model", lat=xu_cif.SGLattice())
    tmp_cif.close()

    twotheta = np.arange(5, 40.01, 0.02)

    powder = xu.simpack.Powder(xu_crystal, 1)
    pm = xu.simpack.PowderModel(powder, I0=100)
    intensity = pm.simulate(twotheta)

    pm.close()

    return intensity

# 2つのベクトルのコサイン類似度を計算。
def cosSimilary(v1, v2):
    v1_norm = np.linalg.norm(v1)
    v2_norm = np.linalg.norm(v2)

    return np.dot(v1, v2) / (v1_norm * v2_norm)

# optunaのobjective
def objective(trial, atoms, pos, target):
    min_a = 18
    max_a = 26
    min_b = 18
    max_b = 26
    min_c = 3
    max_c = 6
    min_angle = 60
    max_angle = 120

    params = {
            "a": trial.suggest_float("a", min_a, max_a),
            "b": trial.suggest_float("b", min_b, max_b),
            "c": trial.suggest_float("c", min_c, max_c),
            "alpha": trial.suggest_float("alpha", min_angle, max_angle),
            "beta": trial.suggest_float("beta", min_angle, max_angle),
            "gamma": trial.suggest_float("gamma", min_angle, max_angle)
    }

    sim = simXRD(params, atoms, pos)
    score = cosSimilary(sim, target)
    return score

def main():
    n_trials = int(sys.argv[1])

    atoms, pos = getAtomsPos("structure.cif")
    target = csv2target("target.csv")

    study = optuna.create_study(study_name="structure", storage="sqlite:///structure.db", load_if_exists=True, direction="maximize")
    study.optimize(lambda trial: objective(trial, atoms, pos, target), n_trials=n_trials, gc_after_trial=True)

    print("Best Score", study.best_params, sep="\n")

if __name__ == "__main__":
    main()
