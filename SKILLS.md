# AI Skills for Lanthanide Multireference Workflows

Purpose: equip the agent to emit concise, runnable Python that leverages **PySCF** and the lab's **embed_sim** / DMET tooling for lanthanide (4f–5d) electronic-structure tasks.

## Style and expectations
- Use library functions; do not re-implement SCF, integral generation, or solvers.
- Prefer small helper functions (build molecule, run SCF, active-space solver, post-CAS correction).
- Keep comments minimal and focused on intent. Favor descriptive variable names instead.
- Return key objects (mol, mf, mc, embedding object) so downstream steps can reuse orbitals and checkpoints.

## Inputs the agent should capture
- Geometry in Å (list of `(element, (x, y, z))` tuples), charge, and spin (`mol.spin = 2S`).
- Choice of basis/ECP suitable for lanthanides: e.g., `def2-TZVPP` + `def2-ECP`, `SARC-DKH-TZ`, or `ANO-RCC...` (use relativistic-consistent options; add diffuse functions for excited states).
- Symmetry: enable when using CASSCF/DMRG (`symmetry=True`) to keep irreps well-labeled.
- Active space hints: include the 4f shell (7 orbitals) and add 5d/6s as needed for 4f→5d excitations; set electron count per oxidation state.
- Roots/state averaging: number of states and weights, especially for spin–orbit / crystal-field analyses.

## Reusable patterns

### Build molecule and SCF
```python
from pyscf import gto, scf

def build_scf(atom, charge=0, spin=0, basis="def2-TZVPP", ecp="def2-tzvpp", symmetry=True):
    mol = gto.Mole()
    mol.build(atom=atom, charge=charge, spin=spin, basis=basis, ecp=ecp,
              unit="Angstrom", symmetry=symmetry, verbose=4)
    mf = scf.ROHF(mol).density_fit()
    mf.conv_tol = 1e-8
    mf.kernel()
    return mol, mf
```

### CASSCF (or DMRG-CASSCF) with state averaging
```python
from pyscf import mcscf
from pyscf.dmrgscf import DMRGSCF

def run_casscf(mf, ncas, nelecas, roots=1, weights=None, use_dmrg=False, maxM=2000):
    solver = DMRGSCF(mf, ncas, nelecas) if use_dmrg else mcscf.CASSCF(mf, ncas, nelecas)
    if use_dmrg:
        solver.fcisolver.maxM = maxM
    if roots > 1:
        solver = solver.state_average_(weights or [1 / roots] * roots)
    solver.conv_tol = 1e-8
    solver.kernel()
    return solver
```

### Dynamic correlation (NEVPT2 / PEVPT2)
```python
from pyscf import mrpt

def add_nevpt2(mc, root=0):
    nevpt = mrpt.NEVPT(mc, root=root)
    return nevpt.kernel()
```

### DMET / embedding notes
- Use SCF orbitals as the low-level bath builder; keep density-fitting to accelerate Coulomb terms.
- Define fragments centered on the lanthanide plus coordinating ligands; include valence 4f/5d functions in impurity projectors.
- When using `embed_sim`, prefer its high-level drivers (e.g., a DMET/DMRG class) and pass a solver factory that returns a CASSCF or DMRG-CASSCF object for each fragment.
- For PySCF's built-in DMET, a common pattern is:
```python
from pyscf.dmet import sdmfet

def run_dmet(mf, fragments, ncas, nelecas):
    emb = sdmfet.SDMFET(mf, fragments=fragments, impbasis=None,
                        solver='CASSCF', cas_norb=ncas, cas_nelec=nelecas,
                        max_cycle=30, conv_tol=1e-7)
    emb.kernel()
    return emb
```
- Save fragment orbitals and bath information so subsequent excited-state or property calculations can reuse them.

## Minimal end-to-end template (4f→5d excitation)
```python
atom = [
    ("Ce", (0.0, 0.0, 0.0)),
    ("Cl", (0.0, 0.0, 2.6)),
    ("Cl", (0.0, 0.0, -2.6)),
]

mol, mf = build_scf(atom, charge=3, spin=1, basis="def2-TZVPP", ecp="def2-tzvpp")
mc = run_casscf(mf, ncas=12, nelecas=8, roots=3, weights=[0.34, 0.33, 0.33])
e_nevpt2 = add_nevpt2(mc)
print("SA-CASSCF energies:", mc.e_states)
print("NEVPT2 correction:", e_nevpt2)
```

## When information is missing
- Ask for: oxidation state (charge/spin), target states and multiplicities, chosen basis/ECP, and whether to use DMET vs direct CASSCF/DMRG.
- If unspecified, default to ROHF + CASSCF + NEVPT2 with `def2-TZVPP/def2-ECP` and a 4f-only active space, noting assumptions in a brief leading comment.
