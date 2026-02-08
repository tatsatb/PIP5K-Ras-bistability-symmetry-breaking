# Main Files

All the main files are named `run_for_(phenotype)_kymo.m`.  
These filenames are self-explanatory: For example, to obtain the kymograph for the WT phenotype, run `run_for_WT_kymo.m`.

---

## Auxiliary files

- `magma.m` – provides the magma colormap.
- `nabla_SDE.m` – discretization of the spatial derivative.
- `nullclines.m` – computes the nullclines and equilibrium values of the system.
- `kymo_sdefile.m`, `kymo2_sdefile.m`, and `SDE_eduler_deb.m` – auxiliary files for the SDE toolbox.

---

## Difference between kymo_sdefile and kymo2_sdefile

The only difference is in the definition of the J₂_Ras term. In the first version:

```matlab
J2_Ras = (p.a3./((p.a4.^2*PIP2.^2 + 1)) + p.a5).*(1 + p.a*(Actin-p.Tmem*Tmem));
```

In the second version, the expression is constrained to be non-negative (which is also mostly true in the first case):

```matlab
J2_Ras = (p.a3./((p.a4.^2*PIP2.^2 + 1)) + p.a5).*max((1 + p.a*Actin - p.Tmem*Tmem), 0);
```

Thus, in `kymo2_sdefile`, the term is explicitly forced to remain positive.
