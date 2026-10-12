# Certificado ATK v2 — ShantiDraconis/m-core-openai-ns-audit

- Commit verificado: `ce4b23683c03572002092b69a3eb5660fa32e256` (tag `atk-kernel-cert-v2`)
- Gerado em: 2026-10-12T00:46:28+00:00
- Teoremas aceitos pelo kernel: **6/6**

## ATK kernel check
- Execução: https://github.com/ShantiDraconis/m-core-openai-ns-audit/actions/runs/38188856695
- Lean (version 4.34.1, x86_64-unknown-linux-gnu, commit 5045d0056413266e57c625dcd7c365b10e377c52, Release)
- SHA-256 das fontes: `92e3b2ef9b10eb7ba602fe7f0c2bb1d89190cf59d8a45cb7262d8b58dc74ecec`
- Job “Navier–Stokes — Taylor–Green, razão S/(νz)”: success
- Job “Lean 4 + Mathlib — compilação e auditoria de axiomas”: success

| Teorema | Veredito | Axiomas |
|---|---|---|
| `ATK.NavierStokes.absorption` | ✅ CERTIFIED  | propext, Classical.choice, Quot.sound |
| `ATK.NavierStokes.absorption_margin_pos` | ✅ CERTIFIED  | propext, Classical.choice, Quot.sound |
| `ATK.NavierStokes.riccati_comparison` | ✅ CERTIFIED  | propext, Classical.choice, Quot.sound |
| `ATK.NavierStokes.riccati_upper_bound` | ✅ CERTIFIED  | propext, Classical.choice, Quot.sound |
| `ATK.NavierStokes.absorption_without_remainder_forces_decay` | ✅ CERTIFIED  | propext, Classical.choice, Quot.sound |
| `ATK.NSControl.weighted_energy_control` | ✅ CERTIFIED  | propext, Classical.choice, Quot.sound |

## Resultados numéricos registrados
- **NS Taylor–Green Re=400**: max S/(nu z) = 11.457533 em t = 1.35; S > nu z em 75.1% do tempo; C_min(0.7) = 6.876468e-04; C_min(0.9) = 6.719851e-04
- **NS Taylor–Green Re=200**: max S/(nu z) = 5.672057 em t = 1.35; S > nu z em 58.2% do tempo; C_min(0.7) = 6.414571e-04; C_min(0.9) = 6.108504e-04
- **NS Taylor–Green Re=100**: max S/(nu z) = 2.776325 em t = 1.40; S > nu z em 45.3% do tempo; C_min(0.7) = 5.478257e-04; C_min(0.9) = 4.881038e-04
- **NS Taylor–Green Re=100**: max S/(nu z) = 2.776325 em t = 1.40; S > nu z em 88.5% do tempo; C_min(0.7) = 5.478257e-04; C_min(0.9) = 4.881038e-04

Lemas auxiliares e cálculos concretos. Nenhuma das sete conjecturas foi provada; OPEN_MATH=true.
