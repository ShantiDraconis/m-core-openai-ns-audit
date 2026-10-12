"""Parâmetros NS 49–53 sobre as soluções de Taylor–Green (64³): integral do resto e estados degenerados.

y = ‖ω‖² (enstrofia), z = ‖∇ω‖² (palinstrofia), S = ∫ω·(∇u)ω, R_min(θ) = max(0, S − θνz)
(o menor resto compatível com um orçamento bulk + tail ≤ θνz). Integrais pela regra do trapézio (Δt = 0.05).
uso: python3 remainder.py
"""
import csv
for Re in (100, 200, 400):
    nu = 1 / Re
    R = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(f"tg_N64_Re{Re}.csv"))]
    t = [r["t"] for r in R]; y = [r["y_enstrofia"] for r in R]
    z = [r["z_palinstrofia"] for r in R]; S = [r["S_stretch"] for r in R]
    trap = lambda f: sum((t[i + 1] - t[i]) * (f[i] + f[i + 1]) / 2 for i in range(len(t) - 1))
    Rmin = lambda th: [max(0.0, S[i] - th * nu * z[i]) for i in range(len(t))]
    print(f"Re={Re}: T={t[-1]}  min y={min(y):.4g}  min νz={nu * min(z):.4g}  "
          f"∫y² dt={trap([v * v for v in y]):.6g}  ∫R_min(0.7) dt={trap(Rmin(0.7)):.6g}  "
          f"∫R_min(0.9) dt={trap(Rmin(0.9)):.6g}  max S/(νz)={max(S[i] / (nu * z[i]) for i in range(len(t))):.4f}")
