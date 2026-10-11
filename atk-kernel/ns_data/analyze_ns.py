"""Resumo das medidas de Taylor–Green: razão S/(ν z) e constante mínima do resto.

Para uma divisão S = bulk + tail + R com bulk + tail ≤ θ ν z (θ = c_B + c_T < 1), o resto
precisa satisfazer R ≥ S − θ ν z. Se o programa exige R ≤ C y², a menor constante
compatível com esta solução é
    C_min(θ) = max_t max(0, S − θ ν z) / y².
Uso: python3 analyze_ns.py Re arquivo.csv [mais pares Re arquivo.csv ...]
"""
import csv, os, sys

args = sys.argv[1:]
lines = ["| Re | max S/(νz) | t do máximo | % do tempo com S > νz | C_min(θ=0.7) | C_min(θ=0.9) |",
         "|---|---|---|---|---|---|"]
notices = []
for Re, path in zip(args[0::2], args[1::2]):
    nu = 1 / float(Re)
    rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(path))]
    rat = [(r["razao_S_sobre_nu_z"], r["t"]) for r in rows]
    mx, tmx = max(rat)
    frac = 100 * sum(1 for v, _ in rat if v > 1) / len(rat)
    def cmin(theta):
        return max(max(0.0, r["S_stretch"] - theta * nu * r["z_palinstrofia"]) / r["y_enstrofia"] ** 2 for r in rows)
    lines.append(f"| {Re} | {mx:.4f} | {tmx:.2f} | {frac:.1f}% | {cmin(0.7):.3e} | {cmin(0.9):.3e} |")
    notices.append(f"::notice title=NS Taylor–Green Re={Re}::max S/(nu z) = {mx:.6f} em t = {tmx:.2f}%0A"
                   f"S > nu z em {frac:.1f}% do tempo%0AC_min(0.7) = {cmin(0.7):.6e}%0AC_min(0.9) = {cmin(0.9):.6e}")
out = "\n".join(["## Navier–Stokes — Taylor–Green: S/(νz) numa solução real", ""] + lines + ["",
      "Se max S/(νz) > 1, nenhuma divisão bulk + tail = S satisfaz c_B + c_T < 1 nesse instante:",
      "as hipóteses de absorção sem resto implicam enstrofia decrescente, o que esta solução contradiz."])
print(out)
if os.environ.get("GITHUB_ACTIONS"):
    print("\n".join(notices))
    open(os.environ["GITHUB_STEP_SUMMARY"], "a").write(out + "\n")
