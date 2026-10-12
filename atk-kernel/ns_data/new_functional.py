"""ATK — FICHA, passo 6 (teste numérico obrigatório antes do Lean) para funcionais candidatos F.

Campo semente: v = rot ψ, ψ = gaussiana × polinômio (σ = 0.5), divergência nula, localizado.
Para cada F mede, NO CAMPO (não por fórmula):
  * expoente de amplitude  (v ↦ a·v, a = 10)
  * expoente de dilatação  (v ↦ λ·v(λx), λ = 1.5 e 2, via ψ ↦ ψ(λx))
  * degeneração: |F| / (y^{3/4} z^{3/4})  (F ≡ 0 não carrega informação)
  * tautologia: F / S  (F proporcional a S é "target laundering": controlar F é controlar S)
Classificação (ficha ATK):
  DEGENERADO        F ≈ 0 em todo campo (ex.: lagrangiano nulo)
  TAUTOLOGICO       F = c·S identicamente
  REJEITADO_ESCALA  vetor de escala de F fora da reta {λ = 3, a ≥ 3} (a cota S ≤ θνz + C·F quebra)
  PASSA_ESCALA      condição necessária ok; falta CONTROLE DINÂMICO de F ao longo da solução (aberto)
uso: python3 new_functional.py [n]   (padrão 128)
"""
import sys
import numpy as np

n = int(sys.argv[1]) if len(sys.argv) > 1 else 128
L = 2 * np.pi
x = np.arange(n) * L / n - np.pi
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
k = np.fft.fftfreq(n, d=L / n) * 2 * np.pi
KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
K3 = (KX, KY, KZ)
D = lambda f, K: np.real(np.fft.ifftn(1j * K * np.fft.fftn(f)))
dV = (L / n) ** 3
I = lambda f: float(f.sum() * dV)


def state(lam=1.0, amp=1.0, sig=0.5):
    Xl, Yl, Zl = lam * X, lam * Y, lam * Z
    g = np.exp(-(Xl**2 + Yl**2 + Zl**2) / sig**2)
    psi = [g * (1 + Xl + 0.5 * Yl), g * (Yl - 0.3 * Zl + Xl * Zl), g * (0.7 * Zl + Xl * Yl)]
    u = [amp * (D(psi[2], KY) - D(psi[1], KZ)), amp * (D(psi[0], KZ) - D(psi[2], KX)),
         amp * (D(psi[1], KX) - D(psi[0], KY))]
    G = [[D(u[i], K) for K in K3] for i in range(3)]          # G[i][j] = ∂_j u_i
    w = [G[2][1] - G[1][2], G[0][2] - G[2][0], G[1][0] - G[0][1]]
    lap = [sum(D(D(u[i], K), K) for K in K3) for i in range(3)]
    y = I(sum(w[i] ** 2 for i in range(3)))
    z = I(sum(D(w[i], K) ** 2 for i in range(3) for K in K3))
    S = I(sum(w[i] * G[i][j] * w[j] for i in range(3) for j in range(3)))
    conv = [sum(u[j] * G[i][j] for j in range(3)) for i in range(3)]      # (u·∇)u
    F = {}
    F["S_alt = -∫(u·∇)u·Δu"] = -I(sum(conv[i] * lap[i] for i in range(3)))
    det = (G[0][0] * (G[1][1] * G[2][2] - G[1][2] * G[2][1])
           - G[0][1] * (G[1][0] * G[2][2] - G[1][2] * G[2][0])
           + G[0][2] * (G[1][0] * G[2][1] - G[1][1] * G[2][0]))
    F["F_det = ∫det(∇u)"] = I(det)   # lagrangiano nulo: ∫det(∇u) = 0 para todo campo periódico/compacto
    F["F_vs = ∫ω·Sym(∇u)·ω"] = I(sum(w[i] * 0.5 * (G[i][j] + G[j][i]) * w[j] for i in range(3) for j in range(3)))
    F["F_L3 = ∫|∇u|³"] = I(sum(G[i][j] ** 2 for i in range(3) for j in range(3)) ** 1.5)
    F["F_crit = ‖u‖_L3·z"] = I(sum(u[i] ** 2 for i in range(3)) ** 1.5) ** (1 / 3) * z
    F["GN = y^(3/4)·z^(3/4)"] = y ** 0.75 * z ** 0.75
    F["T_log = GN/log(e+z/y)"] = y ** 0.75 * z ** 0.75 / np.log(np.e + z / y)
    F["F_log = GN/√log(e+z/y)"] = y ** 0.75 * z ** 0.75 / np.sqrt(np.log(np.e + z / y))
    F["y² (rota 19, falsa)"] = y ** 2
    divmax = float(np.abs(sum(G[i][i] for i in range(3))).max())
    return y, z, S, F, divmax


y0, z0, S0, F0, div0 = state()
print(f"n={n}  semente: y={y0:.6g} z={z0:.6g} S={S0:.6g}  |div u|max={div0:.1e}")
ya, za, Sa, Fa, _ = state(amp=10.0)
dil = {lam: state(lam=lam) for lam in (1.5, 2.0)}
GN0 = F0["GN = y^(3/4)·z^(3/4)"]

print(f"\n{'funcional':<26}{'valor':>14}{'exp. a':>9}{'exp. λ(1.5)':>13}{'exp. λ(2)':>11}{'|F|/GN':>11}{'F/S':>11}  veredito")
for name, v in F0.items():
    def ex(v1, v0, s):
        return np.log(abs(v1) / abs(v0)) / np.log(s) if abs(v0) > 1e-12 * GN0 and v1 != 0 else float("nan")
    ea = ex(Fa[name], v, 10.0)
    el15 = ex(dil[1.5][3][name], v, 1.5)
    el2 = ex(dil[2.0][3][name], v, 2.0)
    rel = abs(v) / GN0
    tau = v / S0
    if rel < 1e-10:
        ver = "DEGENERADO (F ≡ 0)"
    elif abs(tau - round(tau)) < 1e-8 and round(tau) != 0 and name.startswith(("F_vs", "S_alt")):
        ver = f"TAUTOLOGICO (F = {round(tau):+d}·S)"
    elif abs(el2 - 3) < 1e-3 and abs(el15 - 3) < 1e-3 and ea >= 3 - 1e-6:
        ver = "PASSA_ESCALA (falta controle dinâmico)"
    else:
        ver = "REJEITADO_ESCALA"
    print(f"{name:<26}{v:>14.6g}{ea:>9.4f}{el15:>13.4f}{el2:>11.4f}{rel:>11.3e}{tau:>11.4g}  {ver}")

# termo com log: expoente de dilatação local < 3 ⇒ (S − θνz)/T_log cresce como log λ
nu, th = 1.0, 0.9
astar = 2 * th * nu * z0 / S0
Dst = astar**3 * S0 - th * nu * astar**2 * z0
GNst = (astar**2 * y0) ** 0.75 * (astar**2 * z0) ** 0.75
print("\nCandidato 1, termo C·GN/log(e+z/y): (S − θνz)/T_log ao longo da dilatação (a = a*):")
for lam in (1.0, 1e3, 1e9, 1e30, 1e100):
    r = Dst / GNst * np.log(np.e + lam**2 * z0 / y0)   # λ³ cancela; sobra o log
    print(f"  λ={lam:.0e}: razão = {r:.4e}")
print("  → cresce como log λ: nenhum C serve (mesma falha da rota 19, mais lenta).")
print("\nCorpo 1 (F_log = GN/√log(e+z/y)): (S − θνz)/F_log ao longo da dilatação (a = a*):")
for lam in (1.0, 1e3, 1e9, 1e30, 1e100):
    r = Dst / GNst * np.sqrt(np.log(np.e + lam**2 * z0 / y0))
    print(f"  λ={lam:.0e}: razão = {r:.4e}")
print("  → cresce como √log λ: dividir por qualquer log que cresça com z/y quebra a escala λ³.")
print(f"\nz/y no campo semente = {z0 / y0:.6g}; sob amplitude z/y não muda (a² / a²); sob dilatação cresce como λ².")
