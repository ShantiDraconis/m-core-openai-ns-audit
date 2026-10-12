"""ATK — correção NS-SCALING-FALSE (12/10/2026): itens 19, 20, 21 falsos; 52 cai por dependência.

Campo semente de divergência nula, localizado: v = rot ψ, ψ = gaussiana × polinômio (σ = 0.5).
Mede y = ‖ω‖², z = ‖∇ω‖², S = ∫ω·(∇v)ω por FFT e VERIFICA NUMERICAMENTE as duas leis de escala:
  amplitude  v ↦ a·v          : y·a², z·a², S·a³
  dilatação  v ↦ λ·v(λx)      : y·λ,  z·λ³, S·λ³    (implementada como ψ ↦ ψ(λx), pois rot[ψ(λ·)] = λ(rot ψ)(λ·))
e as consequências:
  S/(νz) = a·S₀/(νz₀) → ∞                          ⇒ c_B* + c_T* = +∞   (itens 20, 21)
  (S − θνz)/y² → 0 em amplitude (máximo em a* = 2θνz₀/S₀), mas ∝ λ em dilatação
                                                    ⇒ nenhum C com S ≤ θνz + C·y²   (item 19)
uso: python3 scaling_counterexample.py [n]   (n = 64, 96 ou 128; padrão 128)
"""
import sys
import numpy as np

n = int(sys.argv[1]) if len(sys.argv) > 1 else 128
L = 2 * np.pi
x = np.arange(n) * L / n - np.pi
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
k = np.fft.fftfreq(n, d=L / n) * 2 * np.pi
KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
D = lambda f, K: np.real(np.fft.ifftn(1j * K * np.fft.fftn(f)))
dV = (L / n) ** 3


def functionals(lam=1.0, amp=1.0, sig=0.5):
    Xl, Yl, Zl = lam * X, lam * Y, lam * Z
    g = np.exp(-(Xl**2 + Yl**2 + Zl**2) / sig**2)
    psi = [g * (1 + Xl + 0.5 * Yl), g * (Yl - 0.3 * Zl + Xl * Zl), g * (0.7 * Zl + Xl * Yl)]
    v = [amp * (D(psi[2], KY) - D(psi[1], KZ)), amp * (D(psi[0], KZ) - D(psi[2], KX)),
         amp * (D(psi[1], KX) - D(psi[0], KY))]
    w = [D(v[2], KY) - D(v[1], KZ), D(v[0], KZ) - D(v[2], KX), D(v[1], KX) - D(v[0], KY)]
    grad = [[D(v[i], K) for K in (KX, KY, KZ)] for i in range(3)]
    div = float(np.abs(sum(grad[i][i] for i in range(3))).max())
    y = sum((w[i] ** 2).sum() for i in range(3)) * dV
    z = sum((D(w[i], K) ** 2).sum() for i in range(3) for K in (KX, KY, KZ)) * dV
    S = sum((w[i] * grad[i][j] * w[j]).sum() for i in range(3) for j in range(3)) * dV
    edge = max(float(np.abs(v[i][0]).max()) for i in range(3))
    return y, z, S, div, edge


y0, z0, S0, div0, edge0 = functionals()
print(f"n={n}  semente: y={y0:.6g} z={z0:.6g} S={S0:.6g}  |div v|max={div0:.1e}  |v| na borda={edge0:.1e}")

ok = True
print("\nLei de amplitude (medida):")
for a in (10.0, 100.0):
    y, z, S, *_ = functionals(amp=a)
    ey, ez, eS = (np.log(y / y0) / np.log(a), np.log(z / z0) / np.log(a), np.log(S / S0) / np.log(a))
    print(f"  a={a:>5}: expoentes medidos y:{ey:.6f} z:{ez:.6f} S:{eS:.6f}   (esperado 2, 2, 3)")
    ok &= abs(ey - 2) < 1e-6 and abs(ez - 2) < 1e-6 and abs(eS - 3) < 1e-6

print("Lei de dilatação (medida, ψ ↦ ψ(λx)):")
for lam in (1.5, 2.0):
    y, z, S, div, edge = functionals(lam=lam)
    ey, ez, eS = (np.log(y / y0) / np.log(lam), np.log(z / z0) / np.log(lam), np.log(S / S0) / np.log(lam))
    print(f"  λ={lam}: expoentes medidos y:{ey:.6f} z:{ez:.6f} S:{eS:.6f}   (esperado 1, 3, 3)  |v| na borda={edge:.1e}")
    ok &= abs(ey - 1) < 1e-4 and abs(ez - 3) < 1e-4 and abs(eS - 3) < 1e-4

nu, th = 1.0, 0.9
print("\nItens 20/21 — S/(νz) sob amplitude (ν = 1):")
for a in (1.0, 1e4, 1e6):
    print(f"  a={a:.0e}: S/(νz) = {a * S0 / (nu * z0):.6g}")
astar = 2 * th * nu * z0 / S0
print(f"\nItem 19 — (S − θνz)/y² (θ = {th}, ν = 1):")
for a in (2e4, astar, 1e6):
    r = (a**3 * S0 - th * nu * a**2 * z0) / (a**2 * y0) ** 2
    print(f"  amplitude a={a:.6g}: razão = {r:.4e}")
rmax = (astar**3 * S0 - th * nu * astar**2 * z0) / (astar**2 * y0) ** 2
print(f"  máximo em amplitude: a* = 2θνz₀/S₀ = {astar:.6g}, razão = {rmax:.4e} (a amplitude sozinha NÃO diverge)")
for lam in (1.0, 1e3, 1e9, 1e12):
    print(f"  a = a*, dilatação λ={lam:.0e}: razão = {lam * rmax:.4e}   (∝ λ)")
print(f"  para superar uma constante C basta λ > C / {rmax:.4e}")
print("\nLEIS DE ESCALA CONFIRMADAS" if ok else "\nFALHA NA VERIFICAÇÃO DAS LEIS DE ESCALA")
sys.exit(0 if ok else 1)
