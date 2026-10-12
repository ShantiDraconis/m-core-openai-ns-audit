"""ATK — parâmetros NS 19–21 e 52: as cotas universais S ≤ θνz + C·y² e c_B* + c_T* < 1 são falsas.
Campo de divergência nula, suportado (numericamente) numa bola pequena: v = rot ψ, ψ gaussiana × polinômio.
Calcula y = ‖ω‖², z = ‖∇ω‖², S = ∫ω·(∇v)ω por FFT e a escala λ: v_λ(x) = λ v(λx)."""
import numpy as np
n, L = 96, 2 * np.pi
x = np.arange(n) * L / n - np.pi
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
k = np.fft.fftfreq(n, d=L / n) * 2 * np.pi
KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
D = lambda f, K: np.real(np.fft.ifftn(1j * K * np.fft.fftn(f)))
def fields(sig):
    g = np.exp(-(X**2 + Y**2 + Z**2) / sig**2)
    psi = [g * (1 + X + 0.5 * Y), g * (Y - 0.3 * Z + X * Z), g * (0.7 * Z + X * Y)]
    v = [D(psi[2], KY) - D(psi[1], KZ), D(psi[0], KZ) - D(psi[2], KX), D(psi[1], KX) - D(psi[0], KY)]
    w = [D(v[2], KY) - D(v[1], KZ), D(v[0], KZ) - D(v[2], KX), D(v[1], KX) - D(v[0], KY)]
    dV = (L / n) ** 3
    grad = [[D(v[i], K) for K in (KX, KY, KZ)] for i in range(3)]
    div = max(abs(sum(grad[i][i] for i in range(3))).ravel())
    y = sum((w[i] ** 2).sum() for i in range(3)) * dV
    z = sum((D(w[i], K) ** 2).sum() for i in range(3) for K in (KX, KY, KZ)) * dV
    S = sum((w[i] * grad[i][j] * w[j]).sum() for i in range(3) for j in range(3)) * dV
    edge = max(abs(v[i][0]).max() for i in range(3))
    return y, z, S, div, edge
for sig in (0.5, 0.35):
    y, z, S, div, edge = fields(sig)
    print(f"σ={sig}: y={y:.6g} z={z:.6g} S={S:.6g} |div v|max={div:.2e} |v| na borda={edge:.1e}")
    # escala: S, z ~ λ³, y ~ λ (ℝ³); amplitude a: S ~ a³, z, y ~ a²
    nu, th = 1.0, 0.9
    a = 2 * th * nu * z / abs(S)
    sgn = np.sign(S)
    Sa, za, ya = sgn * a**3 * S, a**2 * z, a**2 * y   # campo sgn·a·v tem S>0
    for lam in (1, 10, 100):
        print(f"   a={a:.3g}, λ={lam}: (S−θνz)/y² = {lam**3*(Sa - th*nu*za) / (lam*ya)**2:.4g},  S/(νz) = {Sa/(nu*za):.4g}")
