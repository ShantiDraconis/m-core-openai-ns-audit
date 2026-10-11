"""
ATK — Navier–Stokes: medida da razão stretch/(ν z) numa solução real (vórtice de Taylor–Green).

Equação incompressível em [0, 2π]³ periódico, sem força:
    ∂t u + (u·∇)u + ∇p = ν Δu,   ∇·u = 0,
    u(0) = (sin x cos y cos z, −cos x sin y cos z, 0).

Grandezas do balanço de enstrofia usado no M-Core (EnstrophyBalance):
    y = ‖ω‖²,  z = ‖∇ω‖²,  S = ∫ ω·(∇u)ω dx,    y'/2 + ν z = S.

As hipóteses BulkAbsorption + ExceptionalTailControl + SubcriticalBudget
(bulk + tail = S, bulk ≤ c_B ν z, tail ≤ δ ν z, c_B + δ < 1) implicam S < ν z sempre que z > 0,
isto é, y' < 0. Se a solução real tiver S/(ν z) > 1 em algum instante, essas hipóteses são FALSAS
para ela, qualquer que seja a divisão bulk/tail escolhida.

Método: pseudo-espectral (rfft), desaliasing 2/3, RK4. Confere o próprio balanço numericamente.
Uso: python3 taylor_green.py N Re T dt saida.csv
"""
import sys, time
import numpy as np

N, Re, T, dt, out = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
nu = 1.0 / Re
L = 2 * np.pi
x = np.arange(N) * L / N
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

k1 = np.fft.fftfreq(N, 1.0 / N)
kz1 = np.fft.rfftfreq(N, 1.0 / N)
KX, KY, KZ = np.meshgrid(k1, k1, kz1, indexing="ij")
K2 = KX**2 + KY**2 + KZ**2
K2safe = np.where(K2 == 0, 1.0, K2)
kmax = N // 3
dealias = (np.abs(KX) <= kmax) & (np.abs(KY) <= kmax) & (np.abs(KZ) <= kmax)
Kvec = (KX, KY, KZ)

fft = lambda f: np.fft.rfftn(f)
ifft = lambda F: np.fft.irfftn(F, s=(N, N, N), axes=(0, 1, 2))


def project(Uh):
    div = (KX * Uh[0] + KY * Uh[1] + KZ * Uh[2]) / K2safe
    return np.array([Uh[i] - Kvec[i] * div for i in range(3)])


def curl_hat(Uh):
    return np.array([1j * (KY * Uh[2] - KZ * Uh[1]),
                     1j * (KZ * Uh[0] - KX * Uh[2]),
                     1j * (KX * Uh[1] - KY * Uh[0])])


def rhs(Uh):
    u = np.array([ifft(Uh[i]) for i in range(3)])
    w = np.array([ifft(c) for c in curl_hat(Uh)])
    nl = np.array([u[1] * w[2] - u[2] * w[1],
                   u[2] * w[0] - u[0] * w[2],
                   u[0] * w[1] - u[1] * w[0]])            # u × ω
    NLh = np.array([fft(nl[i]) for i in range(3)]) * dealias
    return project(NLh) - nu * K2 * Uh


def norm2_hat(Fh):
    """∫ |f|² dx a partir de rfft (Parseval, contando a metade espelhada)."""
    w = np.full(Fh.shape, 2.0)
    w[..., 0] = 1.0
    if N % 2 == 0:
        w[..., -1] = 1.0
    return float(np.sum(w * np.abs(Fh) ** 2)) * (L**3) / N**6


def diagnostics(Uh):
    Wh = curl_hat(Uh)
    y = sum(norm2_hat(Wh[i]) for i in range(3))
    z = sum(norm2_hat(1j * Kvec[j] * Wh[i]) for i in range(3) for j in range(3))
    w = np.array([ifft(Wh[i]) for i in range(3)])
    S = 0.0
    for i in range(3):
        for j in range(3):
            dju_i = ifft(1j * Kvec[j] * Uh[i])
            S += np.mean(w[i] * dju_i * w[j])
    S *= L**3
    E = 0.5 * sum(norm2_hat(Uh[i]) for i in range(3))
    return E, y, z, S


u0 = np.array([np.sin(X) * np.cos(Y) * np.cos(Z), -np.cos(X) * np.sin(Y) * np.cos(Z), 0 * X])
Uh = project(np.array([fft(u0[i]) for i in range(3)]) * dealias)

rows, t, step = [], 0.0, 0
every = max(1, int(round(0.05 / dt)))
t0 = time.time()
nsteps = int(round(T / dt))
for step in range(nsteps + 1):
    if step % every == 0:
        E, y, z, S = diagnostics(Uh)
        rows.append((t, E, y, z, S, S / (nu * z)))
    if step == nsteps:
        break
    k1_ = rhs(Uh)
    k2_ = rhs(Uh + 0.5 * dt * k1_)
    k3_ = rhs(Uh + 0.5 * dt * k2_)
    k4_ = rhs(Uh + dt * k3_)
    Uh = Uh + dt / 6 * (k1_ + 2 * k2_ + 2 * k3_ + k4_)
    t += dt

arr = np.array(rows)
# verificação do balanço: y'/2 (diferença central) vs S - ν z
tt, yy, zz, SS = arr[:, 0], arr[:, 2], arr[:, 3], arr[:, 4]
dy = np.gradient(yy, tt)
resid = np.max(np.abs(dy[2:-2] / 2 - (SS[2:-2] - nu * zz[2:-2]))) / np.max(np.abs(SS))
with open(out, "w") as f:
    f.write("t,E,y_enstrofia,z_palinstrofia,S_stretch,razao_S_sobre_nu_z\n")
    for r in rows:
        f.write(",".join(f"{v:.10g}" for v in r) + "\n")
imax = int(np.argmax(arr[:, 5]))
print(f"N={N} Re={Re:g} nu={nu:g} dt={dt:g} T={T:g} tempo={time.time()-t0:.0f}s")
print(f"max S/(nu z) = {arr[imax, 5]:.6f} em t = {arr[imax, 0]:.2f}")
print(f"S/(nu z) > 1 em {np.mean(arr[:, 5] > 1) * 100:.1f}% dos instantes amostrados")
print(f"pico de enstrofia y em t = {arr[int(np.argmax(yy)), 0]:.2f}  (y(0) = {yy[0]:.6g}, y_max = {yy.max():.6g})")
print(f"energia: E(0) = {arr[0,1]:.6f}  E(T) = {arr[-1,1]:.6f}")
print(f"residuo relativo do balanco y'/2 = S - nu z: {resid:.2e}")
print(f"resolucao: k_max*eta_min ~ {kmax * (nu**3 / (nu*arr[:,2].max()/L**3))**0.25:.2f}")
