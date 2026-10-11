# Navier–Stokes — medida em solução real (Taylor–Green)

`taylor_green.py` resolve Navier–Stokes 3D periódico (pseudo-espectral, desaliasing 2/3, RK4) a partir do
vórtice de Taylor–Green e mede, ao longo do tempo, as grandezas do balanço de enstrofia usado no M-Core:

y = ‖ω‖², z = ‖∇ω‖², S = ∫ ω·(∇u)ω,  com  y'/2 + ν z = S  (conferido numericamente, resíduo ≲ 5·10⁻⁴).

| Re | max S/(νz) | % do tempo com S > νz | C_min(θ=0.7) |
|---|---|---|---|
| 100 | 2.776 (t = 1.40) | 45.3% | 5.48e-4 |
| 200 | 5.672 (t = 1.35) | 58.2% | 6.42e-4 |
| 400 | 11.458 (t = 1.35) | 75.1% | 6.88e-4 |

(64³; Re = 400 fica sub-resolvido perto do pico de enstrofia, mas o máximo da razão ocorre cedo, em t ≈ 1.35.)

**Consequência para o G1 do M-Core:** `BulkAbsorption ∧ ExceptionalTailControl ∧ SubcriticalBudget`, com
bulk + tail = S, implicam S ≤ (c_B + δ) ν z < ν z, logo y' < 0 sempre que z > 0. Esta solução real tem
S > νz em boa parte do tempo, então essas hipóteses são **falsas** para ela, qualquer que seja a divisão
bulk/tail. Um resto R com R ≥ S − θνz é inevitável; C_min(θ) é a menor constante em R ≤ C y² compatível
com estes dados (normalização: integrais sobre o cubo [0, 2π]³).
