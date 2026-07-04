import numpy as np
import matplotlib.pyplot as plt
import os, time

N = 1000
L = 25.0
epsilon = 1.0
sigma = 1.0
r_cut = 2.8 * sigma
r_cut2 = r_cut * r_cut

max_iter = 3000
force_tol = 1e-5
alpha0 = 0.02
random_seed = 42

nl_rebuild_freq = 8

np.random.seed(random_seed)
os.makedirs('results', exist_ok=True)

def minimum_image(vec, box):
    return vec - box * np.round(vec / box)

def build_neighbor_list(pos, box, rcut):
    n = pos.shape[0]
    ncell = max(1, int(box // rcut))
    cell_size = box / ncell
    cells = {}
    for i in range(n):
        cx = int(np.floor(pos[i,0] / cell_size)) % ncell
        cy = int(np.floor(pos[i,1] / cell_size)) % ncell
        cells.setdefault((cx,cy), []).append(i)
    neighbors = [[] for _ in range(n)]
    for cx in range(ncell):
        for cy in range(ncell):
            bucket = cells.get((cx,cy), [])
            for dx in (-1,0,1):
                for dy in (-1,0,1):
                    ncx = (cx+dx) % ncell
                    ncy = (cy+dy) % ncell
                    neigh_bucket = cells.get((ncx,ncy), [])
                    for i in bucket:
                        for j in neigh_bucket:
                            if j <= i:
                                continue
                            rij = minimum_image(pos[j]-pos[i], box)
                            if rij[0]*rij[0] + rij[1]*rij[1] <= rcut*rcut:
                                neighbors[i].append(j)
    return neighbors

def compute_E_F(pos, neighbors, box):
    n = pos.shape[0]
    F = np.zeros_like(pos)
    E = 0.0
    for i in range(n):
        for j in neighbors[i]:
            rvec = minimum_image(pos[j]-pos[i], box)
            r2 = rvec[0]*rvec[0] + rvec[1]*rvec[1]
            if r2 == 0.0:
                continue
            inv_r2 = (sigma*sigma) / r2
            inv_r6 = inv_r2**3
            inv_r12 = inv_r6**2
            E += 4.0 * epsilon * (inv_r12 - inv_r6)
            f_scalar = 24.0 * epsilon * (2.0*inv_r12 - inv_r6) / r2
            fvec = f_scalar * rvec
            F[i] += fvec
            F[j] -= fvec
    return E, F

def init_positions(n, box):
    side = int(np.ceil(np.sqrt(n)))
    margin = 0.3
    xs = np.linspace(margin, box - margin, side)
    ys = np.linspace(margin, box - margin, side)
    grid = np.array(np.meshgrid(xs, ys)).T.reshape(-1,2)
    pos = grid[:n].copy()
    spacing = xs[1] - xs[0] if side > 1 else 1.0
    jitter = 0.2 * spacing
    pos += jitter * (np.random.rand(n,2) - 0.5)
    return pos

def plot_positions(pos, fname, box, title='', dot_size=6):
    plt.figure(figsize=(5,5))
    plt.scatter(pos[:,0], pos[:,1], s=dot_size)
    plt.xlim(0, box); plt.ylim(0, box)
    plt.gca().set_aspect('equal','box')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(fname, dpi=150); plt.close()

def plot_overlay(initial, final, fname, box, arrow_scale=1.0, max_arrows=1000):
    n = initial.shape[0]
    idx = np.arange(n)
    if n > max_arrows:
        idx = np.random.choice(n, max_arrows, replace=False)
    disp = final - initial
    for k in idx:
        disp[k] = minimum_image(disp[k], box)
    lengths = np.linalg.norm(disp[idx], axis=1)
    nonzero = lengths > 0
    plt.figure(figsize=(6,6))
    plt.scatter(initial[idx,0], initial[idx,1], s=8, alpha=0.6, label='initial')
    plt.quiver(initial[idx,0], initial[idx,1],
               disp[idx,0]*arrow_scale, disp[idx,1]*arrow_scale,
               angles='xy', scale_units='xy', scale=1, width=0.002, alpha=0.8, color='r', label='displacement')
    plt.xlim(0, box); plt.ylim(0, box)
    plt.gca().set_aspect('equal','box')
    plt.legend()
    plt.tight_layout()
    plt.savefig(fname, dpi=200); plt.close()

pos0 = init_positions(N, L)
np.savetxt('results/initial_positions.txt', pos0)
plot_positions(pos0, 'results/initial.png', L, title='Initial positions')

mpd = 1e12
if N <= 2000:
    for i in range(N-1):
        for j in range(i+1,N):
            d2 = np.sum(minimum_image(pos0[j]-pos0[i], L)**2)
            if d2 < mpd:
                mpd = d2
else:
    idx = np.random.choice(N, min(2000, N), replace=False)
    for i in range(len(idx)-1):
        for j in range(i+1, len(idx)):
            d2 = np.sum(minimum_image(pos0[idx[j]]-pos0[idx[i]], L)**2)
            if d2 < mpd:
                mpd = d2
mpd = np.sqrt(mpd)
print("Initial min pair distance (approx) =", mpd)

neighbors = build_neighbor_list(pos0, L, r_cut)
E, F = compute_E_F(pos0, neighbors, L)
alpha = alpha0

energies = [E]
max_forces = [np.max(np.linalg.norm(F, axis=1))]

pos = pos0.copy()
start = time.time()

for it in range(1, max_iter+1):
    trial = pos - alpha * F
    if it % nl_rebuild_freq == 0:
        neighbors_trial = build_neighbor_list(trial, L, r_cut)
    else:
        neighbors_trial = neighbors
    E_trial, F_trial = compute_E_F(trial, neighbors_trial, L)
    bt = 0
    while E_trial > E and alpha > 1e-12:
        alpha *= 0.5
        trial = pos - alpha * F
        if it % nl_rebuild_freq == 0:
            neighbors_trial = build_neighbor_list(trial, L, r_cut)
        E_trial, F_trial = compute_E_F(trial, neighbors_trial, L)
        bt += 1
        if bt > 30:
            break
    pos = trial
    neighbors = neighbors_trial
    E = E_trial
    F = F_trial
    energies.append(E)
    max_forces.append(np.max(np.linalg.norm(F, axis=1)))
    if bt == 0 and alpha < alpha0:
        alpha = min(alpha0, alpha * 1.05)
    if it % 50 == 0:
        elapsed = time.time() - start
        print(f"it {it:5d} E {E:.6f} maxF {max_forces[-1]:.3e} alpha {alpha:.1e} elapsed {elapsed:.1f}s")
    if max_forces[-1] < force_tol:
        print(f"Converged at iter {it}, maxF {max_forces[-1]:.3e}, E {E:.6f}")
        break

end = time.time()
print("Total time (s):", end - start)

np.savetxt('results/final_positions.txt', pos)
plot_positions(pos, 'results/final.png', L, title='Final positions')
plt.figure(); plt.plot(energies); plt.xlabel('Iteration'); plt.ylabel('Energy'); plt.tight_layout(); plt.savefig('results/energy.png'); plt.close()
plt.figure(); plt.semilogy(max_forces); plt.xlabel('Iteration'); plt.ylabel('Max Force'); plt.tight_layout(); plt.savefig('results/maxforce.png'); plt.close()

disp = pos - pos0
for i in range(N):
    disp[i] = minimum_image(disp[i], L)
dists = np.linalg.norm(disp, axis=1)
rms = np.sqrt(np.mean(dists**2))
mmax = np.max(dists)
count_moved = np.sum(dists > 1e-6)
print(f"RMS displacement = {rms:.6f}, max displacement = {mmax:.6f}, particles moved = {count_moved}/{N}")

plot_overlay(pos0, pos, 'results/overlay.png', L, arrow_scale=1.0, max_arrows=1500)

print("Saved: initial/final positions + initial.png, final.png, overlay.png, energy.png, maxforce.png")