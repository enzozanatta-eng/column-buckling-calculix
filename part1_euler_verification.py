"""
Part 1: Linear buckling of a column, verified against Euler
=============================================================
Gmsh (hexahedral mesh) + CalculiX (*BUCKLE) + Python

The same steel column is solved with the four classical end conditions.
Each critical load is compared with Euler's formula, each buckled shape
with the analytical mode shape, and the result is checked for mesh
independence.

Units: mm, N, MPa
"""
import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gmsh

# ============================ SETTINGS ============================
CCX = r"C:\Users\enzoz\Downloads\PrePoMax\PrePoMax v3.0.0\Solver\ccx_dynamic.exe"

L, B, H = 400.0, 20.0, 10.0      # length (x), width (y), thickness (z) [mm]
E, NU = 210000.0, 0.3            # steel [MPa]
MESH = (40, 4, 4)                # elements in x, y, z (y and z must be even)
OUT = "results_part1"
# ==================================================================

I_WEAK = B * H ** 3 / 12         # buckling happens about the weak axis
P_EULER = np.pi ** 2 * E * I_WEAK / L ** 2
os.makedirs(OUT, exist_ok=True)

# Effective length factor K and analytical buckled shape for each case
# (x measured from the bottom, shapes normalized to a maximum of 1)
KL_FIXED_PINNED = 4.49341        # root of tan(kL) = kL


def shape_fixed_free(x):
    return 1 - np.cos(np.pi * x / (2 * L))


def shape_pinned_pinned(x):
    return np.sin(np.pi * x / L)


def shape_fixed_pinned(x):
    k = KL_FIXED_PINNED / L
    return -np.sin(k * x) / k + L * np.cos(k * x) + x - L


def shape_fixed_fixed(x):
    return 1 - np.cos(2 * np.pi * x / L)


CASES = [
    # name, bottom, top, K, analytical shape
    ("Fixed-free", "fixed", "free", 2.0, shape_fixed_free),
    ("Pinned-pinned", "pinned", "pinned", 1.0, shape_pinned_pinned),
    ("Fixed-pinned", "fixed", "pinned", np.pi / KL_FIXED_PINNED, shape_fixed_pinned),
    ("Fixed-fixed", "fixed", "guided", 0.5, shape_fixed_fixed),
]


# ------------------------------------------------------------------
# MESH
# ------------------------------------------------------------------
def build_mesh(nx, ny, nz):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("column")
    p = gmsh.model.geo.addPoint(0, 0, 0)
    line = gmsh.model.geo.extrude([(0, p)], L, 0, 0, [nx])
    surf = gmsh.model.geo.extrude([line[1]], 0, B, 0, [ny], recombine=True)
    gmsh.model.geo.extrude([surf[1]], 0, 0, H, [nz], recombine=True)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = {int(t): np.array(coords[3 * i:3 * i + 3]) for i, t in enumerate(tags)}
    _, etags, enodes = gmsh.model.mesh.getElements(3)
    elements = []
    for i, t in enumerate(etags[0]):
        conn = [int(n) for n in enodes[0][8 * i:8 * i + 8]]
        p0, p1, p3, p4 = (nodes[conn[j]] for j in (0, 1, 3, 4))
        if np.dot(np.cross(p1 - p0, p3 - p0), p4 - p0) < 0:
            conn = [conn[0], conn[3], conn[2], conn[1], conn[4], conn[7], conn[6], conn[5]]
        elements.append((int(t), conn))
    gmsh.finalize()
    return nodes, elements


# ------------------------------------------------------------------
# INPUT DECK
#   fixed : every node of the end face fixed in x, y, z
#   pinned: end face fixed in y and z, only the centroid fixed in x,
#           so the face can rotate freely (a hinge)
#   guided: end face fixed in y and z, and every node forced to move
#           together in x (*EQUATION), so the face cannot rotate but
#           can still shorten under the load
#   free  : no restraint
# ------------------------------------------------------------------
def write_inp(nodes, elements, bottom, top, path, n_modes=3):
    tol = 1e-6
    bot = [n for n, c in nodes.items() if abs(c[0]) < tol]
    tip = [n for n, c in nodes.items() if abs(c[0] - L) < tol]

    def centroid(face):
        return min(face, key=lambda n: (nodes[n][1] - B / 2) ** 2 + (nodes[n][2] - H / 2) ** 2)

    with open(path, "w") as f:
        f.write("*NODE, NSET=NALL\n")
        for n, c in nodes.items():
            f.write(f"{n}, {c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}\n")
        f.write("*ELEMENT, TYPE=C3D8I, ELSET=EALL\n")
        for e, conn in elements:
            f.write(f"{e}, " + ", ".join(map(str, conn)) + "\n")
        for name, face in (("BOT", bot), ("TIP", tip)):
            f.write(f"*NSET, NSET={name}\n")
            for i in range(0, len(face), 10):
                f.write(", ".join(map(str, face[i:i + 10])) + "\n")
        f.write(f"*MATERIAL, NAME=STEEL\n*ELASTIC\n{E}, {NU}\n")
        f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL\n")

        if top == "guided":
            ref = centroid(tip)
            f.write("*EQUATION\n")
            for n in tip:
                if n != ref:
                    f.write(f"2\n{n}, 1, 1.0, {ref}, 1, -1.0\n")

        f.write("*BOUNDARY\n")
        if bottom == "fixed":
            f.write("BOT, 1, 3\n")
        else:
            f.write(f"BOT, 2, 3\n{centroid(bot)}, 1, 1\n")
        if top in ("pinned", "guided"):
            f.write("TIP, 2, 3\n")

        f.write(f"*STEP\n*BUCKLE\n{n_modes}\n*CLOAD\n")
        if top == "guided":
            f.write(f"{centroid(tip)}, 1, -1.0\n")            # unit load, all on the guided face
        else:
            for n in tip:
                f.write(f"{n}, 1, {-1.0 / len(tip):.10f}\n")   # unit load spread over the face
        f.write("*NODE FILE\nU\n*END STEP\n")


def run_ccx(job):
    res = subprocess.run([CCX, job], capture_output=True, text=True)
    if not os.path.exists(job + ".dat"):
        print(res.stdout[-2000:])
        raise RuntimeError("CalculiX failed")


def read_buckling_factors(path):
    factors, reading = [], False
    with open(path) as f:
        for line in f:
            if "B U C K L I N G   F A C T O R" in line:
                reading = True
                continue
            parts = line.split()
            if reading and len(parts) == 2 and parts[0].isdigit():
                factors.append(float(parts[1]))
    return np.array(factors)


def read_first_mode(path):
    """Displacements of buckling mode 1 (the first DISP block is the static pre-load step)."""
    blocks, current = [], None
    with open(path) as f:
        for line in f:
            if line.startswith(" -4  DISP"):
                current = {}
            elif current is not None and line.startswith(" -1"):
                current[int(line[3:13])] = np.array(
                    [float(line[13:25]), float(line[25:37]), float(line[37:49])])
            elif current is not None and line.startswith(" -3"):
                blocks.append(current)
                current = None
    return blocks[1]


def solve(bottom, top, mesh, job):
    nodes, elements = build_mesh(*mesh)
    write_inp(nodes, elements, bottom, top, job + ".inp")
    run_ccx(job)
    return nodes, elements, read_buckling_factors(job + ".dat"), read_first_mode(job + ".frd")


def centerline(nodes, mode):
    """Mean lateral (z) displacement of each cross-section, normalized to max |value| = 1."""
    xs = sorted({round(c[0], 6) for c in nodes.values()})
    w = {x: [] for x in xs}
    for n, c in nodes.items():
        w[round(c[0], 6)].append(mode[n][2])
    w = np.array([np.mean(w[x]) for x in xs])
    return np.array(xs), w / w[np.argmax(np.abs(w))]


# ==================================================================
def main():
    report = []

    def log(text=""):
        print(text)
        report.append(text)

    log(f"Steel column {L:.0f} x {B:.0f} x {H:.0f} mm, buckling about the weak axis "
        f"(I = {I_WEAK:.0f} mm^4) | mesh {MESH}\n")
    log(f"{'End conditions':<15} {'K':>6} {'Euler [N]':>11} {'FEA [N]':>11} {'Diff':>8} {'MAC':>8}")

    results = []
    for name, bottom, top, K, shape in CASES:
        job = "buckle_" + name.lower().replace("-", "_")
        nodes, elements, factors, mode = solve(bottom, top, MESH, job)
        p_th = P_EULER / K ** 2
        x, w = centerline(nodes, mode)
        phi = shape(x)
        phi = phi / phi[np.argmax(np.abs(phi))]
        mac = np.dot(w, phi) ** 2 / (np.dot(w, w) * np.dot(phi, phi))
        diff = (factors[0] - p_th) / p_th * 100
        results.append(dict(name=name, K=K, p_th=p_th, p=factors[0], x=x, w=w, shape=shape, mac=mac))
        log(f"{name:<15} {K:>6.3f} {p_th:>11.0f} {factors[0]:>11.0f} {diff:>+7.2f}% {mac:>8.5f}")

    # ---------------- Buckled shapes ----------------
    fig, axes = plt.subplots(1, 4, figsize=(12, 6), sharey=True)
    xf = np.linspace(0, L, 300)
    for ax, r in zip(axes, results):
        phi = r["shape"](xf)
        phi = phi / phi[np.argmax(np.abs(phi))]
        ax.plot([0, 0], [0, L], color="gray", lw=6, alpha=0.25, solid_capstyle="butt")
        ax.plot(phi, xf, color="tab:blue", lw=2, label="Euler theory")
        ax.plot(r["w"][::2], r["x"][::2], "o", color="tab:red", ms=4, label="FEA (CalculiX)")
        ax.set_title(f"{r['name']}\nK = {r['K']:.3f}\n"
                     f"P = {r['p'] / 1e3:.2f} kN (theory {r['p_th'] / 1e3:.2f})", fontsize=9)
        ax.set_xlim(-1.3, 1.3)
        ax.set_xlabel("Normalized lateral deflection")
        ax.grid(alpha=0.3)
        for y, kind in ((0, BOTTOM_SYMBOL[r["name"]][0]), (L, BOTTOM_SYMBOL[r["name"]][1])):
            ax.text(0, y + (-18 if y == 0 else 8), kind, ha="center", fontsize=8, color="dimgray")
    axes[0].set_ylabel("Height along the column [mm]")
    axes[0].legend(fontsize=8, loc="center left")
    fig.suptitle("First buckling mode for the four classical end conditions", fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{OUT}/buckling_modes.png", dpi=150)
    plt.close(fig)

    # ---------------- Mesh convergence (pinned-pinned) ----------------
    log("\nMesh convergence, pinned-pinned case:")
    meshes = [(10, 2, 2), (20, 2, 2), (40, 4, 4), (80, 4, 4), (160, 8, 8)]
    conv = []
    for m in meshes:
        _, el, fac, _ = solve("pinned", "pinned", m, f"conv_{m[0]}x{m[1]}x{m[2]}")
        conv.append((len(el), fac[0]))
        log(f"  {len(el):>6} elements: P_cr = {fac[0]:9.1f} N  ({(fac[0] - P_EULER) / P_EULER * 100:+.3f}% vs Euler)")

    fig, ax = plt.subplots(figsize=(7, 4.3))
    ax.semilogx([c[0] for c in conv], [c[1] / 1e3 for c in conv], "o-", color="tab:red", label="FEA")
    ax.axhline(P_EULER / 1e3, color="tab:blue", ls="--", label=f"Euler = {P_EULER / 1e3:.3f} kN")
    ax.set_xlabel("Number of elements")
    ax.set_ylabel("Critical load [kN]")
    ax.set_title("Mesh convergence of the critical load (pinned-pinned)")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{OUT}/mesh_convergence.png", dpi=150)
    plt.close(fig)

    with open(f"{OUT}/report.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(report))
    print(f"\nFigures and report saved in the '{OUT}' folder.")


BOTTOM_SYMBOL = {
    "Fixed-free": ("fixed", "free"),
    "Pinned-pinned": ("pinned", "pinned"),
    "Fixed-pinned": ("fixed", "pinned"),
    "Fixed-fixed": ("fixed", "fixed"),
}

if __name__ == "__main__":
    main()
