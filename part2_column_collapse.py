"""
Part 2: Collapse of a steel mezzanine column (geometric + material nonlinear FEA)
===================================================================================
Gmsh (shell mesh) + CalculiX + Python, checked against Eurocode 3 (EN 1993-1-1)

Real-world case
---------------
A storage mezzanine is supported by hot-finished square hollow section columns,
SHS 80x80x4 in S355 steel, 3.0 m high, pinned at both ends. Each column carries
a tributary floor area of 3.0 m x 3.5 m. Will it hold, and how much margin is there?

What the script does
--------------------
1. Hand calculation: section class, Euler load, Eurocode 3 flexural buckling resistance.
2. Linear buckling (*BUCKLE) of the shell model, compared with Euler.
3. Nonlinear collapse analysis (GMNIA): large displacements (NLGEOM), elastic-plastic
   steel, initial bow imperfection of L/1000, loaded by imposed shortening so the
   analysis follows the load past the peak.
4. Slenderness study: the full analysis repeated for lengths from 0.75 m to 6 m and
   plotted on top of the Eurocode 3 buckling curves.
5. Imperfection sensitivity and mesh check.
6. Figures: load-shortening curves, 3D collapse shape with plastic strain, design charts.

Units: mm, N, MPa
"""
import os
import re
import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import gmsh

# ============================ SETTINGS ============================
CCX = r"C:\Users\enzoz\Downloads\PrePoMax\PrePoMax v3.0.0\Solver\ccx_dynamic.exe"
N_JOBS = max(1, (os.cpu_count() or 2) - 1)   # CalculiX runs in parallel

# Section: SHS 80x80x4 (outer width x wall thickness), modeled at the wall midline
WIDTH, T = 80.0, 4.0
A_MID = WIDTH - T                 # midline width = 76 mm
E, NU, FY = 210000.0, 0.3, 355.0  # S355 steel
GAMMA_M1 = 1.0                    # partial factor for member stability (EN 1993-1-1)
ALPHA = 0.21                      # imperfection factor, curve "a" (hot-finished SHS, S355)

# Real-world design case
L_DESIGN = 3000.0                 # column height [mm]
TRIB_AREA = 3.0 * 3.5             # tributary floor area [m^2]
G_K = 1.5                         # permanent load: deck, joists, finishes [kN/m^2]
Q_K = 7.5                         # imposed load: storage, EN 1991-1-1 category E1 [kN/m^2]
N_ED = (1.35 * G_K + 1.5 * Q_K) * TRIB_AREA * 1e3   # design axial force [N]

# Mesh: elements per side of the tube and element length along the column
NS = 4
ELEM_LEN = 30.0
OUT = "results_part2"
# ==================================================================

os.makedirs(OUT, exist_ok=True)

# Section properties (thin-walled, midline, sharp corners)
AREA = 4 * A_MID * T
I_SEC = 2 * (T * A_MID ** 3 / 12) + 2 * (A_MID * T * (A_MID / 2) ** 2)
R_GYR = np.sqrt(I_SEC / AREA)
N_PL = AREA * FY
EPS = np.sqrt(235 / FY)
LAMBDA_1 = np.pi * np.sqrt(E / FY)


def euler_load(length):
    return np.pi ** 2 * E * I_SEC / length ** 2


def slenderness(length):
    return length / R_GYR / LAMBDA_1


def chi_ec3(lam, alpha=ALPHA):
    """Eurocode 3 reduction factor for flexural buckling (EN 1993-1-1, 6.3.1.2)."""
    lam = np.asarray(lam, dtype=float)
    phi = 0.5 * (1 + alpha * (lam - 0.2) + lam ** 2)
    chi = 1 / (phi + np.sqrt(phi ** 2 - lam ** 2))
    return np.minimum(chi, 1.0)


# ------------------------------------------------------------------
# MESH: four walls (shell quads) + a stiff end plate at each end
# ------------------------------------------------------------------
def build_mesh(length, ns=NS, elem_len=ELEM_LEN):
    a = A_MID
    nl = max(4, int(round(length / elem_len)))
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("shs")
    g = gmsh.model.geo
    pts = [g.addPoint(0, y, z) for y, z in [(-a / 2, -a / 2), (a / 2, -a / 2), (a / 2, a / 2), (-a / 2, a / 2)]]
    base = [g.addLine(pts[i], pts[(i + 1) % 4]) for i in range(4)]
    for ln in base:
        g.mesh.setTransfiniteCurve(ln, ns + 1)
    top = [g.extrude([(1, ln)], length, 0, 0, [nl], recombine=True)[0][1] for ln in base]
    for loop in (base, top):
        s = g.addPlaneSurface([g.addCurveLoop(loop)])
        g.mesh.setTransfiniteSurface(s)
        g.mesh.setRecombine(2, s)
    g.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.removeDuplicateNodes()
    tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = {int(t): np.array(coords[3 * i:3 * i + 3]) for i, t in enumerate(tags)}
    _, etags, enodes = gmsh.model.mesh.getElements(2)
    gmsh.finalize()

    walls, plates = [], []
    for i, t in enumerate(etags[0]):
        conn = [int(n) for n in enodes[0][4 * i:4 * i + 4]]
        P = [nodes[n] for n in conn]
        normal = np.cross(P[1] - P[0], P[3] - P[0])
        center = np.mean(P, axis=0)
        if np.ptp([p[0] for p in P]) < 1e-6:            # end plate (lies in a y-z plane)
            if normal[0] < 0:
                conn = conn[::-1]
            plates.append((int(t), conn))
        else:                                           # wall: make every normal point outward
            if np.dot(normal, [0, center[1], center[2]]) < 0:
                conn = conn[::-1]
            walls.append((int(t), conn))

    def find(x, y, z):
        return min(nodes, key=lambda n: np.sum((nodes[n] - [x, y, z]) ** 2))

    special = dict(bot=find(0, 0, 0), top=find(length, 0, 0),
                   twist=find(0, 0, a / 2), mid=find(length / 2, 0, a / 2))
    return nodes, walls, plates, special


# ------------------------------------------------------------------
# INPUT DECK
#   Pinned ends: the centre of each end plate is a hinge. Bottom centre fixed
#   in x, y, z; top centre fixed in y, z and pushed down in x. One extra node
#   stops the column from spinning about its own axis.
# ------------------------------------------------------------------
def write_inp(path, length, analysis, e0=0.0, ns=NS, elem_len=ELEM_LEN, field_output=False):
    nodes, walls, plates, sp = build_mesh(length, ns, elem_len)
    # initial bow imperfection: half sine wave with amplitude e0 at mid-height
    for n in nodes:
        nodes[n] = nodes[n] + [0, 0, e0 * np.sin(np.pi * nodes[n][0] / length)]
    with open(path, "w") as f:
        f.write("*NODE, NSET=NALL\n")
        for n, c in nodes.items():
            f.write(f"{n}, {c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}\n")
        for name, els in (("WALLS", walls), ("PLATES", plates)):
            f.write(f"*ELEMENT, TYPE=S4, ELSET={name}\n")
            for e, conn in els:
                f.write(f"{e}, " + ", ".join(map(str, conn)) + "\n")
        f.write(f"*NSET, NSET=CTOP\n{sp['top']}\n*NSET, NSET=MID\n{sp['mid']}\n")
        f.write(f"*MATERIAL, NAME=S355\n*ELASTIC\n{E}, {NU}\n")
        if analysis == "collapse":
            # elastic-plastic with a very small hardening slope (E/10000)
            f.write(f"*PLASTIC\n{FY}, 0.0\n{FY + 0.1 * E / 10000:.2f}, 0.1\n")
        f.write(f"*MATERIAL, NAME=RIGID\n*ELASTIC\n{10 * E}, {NU}\n")
        f.write(f"*SHELL SECTION, ELSET=WALLS, MATERIAL=S355\n{T}\n")
        f.write("*SHELL SECTION, ELSET=PLATES, MATERIAL=RIGID\n40\n")
        f.write(f"*BOUNDARY\n{sp['bot']}, 1, 3\n{sp['top']}, 2, 3\n{sp['twist']}, 2, 2\n")
        if analysis == "buckle":
            f.write(f"*STEP\n*BUCKLE\n2\n*CLOAD\n{sp['top']}, 1, -1.0\n*END STEP\n")
        else:
            shortening = 1.6 * length * FY / E
            f.write("*STEP, NLGEOM, INC=1000\n*STATIC\n0.02, 1.0, 1e-6, 0.02\n")
            f.write(f"*BOUNDARY\n{sp['top']}, 1, 1, {-shortening:.6f}\n")
            f.write("*NODE PRINT, NSET=CTOP\nRF, U\n*NODE PRINT, NSET=MID\nU\n")
            if field_output:
                f.write("*NODE FILE\nU\n*EL FILE\nS, PEEQ\n")
            f.write("*END STEP\n")


def run_ccx(job):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    res = subprocess.run([CCX, job], capture_output=True, text=True, env=env)
    return res.returncode == 0


def run_all(jobs):
    """Run several CalculiX jobs at the same time."""
    with ThreadPoolExecutor(max_workers=N_JOBS) as pool:
        return dict(zip(jobs, pool.map(run_ccx, jobs)))


# ------------------------------------------------------------------
# READ RESULTS
# ------------------------------------------------------------------
def read_buckling_factor(path):
    reading = False
    with open(path) as f:
        for line in f:
            if "B U C K L I N G   F A C T O R" in line:
                reading = True
                continue
            parts = line.split()
            if reading and len(parts) == 2 and parts[0] == "1":
                return float(parts[1])
    raise RuntimeError(f"No buckling factor found in {path}")


def read_curve(path):
    """Axial force, shortening and mid-height lateral deflection for every increment."""
    with open(path) as f:
        txt = f.read()
    pat = r"{} \(.*?\) for set {} and time\s+(\S+)\s*\n\s*\n\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)"
    rf = re.findall(pat.format("forces", "CTOP"), txt)
    ut = re.findall(pat.format("displacements", "CTOP"), txt)
    um = re.findall(pat.format("displacements", "MID"), txt)
    n = min(len(rf), len(ut), len(um))
    time = np.array([float(r[0]) for r in rf[:n]])
    force = np.array([-float(r[1]) for r in rf[:n]])
    short = np.array([-float(u[1]) for u in ut[:n]])
    lateral = np.array([float(u[3]) for u in um[:n]])
    return dict(time=np.r_[0, time], force=np.r_[0, force], short=np.r_[0, short], lateral=np.r_[0, lateral])


def read_frd(path, wanted_times):
    """Expanded 3D mesh and fields (U, von Mises, PEEQ) at the requested increments."""
    nodes, elements, frames = {}, [], {}
    block, cur_time, data = None, None, None
    with open(path) as f:
        for line in f:
            if line.startswith("    2C"):
                block = "nodes"
                continue
            if line.startswith("    3C"):
                block = "elements"
                continue
            if line.startswith("  100CL"):
                cur_time = float(line[12:24])
                continue
            if line.startswith(" -4"):
                name = line.split()[1]
                keep = any(abs(cur_time - t) < 1e-7 for t in wanted_times)
                block = name if keep else "skip"
                data = {}
                continue
            if line.startswith(" -3"):
                if block in ("DISP", "STRESS", "PE"):
                    frames.setdefault(round(cur_time, 7), {})[block] = data
                block = None
                continue
            if block == "nodes" and line.startswith(" -1"):
                nodes[int(line[3:13])] = np.array([float(line[13:25]), float(line[25:37]), float(line[37:49])])
            elif block == "elements" and line.startswith(" -2"):
                elements.append([int(v) for v in line.split()[1:9]])
            elif block in ("DISP", "STRESS", "PE") and line.startswith(" -1"):
                n = int(line[3:13])
                vals = [float(line[13 + 12 * k:25 + 12 * k]) for k in range((len(line.rstrip()) - 13) // 12)]
                if block == "STRESS":
                    sx, sy, sz, txy, tyz, tzx = vals[:6]
                    data[n] = np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2)
                                      + 3 * (txy ** 2 + tyz ** 2 + tzx ** 2))
                elif block == "DISP":
                    data[n] = np.array(vals[:3])
                else:
                    data[n] = vals[0]
    return nodes, elements, frames


# ------------------------------------------------------------------
# 3D PLOT
# ------------------------------------------------------------------
FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def boundary_faces(elements):
    count = {}
    for conn in elements:
        for fc in FACES:
            quad = tuple(conn[i] for i in fc)
            key = tuple(sorted(quad))
            if key in count:
                count[key][0] += 1
            else:
                count[key] = [1, quad]
    return [v[1] for v in count.values() if v[0] == 1]


def draw_column(ax, nodes, faces, disp, field, norm, cmap, scale, xlim=None):
    polys, cols = [], []
    for q in faces:
        xs = [nodes[n][0] for n in q]
        if xlim and (min(xs) < xlim[0] or max(xs) > xlim[1]):
            continue
        # plot the column standing up: column axis x -> vertical axis
        pts = [nodes[n] + scale * disp[n] for n in q]
        polys.append([(p[2], p[1], p[0]) for p in pts])
        cols.append(cmap(norm(np.mean([field[n] for n in q]))))
    ax.add_collection3d(Poly3DCollection(polys, facecolors=cols, edgecolors=(0, 0, 0, 0.08), linewidths=0.2))
    return polys


def set_view(ax, polys, half_lateral, half_vertical, box):
    allp = np.array([p for poly in polys for p in poly])
    cx, cy, zc = allp[:, 0].mean(), allp[:, 1].mean(), allp[:, 2].mean()
    ax.set_xlim(cx - half_lateral, cx + half_lateral)
    ax.set_ylim(cy - half_lateral, cy + half_lateral)
    ax.set_zlim(zc - half_vertical, zc + half_vertical)
    ax.set_box_aspect(box)
    ax.view_init(elev=12, azim=-55)
    ax.set_axis_off()


# ==================================================================
def main():
    report = []

    def log(text=""):
        print(text)
        report.append(text)

    lam_d = slenderness(L_DESIGN)
    chi_d = float(chi_ec3(lam_d))
    nb_rd = chi_d * N_PL / GAMMA_M1
    ct = (WIDTH - 3 * T) / T

    log("=" * 74)
    log("MEZZANINE COLUMN CHECK: SHS 80x80x4, S355, L = 3.0 m, pinned-pinned")
    log("=" * 74)
    log(f"Design load: (1.35 x {G_K} + 1.5 x {Q_K}) kN/m2 x {TRIB_AREA:.1f} m2 = N_Ed = {N_ED / 1e3:.1f} kN")
    log(f"Section (midline): A = {AREA:.0f} mm2, I = {I_SEC / 1e4:.1f} cm4, i = {R_GYR:.1f} mm")
    log(f"Section class: c/t = {ct:.1f} <= 33*eps = {33 * EPS:.1f}  ->  Class 1 (no local buckling before yield)")
    log(f"N_pl = {N_PL / 1e3:.1f} kN | N_cr (Euler) = {euler_load(L_DESIGN) / 1e3:.1f} kN | "
        f"lambda_bar = {lam_d:.3f}")
    log(f"Eurocode 3, curve a: chi = {chi_d:.3f}  ->  N_b,Rd = {nb_rd / 1e3:.1f} kN  ->  "
        f"utilization = {N_ED / nb_rd:.2f}")

    # ---------------- Build every input file ----------------
    lengths = [750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000]
    imperf = {2000: "L/2000", 1000: "L/1000", 500: "L/500", 250: "L/250"}
    jobs = []
    for Ln in lengths:
        write_inp(f"buckle_L{Ln}.inp", Ln, "buckle")
        write_inp(f"collapse_L{Ln}.inp", Ln, "collapse", e0=Ln / 1000, field_output=(Ln == L_DESIGN))
        jobs += [f"buckle_L{Ln}", f"collapse_L{Ln}"]
    for k in imperf:
        if k != 1000:
            write_inp(f"imperf_{k}.inp", L_DESIGN, "collapse", e0=L_DESIGN / k)
            jobs.append(f"imperf_{k}")
    write_inp("mesh_fine.inp", L_DESIGN, "collapse", e0=L_DESIGN / 1000, ns=2 * NS, elem_len=ELEM_LEN / 2)
    jobs.append("mesh_fine")

    print(f"\nRunning {len(jobs)} CalculiX analyses, {N_JOBS} at a time. This takes a few minutes...")
    ok = run_all(jobs)
    failed = [j for j, s in ok.items() if not s]
    if failed:
        print("Note: these runs stopped early (normal after collapse, the peak is still captured):", failed)

    # ---------------- Linear buckling vs Euler ----------------
    log("\nLINEAR BUCKLING (*BUCKLE) vs EULER")
    log(f"{'L [mm]':>7} {'lambda':>7} {'Euler [kN]':>11} {'FEA [kN]':>10} {'Diff':>8}")
    rows = []
    for Ln in lengths:
        ncr_fea = read_buckling_factor(f"buckle_L{Ln}.dat")
        ncr = euler_load(Ln)
        log(f"{Ln:>7} {slenderness(Ln):>7.3f} {ncr / 1e3:>11.1f} {ncr_fea / 1e3:>10.1f} "
            f"{(ncr_fea - ncr) / ncr * 100:>+7.2f}%")
        rows.append(dict(L=Ln, lam=slenderness(Ln), ncr=ncr, ncr_fea=ncr_fea))

    # ---------------- Nonlinear collapse: slenderness study ----------------
    log("\nNONLINEAR COLLAPSE (GMNIA, e0 = L/1000) vs EUROCODE 3")
    log(f"{'L [mm]':>7} {'lambda':>7} {'N_u FEA':>9} {'N_b EC3':>9} {'chi FEA':>8} {'chi EC3':>8} {'FEA/EC3':>8}")
    curves = {}
    for r in rows:
        c = read_curve(f"collapse_L{r['L']}.dat")
        curves[r["L"]] = c
        r["nu"] = c["force"].max()
        r["chi_fea"] = r["nu"] / N_PL
        r["chi_ec3"] = float(chi_ec3(r["lam"]))
        log(f"{r['L']:>7} {r['lam']:>7.3f} {r['nu'] / 1e3:>9.1f} {r['chi_ec3'] * N_PL / 1e3:>9.1f} "
            f"{r['chi_fea']:>8.3f} {r['chi_ec3']:>8.3f} {r['chi_fea'] / r['chi_ec3']:>8.3f}")

    with open(f"{OUT}/collapse_results.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["L_mm", "lambda_bar", "Ncr_Euler_kN", "Ncr_FEA_kN", "Nu_FEA_kN", "Nb_EC3_kN", "chi_FEA", "chi_EC3"])
        for r in rows:
            w.writerow([r["L"], f"{r['lam']:.4f}", f"{r['ncr'] / 1e3:.2f}", f"{r['ncr_fea'] / 1e3:.2f}",
                        f"{r['nu'] / 1e3:.2f}", f"{r['chi_ec3'] * N_PL / 1e3:.2f}",
                        f"{r['chi_fea']:.4f}", f"{r['chi_ec3']:.4f}"])

    # ---------------- Design case verdict ----------------
    design = curves[L_DESIGN]
    nu_d = design["force"].max()
    i_pk = int(np.argmax(design["force"]))
    log("\nDESIGN CASE VERDICT (L = 3.0 m)")
    log(f"  FEA collapse load          N_u    = {nu_d / 1e3:.1f} kN")
    log(f"  Eurocode 3 resistance      N_b,Rd = {nb_rd / 1e3:.1f} kN")
    log(f"  Design load                N_Ed   = {N_ED / 1e3:.1f} kN")
    log(f"  Utilization (EC3)  = {N_ED / nb_rd:.2f}  ->  {'OK' if N_ED <= nb_rd else 'NOT OK'}")
    log(f"  Safety factor against FEA collapse = {nu_d / N_ED:.2f}")
    log(f"  At collapse: shortening = {design['short'][i_pk]:.2f} mm, "
        f"mid-height deflection = {design['lateral'][i_pk]:.1f} mm "
        f"(initial bow {L_DESIGN / 1000:.1f} mm)")

    # maximum height for this load, by EC3 and by FEA
    Lgrid = np.linspace(500, 6000, 2000)
    nb_grid = chi_ec3(slenderness(Lgrid)) * N_PL
    l_max_ec3 = Lgrid[nb_grid >= N_ED].max()
    Ls = np.array([r["L"] for r in rows])
    nus = np.array([r["nu"] for r in rows])
    l_max_fea = np.interp(-N_ED, -nus, Ls)
    log(f"  Tallest column that still carries N_Ed: {l_max_ec3 / 1e3:.2f} m (EC3), "
        f"{l_max_fea / 1e3:.2f} m (FEA)")

    # ---------------- Imperfection sensitivity ----------------
    log("\nIMPERFECTION SENSITIVITY (L = 3.0 m)")
    imp_curves = {}
    for k, label in imperf.items():
        c = design if k == 1000 else read_curve(f"imperf_{k}.dat")
        imp_curves[label] = c
        log(f"  e0 = {label:<7} ({L_DESIGN / k:5.1f} mm): N_u = {c['force'].max() / 1e3:6.1f} kN  "
            f"({(c['force'].max() - nu_d) / nu_d * 100:+5.1f}% vs L/1000)")

    # ---------------- Mesh check ----------------
    fine = read_curve("mesh_fine.dat")
    log("\nMESH CHECK (L = 3.0 m)")
    log(f"  {NS} elements per side, {ELEM_LEN:.0f} mm long : N_u = {nu_d / 1e3:.1f} kN")
    log(f"  {2 * NS} elements per side, {ELEM_LEN / 2:.0f} mm long : N_u = {fine['force'].max() / 1e3:.1f} kN  "
        f"({(fine['force'].max() - nu_d) / nu_d * 100:+.2f}%)")

    # ================= FIGURES =================
    # 1. Load-shortening and load-deflection for the design case
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    refs = [(N_PL, "Plastic resistance N_pl", "gray", ":"),
            (euler_load(L_DESIGN), "Euler load N_cr", "tab:blue", "--"),
            (nb_rd, "Eurocode 3 N_b,Rd", "tab:green", "-."),
            (N_ED, "Design load N_Ed", "tab:red", "-")]
    for ax, key, xlabel in ((axes[0], "short", "Axial shortening [mm]"),
                            (axes[1], "lateral", "Lateral deflection at mid-height [mm]")):
        ax.plot(design[key], design["force"] / 1e3, color="k", lw=2, label="FEA (GMNIA)")
        ax.plot(design[key][i_pk], nu_d / 1e3, "o", color="k", ms=8)
        ax.annotate(f"Collapse: {nu_d / 1e3:.0f} kN", (design[key][i_pk], nu_d / 1e3),
                    textcoords="offset points", xytext=(10, 8), fontsize=9)
        for val, lab, col, ls in refs:
            ax.axhline(val / 1e3, color=col, ls=ls, lw=1.3, label=f"{lab} = {val / 1e3:.0f} kN")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Axial load [kN]")
        ax.set_ylim(0, N_PL / 1e3 * 1.08)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8, loc="lower right")
    fig.suptitle("SHS 80x80x4 S355 column, L = 3.0 m: nonlinear response up to and beyond collapse", fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{OUT}/load_displacement.png", dpi=150)
    plt.close(fig)

    # 2. FEA vs Eurocode 3 buckling curves + capacity vs height
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    lam = np.linspace(0.01, 2.8, 400)
    for name, al, col in (("a0", 0.13, "tab:purple"), ("a", 0.21, "tab:green"), ("b", 0.34, "tab:orange"),
                          ("c", 0.49, "tab:brown"), ("d", 0.76, "tab:pink")):
        ax.plot(lam, chi_ec3(lam, al), color=col, lw=2.2 if name == "a" else 1, alpha=1 if name == "a" else 0.6,
                label=f"EC3 curve {name}" + (" (this section)" if name == "a" else ""))
    ax.plot(lam[lam > 0.5], 1 / lam[lam > 0.5] ** 2, "k--", lw=1, label="Euler (elastic)")
    ax.plot([r["lam"] for r in rows], [r["chi_fea"] for r in rows], "ko", ms=7, mfc="w", mew=2,
            label="FEA collapse, e0 = L/1000")
    ax.plot(lam_d, nu_d / N_PL, "o", color="tab:red", ms=9, label="Design case (L = 3.0 m)")
    ax.set_xlim(0, 2.8)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel(r"Non-dimensional slenderness $\bar{\lambda}$")
    ax.set_ylabel(r"Reduction factor $\chi = N_u / N_{pl}$")
    ax.set_title("FEA collapse loads on the Eurocode 3 buckling curves")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.plot(Lgrid / 1e3, nb_grid / 1e3, color="tab:green", lw=2, label="Eurocode 3, N_b,Rd")
    ax.plot(Lgrid / 1e3, np.minimum(euler_load(Lgrid), N_PL) / 1e3, "k--", lw=1, label="min(Euler, N_pl)")
    ax.plot(Ls / 1e3, nus / 1e3, "ko", ms=7, mfc="w", mew=2, label="FEA collapse load")
    ax.axhline(N_ED / 1e3, color="tab:red", lw=1.5, label=f"Design load N_Ed = {N_ED / 1e3:.0f} kN")
    ax.axvline(l_max_ec3 / 1e3, color="tab:green", ls=":", lw=1)
    ax.text(l_max_ec3 / 1e3 + 0.05, N_PL / 1e3 * 0.75, f"max height\n(EC3) = {l_max_ec3 / 1e3:.2f} m",
            fontsize=8, color="tab:green")
    ax.axvline(L_DESIGN / 1e3, color="tab:red", ls=":", lw=1)
    ax.set_xlabel("Column height L [m]")
    ax.set_ylabel("Axial capacity [kN]")
    ax.set_title("How tall can this column be?")
    ax.set_ylim(0, N_PL / 1e3 * 1.1)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{OUT}/buckling_curves.png", dpi=150)
    plt.close(fig)

    # 3. Imperfection sensitivity
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    cmap_i = plt.get_cmap("viridis")
    for j, (label, c) in enumerate(imp_curves.items()):
        col = cmap_i(j / (len(imp_curves) - 1))
        axes[0].plot(c["lateral"], c["force"] / 1e3, color=col, lw=2, label=f"e0 = {label}")
        axes[1].plot(c["short"], c["force"] / 1e3, color=col, lw=2, label=f"e0 = {label}")
    for ax, xl in ((axes[0], "Lateral deflection at mid-height [mm]"), (axes[1], "Axial shortening [mm]")):
        ax.axhline(euler_load(L_DESIGN) / 1e3, color="tab:blue", ls="--", lw=1, label="Euler N_cr")
        ax.axhline(nb_rd / 1e3, color="tab:green", ls="-.", lw=1.3, label="Eurocode 3 N_b,Rd")
        ax.set_xlabel(xl)
        ax.set_ylabel("Axial load [kN]")
        ax.grid(alpha=0.3)
    axes[0].set_xlim(0, 150)
    axes[0].legend(fontsize=8)
    fig.suptitle("Imperfection sensitivity: the straighter the column, the closer it gets to Euler", fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{OUT}/imperfection_sensitivity.png", dpi=150)
    plt.close(fig)

    # 4. 3D collapse shape
    t_peak, t_end = design["time"][i_pk], design["time"][-1]
    fnodes, felems, frames = read_frd(f"collapse_L{int(L_DESIGN)}.frd", [t_peak, t_end])
    faces = boundary_faces(felems)
    fp, fe = frames[round(t_peak, 7)], frames[round(t_end, 7)]
    fig = plt.figure(figsize=(13, 7.5))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    norm_s = colors.Normalize(0, FY)
    polys = draw_column(ax1, fnodes, faces, fp["DISP"], fp["STRESS"], norm_s, cm.turbo, scale=5.0)
    set_view(ax1, polys, half_lateral=220, half_vertical=L_DESIGN / 2 + 50, box=(1, 1, 3.2))
    ax1.set_title(f"At collapse load ({nu_d / 1e3:.0f} kN)\nvon Mises stress, deformation x5\n(lateral axes stretched for visibility)", fontsize=10)
    fig.colorbar(cm.ScalarMappable(norm=norm_s, cmap=cm.turbo), ax=ax1, shrink=0.55, pad=0.0,
                 label="von Mises stress [MPa]")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    pe_max = max(fe["PE"].values())
    norm_p = colors.Normalize(0, pe_max)
    zone = (L_DESIGN / 2 - 450, L_DESIGN / 2 + 450)
    polys = draw_column(ax2, fnodes, faces, fe["DISP"], fe["PE"], norm_p, cm.inferno, scale=1.0, xlim=zone)
    set_view(ax2, polys, half_lateral=330, half_vertical=330, box=(1, 1, 1))
    ax2.set_title(f"End of analysis (load fell to {design['force'][-1] / 1e3:.0f} kN)\n"
                  "plastic hinge at mid-height, true scale", fontsize=10)
    fig.colorbar(cm.ScalarMappable(norm=norm_p, cmap=cm.inferno), ax=ax2, shrink=0.55, pad=0.0,
                 label="Equivalent plastic strain [-]")
    fig.suptitle("Collapse mechanism of the 3.0 m mezzanine column", fontsize=13)
    fig.tight_layout()
    fig.savefig(f"{OUT}/collapse_3d.png", dpi=140)
    plt.close(fig)

    with open(f"{OUT}/report.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(report))
    print(f"\nFigures, CSV and report saved in the '{OUT}' folder.")


if __name__ == "__main__":
    main()
