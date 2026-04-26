# Basilar Membrane 3D Reconstruction — Pipeline Update

> **Master Thesis** · Yuan Fang · Supervisor: Lukas Driendl, M.Sc.
>
> Reconstructs patient-specific 3D surface meshes (STL) of the Basilar Membrane (BM)
> from µCT-derived Dragonfly point-cloud annotations, for integration into a COMSOL FEM
> model of the electrically stimulated human cochlea.
>
> This document supersedes [`README_BM.md`](README_BM.md). It documents the full pipeline
> evolution from the original `Basilar_test.py` approach to Lukas's two-script watertight
> pipeline, and the subsequent COMSOL end-cap triangulation fix.

---

## Overview

The pipeline has evolved in three stages:

1. **Original** (`Basilar_test.py`) — 9-step spine extraction; produces an open ribbon surface; confirmed as topologically insufficient for COMSOL solid meshing.
2. **Lukas's pipeline** (`interpolate_points_and_slices.py`) — 14-step closed-contour loft; produces a watertight solid STL; validated in Dragonfly (2026-04-07).
3. **Modified pipeline** (`interpolate_points_and_slices_modified.py`) — identical to Stage 2 except for end-cap triangulation, which is replaced with a 2D Delaunay + ray-casting filter to eliminate self-intersecting triangles that caused COMSOL meshing failures.

The **active script** for COMSOL submission is `interpolate_points_and_slices_modified.py`.

---

## Biomedical Background

### Cochlear Anatomy

The cochlea is a ~2.5-turn spiral bony structure in the inner ear. Two thin membranes divide its interior into three fluid-filled compartments:

```
┌──────────────────────────────────────┐
│         Scala Vestibuli              │  ← perilymph (low K⁺)
├──────────[Reissner's Membrane]───────┤
│         Scala Media                  │  ← endolymph (~80 mV EP, high K⁺)
├──────────[Basilar  Membrane]─────────┤
│         Scala Tympani                │  ← perilymph (low K⁺)
└──────────────────────────────────────┘
```

Both membranes act as **frequency-dependent electrical barriers**: primarily resistive at low frequencies, with capacitive displacement currents becoming relevant at higher CI stimulation frequencies.

### Why Closed Contours Are Required

The BM is modelled as a **thin elastic solid** (~80 µm thick) in the COMSOL FEM simulation. COMSOL requires the geometry to be a **closed solid body** (watertight STL) rather than an open surface, so that:

- the mesh generator can compute a volumetric tetrahedral mesh inside the BM domain,
- material properties (resistivity, density) can be assigned to the solid region,
- boundary conditions at BM–perilymph and BM–endolymph interfaces are well-defined.

An open ribbon surface (Stage 1) satisfies none of these requirements.

### Wimmer Algorithm — Role in This Project

> W. Wimmer, *Robust Cochlear Modiolar Axis Detection in CT*, Frontiers in Neuroscience, 2019.

The Wimmer algorithm fits a screw motion model `m = {r, c, γ}` to cochlear surface point clouds via a 7D Plücker coordinate eigenvalue problem:

```
7D feature vector:      φ(p, n) = [p × n,  n,  |p|² · n]
Spiral velocity field:  v(p)    = r × p  +  c  +  γ · p
```

`r, c, γ` were pre-fitted by K. Schroll (2025) from 2.74 M-point Scalae data.
In **Stage 1** they were used to anchor the modiolar end and sort slices base→apex.
In **Stages 2–3 (Lukas's pipeline)** they are **not used** — sorting is achieved via SVD on cluster centroids, making the pipeline fully self-contained.

---

## Algorithm Pipeline

### Stage 1 (Reference): `Basilar_test.py` — Open Ribbon (Superseded)

```
P1.1  Load & Preprocess
      ├── Read BM CSV  (units: metres → ×1000 → mm)
      └── Subtract Scalae_Slicer centroid  →  align to r, c, γ coordinate frame

P1.2  DBSCAN Clustering  (eps = 0.10 mm)
      └── Separate continuous point cloud into individual Dragonfly cross-sections

P1.3  PCA Spine Extraction  (25 pts / cluster)
      ├── PC1 = long axis of the membrane strip  (modiolar → lateral direction)
      ├── Sample 25 evenly-spaced spine points along PC1
      └── Anchor spine[0] = modiolar end  (min radial distance to modiolar axis)

P1.4  MST Sorting  (base → apex)
      ├── Radial distance R of each centroid from modiolar axis
      ├── MST on centroid Euclidean distances
      └── DFS from highest-R leaf  →  base→apex order

P1.5  Open-Curve Matching  (between consecutive slices)
      └── Forward vs. reversed diagonal cost  →  no longitudinal track crossings

P1.6  B-spline Track Smoothing
      └── Periodic cubic B-spline per track  →  300-pt resampling  (scipy splprep/splev)

P1.7  Restore Coordinates
      └── Add back Scalae centroid  →  original mm coordinate system

P1.8  Build Open Ribbon Mesh
      └── Quad-strip triangulation of adjacent spine tracks  (NOT cyclic, NOT watertight)

P1.9  Export STL
```

**Known limitation:** An open ribbon cannot be imported as a solid body in COMSOL, and shows only intersection-line wireframes (not filled cross-sections) in Dragonfly. Superseded by Stage 2.

---

### Stage 2 (Current Baseline): Two-Script Pipeline — Watertight Solid

```
BM raw point cloud CSV  (Dragonfly export, native metres)
    │
    ▼  group_points_for_slices.py
grouped CSV  (added slice_id column, ordered base → apex)
    │
    ▼  interpolate_points_and_slices_modified.py
BM_*_watertight_modified.stl  (watertight = True, Euler = 2)
```

#### Script 1: `group_points_for_slices.py`

```
Step 1  Load CSV  (x, y, z columns required)
        └── Print coordinate bounds

Step 2  Auto-estimate DBSCAN eps
        └── P95 of k=2 nearest-neighbour distances × 2.5  (sample_size = 5000 pts)

Step 3  DBSCAN Clustering
        └── min_samples=8, drop clusters < min_cluster_size=200

Step 4  Base→apex reordering
        ├── SVD on all cluster centroids  →  smallest singular vector ≈ modiolar axis
        ├── Project centroids onto perpendicular plane  →  radial distances
        ├── MST on centroid Euclidean distances
        └── DFS from highest-radius leaf  →  base→apex order

Step 5  Save grouped CSV  (<input>_grouped.csv)
        └── Adds integer slice_id column (0 = base, N−1 = apex)
```

#### Key Functions — `group_points_for_slices.py`

| Function | Location | Purpose |
|----------|----------|---------|
| `load_points()` | [`group_points_for_slices.py:69`](../Lukas/BM/group_points_for_slices.py#L69) | Loads raw CSV; validates required columns; prints bounds |
| `estimate_dbscan_eps()` | [`group_points_for_slices.py:96`](../Lukas/BM/group_points_for_slices.py#L96) | k-NN P95 × 2.5 automatic eps estimation |
| `build_mst_adjacency()` | [`group_points_for_slices.py:115`](../Lukas/BM/group_points_for_slices.py#L115) | Prim's MST on centroid point array |
| `renumber_clusters_base_to_apex()` | [`group_points_for_slices.py:142`](../Lukas/BM/group_points_for_slices.py#L142) | SVD projection + DFS traversal to assign ordered slice_id |
| `group_points_by_slice()` | [`group_points_for_slices.py:186`](../Lukas/BM/group_points_for_slices.py#L186) | Orchestrates DBSCAN + size filter + reordering |

---

#### Script 2 (Modified): `interpolate_points_and_slices_modified.py`

```
Step 1   Load & validate grouped CSV  (requires x, y, z, slice_id columns)

Step 2   SVD plane fit per slice
         └── Centroid + SVD → local 2D basis (u, v) + slice normal

Step 3   Estimate in-plane point spacing
         └── cKDTree P80 nearest-neighbour distance  →  alpha_radius = alpha_scale × spacing

Step 4   Alpha-shape boundary extraction in 2D
         ├── Delaunay + circumradius filter  →  boundary edges
         ├── 6 progressive alpha radii: 1.0 / 1.5 / 2.0 / 3.0 / 4.0 / 6.0 × alpha_radius
         └── Convex hull fallback if all radii fail

Step 5   Periodic cubic spline resampling
         └── Deduplicate → arc-length parameterisation → CubicSpline(periodic)
             → contour_points = 120 uniformly-spaced closed contour points

Step 6   Lift 2D contour back to 3D
         └── 3D point = centroid + u·basis[0] + v·basis[1]

Step 7   Align neighbouring contours
         └── For each contour: test 2 winding directions × contour_points cyclic shifts
             → keep shift+direction with minimum sum-of-point-distances to previous

Step 8   Per-vertex arc-length spline interpolation
         ├── Cumulative centroid arc-length as spline parameter
         ├── make_interp_spline (order 3) per vertex across all slices
         └── Evaluate at dense_slices = 250 uniformly-spaced parameter values

Step 9   Build lateral triangle mesh  (UNCHANGED from original)
         └── Quad-strip: contour[i][j] → contour[i][j+1] → contour[i+1][j] → ...
             → 2 triangles per quad, cyclic modulo contour_points

Step 10  Cap first and last contours  ← MODIFIED (see COMSOL End-Cap Fix below)
         └── triangulate_cap_2d(): 2D Delaunay + ray-casting interior filter
             replaces ear-clipping  →  correct for concave BM cross-sections

Step 11  Global face winding consistency
         └── Adjacency graph traversal → flip neighbours as needed
             → signed volume check → flip all if volume < 0

Step 12  trimesh cleanup
         └── Remove duplicate / degenerate faces; merge vertices

Step 13  Watertight validation
         └── Assert: is_watertight, is_winding_consistent, euler_number == 2, volume > 0

Step 14  Export STL  (+ optional surface-point CSV + PyVista preview)
```

#### Key Functions — `interpolate_points_and_slices_modified.py`

| Function | Location | Purpose |
|----------|----------|---------|
| `load_grouped_points()` | [`interpolate_points_and_slices_modified.py:120`](../Lukas/BM/interpolate_points_and_slices_modified.py#L120) | Loads grouped CSV; validates columns; sorts by slice_id |
| `fit_plane()` | [`interpolate_points_and_slices_modified.py:152`](../Lukas/BM/interpolate_points_and_slices_modified.py#L152) | SVD best-fit plane → local 2D basis + projected points |
| `alpha_shape_boundary()` | [`interpolate_points_and_slices_modified.py:264`](../Lukas/BM/interpolate_points_and_slices_modified.py#L264) | Delaunay + circumradius filter → closed boundary loop in 2D |
| `resample_closed_curve()` | [`interpolate_points_and_slices_modified.py:310`](../Lukas/BM/interpolate_points_and_slices_modified.py#L310) | Periodic cubic spline → N uniformly-spaced contour points |
| `align_contour_to_reference()` | [`interpolate_points_and_slices_modified.py:385`](../Lukas/BM/interpolate_points_and_slices_modified.py#L385) | Exhaustive cyclic shift × 2 windings → minimum-cost alignment |
| `interpolate_contours()` | [`interpolate_points_and_slices_modified.py:417`](../Lukas/BM/interpolate_points_and_slices_modified.py#L417) | Per-vertex arc-length spline → 250 dense interpolated contours |
| `triangulate_cap_2d()` | [`interpolate_points_and_slices_modified.py:558`](../Lukas/BM/interpolate_points_and_slices_modified.py#L558) | **NEW** — 2D Delaunay + ray-casting filter for end-cap triangulation |
| `contour_to_faces()` | [`interpolate_points_and_slices_modified.py:598`](../Lukas/BM/interpolate_points_and_slices_modified.py#L598) | Assembles all triangle faces (side loft + two end caps) |
| `build_watertight_mesh()` | [`interpolate_points_and_slices_modified.py:699`](../Lukas/BM/interpolate_points_and_slices_modified.py#L699) | Builds and cleans trimesh; returns validated watertight mesh |

---

## COMSOL End-Cap Fix

### Problem

Importing the original `_watertight.stl` into COMSOL Multiphysics 6.3.0 and running Free Tetrahedral meshing produced:

```
Error: Intersecting edge elements. (Edges 140, 156)
Coordinates: -0.00382639, -0.00688259, 0.00849347
Failed to generate mesh for face.
Failure due to incomplete boundary mesh.
```

The error was localised to the **Apex end cap**. The BM Apex cross-section is a narrow, elongated, **concave (non-convex) polygon**. Ear-clipping triangulation falls back to a fan triangulation at non-convex vertices, producing self-intersecting triangles that pass trimesh's topological watertight test but are rejected by COMSOL's geometric mesh checker.

### Fix Iteration History

| Attempt | Algorithm | trimesh result | COMSOL result |
|---------|-----------|----------------|---------------|
| Original | Ear-clipping, 120 pts | `watertight=True` ✓ | Intersecting edges 140/156 ✗ |
| Attempt 1 | Centroid fan, 40 pts (downsampled cap) | `Euler = −78` ✗ (open seam at cap/side junction) | — |
| Attempt 2 | Centroid fan, 120 pts | `watertight=True` ✓ | 7 intersecting edges ✗ |
| **Current** | **2D Delaunay + ray-casting, 120 pts** | `watertight=True` ✓ | **Pending validation** |

### Current Algorithm: `triangulate_cap_2d()`

```python
# Step 1 — project 3D cap contour to local SVD best-fit 2D plane
_, _, pts_2d = fit_plane(contour_3d)

# Step 2 — Delaunay triangulates all 120 boundary vertices in 2D
#           (covers convex hull; includes triangles inside and outside the polygon)
tri = Delaunay(pts_2d, qhull_options="QJ")

# Step 3 — ray-casting filter: keep only triangles whose centroid is inside polygon
valid = [s for s in tri.simplices
         if _point_in_polygon(pts_2d[s].mean(axis=0), pts_2d)]

# Step 4 — return triangle indices (same vertex ordering as contour_3d)
return np.array(valid, dtype=int)
```

**Why this is correct for concave polygons:**

| Property | Guarantee |
|----------|-----------|
| Delaunay triangles do not intersect each other in 2D | No self-intersecting triangles in the plane |
| The 2D→3D lift is a linear map (SVD basis) | Linear maps preserve the non-intersection property |
| All 120 original boundary vertices are used | No T-junctions at the cap/side-face seam → trimesh stays watertight |
| Ray-casting handles re-entrant (concave) corners correctly | Correct interior classification for the narrow Apex cross-section |

**Scope of modification — only end caps changed:**

| Region | Original | Modified |
|--------|----------|---------|
| Side-face loft (Steps 1–9, 11–14) | ear-clipping | **Unchanged** |
| Basal end cap (Step 10) | Ear-clipping, 120 pts | 2D Delaunay + ray-casting, 120 pts |
| Apex end cap (Step 10) | Ear-clipping, 120 pts | 2D Delaunay + ray-casting, 120 pts |

---

## Key Design Decision — Why Closed Contours, Not Open Spine?

| Issue | Open spine (`Basilar_test.py`) | Closed contour (Lukas's pipeline) |
|-------|-------------------------------|----------------------------------|
| Dragonfly cross-section view | Shows only intersection lines (wireframe) | Renders as filled solid region |
| COMSOL import type | Surface body (complex setup) | Solid body (direct import) |
| Physical BM thickness | Manual ±40 µm normal offset | Directly encoded in contour cross-section area |
| Normal vector dependency | Required (P1.8b was numerically fragile) | Not needed — winding direction from signed volume |
| `r, c, γ` parameters | Required for sorting + modiolar anchor | Not used — SVD self-contained |
| Slice alignment robustness | Forward/reverse only | Exhaustive cyclic shift × 2 windings |

> **Biomedical note:** The physical thickness of the BM (~80 µm) is encoded in the cross-sectional area of each extracted contour. Lukas's approach reads this geometry directly from the Dragonfly annotations and is physically more faithful than a manually prescribed thickness offset.

---

## Repository Structure

```
Master_Thesis/
│
├── GitHub/
│   ├── Basilar_test.py                  ← Stage 1: original 9-step spine pipeline (reference)
│   ├── Basilar_test_EN.md               ← Stage 1 session log (English)
│   ├── Basilar_test_CN.md               ← Stage 1 session log (Chinese)
│   ├── README_BM.md                     ← previous README (Stage 1)
│   └── README_modified.md               ← this file (Stages 1–3)
│
├── Lukas/
│   └── BM/
│       ├── group_points_for_slices.py               ← Script 1: DBSCAN + base→apex sort
│       ├── interpolate_points_and_slices.py          ← Script 2: original (Lukas's version)
│       ├── interpolate_points_and_slices_modified.py ← Script 2: COMSOL cap fix  ← ACTIVE
│       ├── BM_Lukas_EN.md                           ← Lukas pipeline session log (English)
│       ├── BM_Lukas_CN.md                           ← Lukas pipeline session log (Chinese)
│       ├── BM_single_slice_watertight_pipeline.md   ← Lukas's official pipeline documentation
│       └── BM_Dragonfly_final/
│           ├── BM_Slices_24Apr_test7.csv                            ← raw BM point cloud (current)
│           ├── BM_Slices_24Apr_test7_grouped.csv                    ← Script 1 output
│           ├── BM_Slices_24Apr_test7_grouped_watertight.stl         ← Script 2 original output
│           └── BM_Slices_24Apr_test7_grouped_watertight_modified.stl ← Script 2 modified  ← ACTIVE
│
├── BM_Apex_modified/
│   └── BM_Apex_Sparse_verylong.csv   ← Stage 1 BM input (re-annotated apex, 2026-03-25)
│
├── Scalae_Slicer.csv                 ← 2.74 M pts, mm — source of r, c, γ (Stage 1 only)
└── Reissner_Membrane.stl             ← RM input for Phase 3 (not yet processed)
```

---

## Requirements

**Python 3.9.6** (virtual environment recommended)

```bash
pip install numpy scipy pandas scikit-learn pyvista trimesh open3d \
            shapely alphashape matplotlib vtk
```

All packages are pre-installed in the project virtual environment at `.venv/`.

PyVista is configured for **windowed (non-notebook) rendering** — running inside a Jupyter notebook requires setting `pv.global_theme.notebook = True`.

---

## How to Run

```bash
# Activate the virtual environment (from workspace root)
cd /Users/huoyu/Documents/Visual_Studio_Code
source .venv/bin/activate
```

### Step 1 — Group point cloud into ordered slices

```bash
python Master_Thesis/Lukas/BM/group_points_for_slices.py \
  --csv Master_Thesis/Lukas/BM/BM_Dragonfly_final/BM_Slices_24Apr_test7.csv
```

Output: `BM_Slices_24Apr_test7_grouped.csv`

### Step 2 — Watertight STL reconstruction (COMSOL-ready)

```bash
# Default parameters
python Master_Thesis/Lukas/BM/interpolate_points_and_slices_modified.py

# With PyVista preview
python Master_Thesis/Lukas/BM/interpolate_points_and_slices_modified.py --show
```

Output: `BM_Slices_24Apr_test7_grouped_watertight_modified.stl`

The `--show` flag opens a PyVista window displaying:

| Layer | Colour | Content |
|-------|--------|---------|
| Point cloud | Tomato (12% opacity) | Original input points |
| Black lines | Black | Extracted slice contours (one per annotated slice) |
| Dense contours | Royal blue | 250 interpolated contours |
| Surface | Steel blue (70% opacity) | Final reconstructed BM solid |

### Tunable Parameters

| Parameter | CLI flag | Default | When to adjust |
|-----------|----------|---------|----------------|
| `contour_points` | `--contour-points` | `120` | Increase for smoother caps; decrease to reduce COMSOL mesh complexity |
| `dense_slices` | `--dense-slices` | `250` | Increase for smoother longitudinal surface |
| `alpha_scale` | `--alpha-scale` | `8.0` | Decrease if boundary over-smooths concave features; increase if boundary fragments |
| `eps` (Script 1) | `--eps` | auto | Override if DBSCAN produces wrong cluster count |
| `min_cluster_size` (Script 1) | `--min-cluster-size` | `200` | Decrease if valid apex slices are being dropped |

---

## Data Reference

### Input / Output Files

| File | Path | Units | Notes |
|------|------|-------|-------|
| BM raw point cloud | `Lukas/BM/BM_Dragonfly_final/BM_Slices_24Apr_test7.csv` | **metres** | Current Dragonfly export |
| BM grouped (Script 1 output) | `Lukas/BM/BM_Dragonfly_final/BM_Slices_24Apr_test7_grouped.csv` | metres | Contains `slice_id` column |
| BM STL — original | `Lukas/BM/BM_Dragonfly_final/BM_Slices_24Apr_test7_grouped_watertight.stl` | metres | Lukas's baseline output |
| BM STL — COMSOL fix | `Lukas/BM/BM_Dragonfly_final/BM_Slices_24Apr_test7_grouped_watertight_modified.stl` | metres | **Active output** |
| Scalae point cloud | `Master_Thesis/Scalae_Slicer.csv` | mm | 2.74 M pts; `r, c, γ` source (Stage 1 only) |
| Reissner Membrane | `Master_Thesis/Reissner_Membrane.stl` | metres | Phase 3 input; not yet processed |

> **Unit note:** All STL files use native Dragonfly **metre-scale coordinates**. When importing into Dragonfly, select `Meters` in the Mesh Import dialog (the default `Millimeters` produces a 1000× scale error). When importing into COMSOL, the unit must also be set to `m`.

### Spiral Parameters (Stage 1 only — not used in Stages 2–3)

```python
r     = np.array([ 0.15667362, -0.69239267, -0.62562308])  # modiolar axis direction
c     = np.array([-0.01322494,  0.07648561,  0.30903837])  # translation component
gamma = 0.05578412602375563                                  # scaling factor
# p₀  ≈ [0.592, -3.481, -3.116] mm  (modiolar axis anchor, centred coordinates)
```

Fitted by K. Schroll (2025) from `Scalae_Slicer.csv`. Valid for the whole cochlea; do **not** refit from BM/RM data.

---

## Project Status

| Phase | Description | Status |
|-------|-------------|--------|
| Phase 0 | `r, c, γ` pre-computed from `Scalae_Slicer.csv` (K. Schroll) | Complete |
| Phase 1 | `Basilar_test.py` — open ribbon BM reconstruction | Complete (reference; superseded) |
| Phase 2a | Lukas pipeline — watertight STL, local trimesh validation | Complete |
| Phase 2b | Dragonfly import + size / position validation (unit fix) | Complete (2026-04-07) |
| Phase 2c | COMSOL end-cap fix — `interpolate_points_and_slices_modified.py` | Implemented; **COMSOL validation pending** |
| Phase 2d | Dragonfly µCT slice-by-slice boundary accuracy overlay | Pending |
| Phase 3 | Reissner Membrane reconstruction | Not started |
| Final | COMSOL FEM integration + simulation | Not started |

### Phase 2c Validation Checklist

- [x] `triangulate_cap_2d()` implemented — 2D Delaunay + ray-casting interior filter
- [x] trimesh reports `watertight=True`, `winding_consistent=True`, `Euler=2`, `volume>0`
- [x] NumPy 2.0 deprecation warning in `triangle_circumradius` fixed
- [ ] COMSOL Free Tetrahedral meshing — no `Intersecting edge elements` on Apex cap
- [ ] COMSOL Free Tetrahedral meshing — no `Intersecting edge elements` on Basal cap
- [ ] COMSOL mesh quality: minimum element quality > 0.1

### Phase 2d Validation Checklist

- [x] STL imports into Dragonfly with correct size and position (Meters unit)
- [ ] Slice-by-slice µCT cross-section overlay — BM boundary traces anatomy at Base
- [ ] Slice-by-slice µCT cross-section overlay — BM boundary traces anatomy at Mid-turn
- [ ] Slice-by-slice µCT cross-section overlay — BM boundary traces anatomy at Apex

---

## Next Steps

### Phase 2c — COMSOL Validation

```
□  Import BM_Slices_24Apr_test7_grouped_watertight_modified.stl into COMSOL 6.3.0
□  Geometry → Form Union
□  Mesh → Free Tetrahedral → confirm no "Intersecting edge elements"
□  If error persists: inspect which edges → report coordinates back for diagnosis
□  If successful: run BM material assignment + boundary condition setup
```

### Phase 2d — Dragonfly Boundary Accuracy

```
□  Open µCT dataset in Dragonfly
□  Import modified STL with unit = Meters
□  Overlay with 2D cross-sections at Base, Mid-turn, Apex
□  Assess BM boundary alignment: does contour trace the µCT anatomy?
□  Optional: measure Hausdorff distance between STL boundary and Dragonfly annotation
```

### Phase 3 — Reissner Membrane

```
□  Confirm RM Dragonfly annotation format: x, y, z per cross-section (same as BM)?
□  Export RM point cloud as CSV from Dragonfly
□  Run group_points_for_slices.py  (adjust eps and min_cluster_size for RM geometry)
□  Run interpolate_points_and_slices_modified.py  (adjust alpha_scale if needed)
□  Validate RM STL in Dragonfly and COMSOL
```

### Open Questions for Discussion

Two integration paths for COMSOL are under evaluation — decision to be made with Siwei and Lukas:

| Path | Description | Trade-offs |
|------|-------------|-----------|
| **(a)** Insert into existing model | Add BM/RM STLs as solid bodies to the current COMSOL model | Lower effort; coordinate alignment must be verified |
| **(b)** New µCT pipeline | Start fresh from a new µCT scan with all 3 scalae segmented | Cleaner geometry; higher upfront effort |

Membrane electrical properties (resistivity values for COMSOL boundary conditions):
→ Source: https://www.mdpi.com/2076-3417/14/22/10408

---

## References

1. **W. Wimmer** (2019). *Robust Cochlear Modiolar Axis Detection in CT*. Frontiers in Neuroscience.
   Screw motion parameterisation used in Stage 1. Detection accuracy: distance error 0.13 mm, angular error 2.4° (n=23 specimens).

2. **K. Schroll** (2025). *3D Parametrization of Cochlear Scalae* (Master Thesis, TUM).
   Source of `r, c, γ` parameters used in Stage 1. Implements Wimmer algorithm on Scala Tympani, Media, and Vestibuli.

3. **L. Driendl** (2026). *interpolate_points_and_slices.py* — internal GitLab repository.
   Two-script watertight BM reconstruction pipeline adopted as Stage 2 baseline in this project.

4. **Membrane electrical properties** (2024). *Applied Sciences* 14(22), 10408.
   https://www.mdpi.com/2076-3417/14/22/10408
   Resistivity values for BM and RM needed for COMSOL material assignment.
