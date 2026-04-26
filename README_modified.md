# Basilar Membrane 3D Reconstruction — Pipeline Update

> **Master Thesis** · Yuan Fang · Supervisor: Lukas Driendl, M.Sc.
>
> This document is an updated version of `README_BM.md`, written for GitLab communication.
> It documents the full pipeline evolution: from the original `Basilar_test.py` approach
> to Lukas's two-script watertight pipeline, and the subsequent COMSOL end-cap fix.
>
> _Last updated: 2026-04-26_

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

Both membranes act as **frequency-dependent electrical barriers** that significantly modulate how stimulation current spreads from cochlear implant (CI) electrodes through the cochlea.

**Thesis objective:** Reconstruct patient-specific 3D surface meshes (STL) of the Basilar Membrane (BM) and Reissner's Membrane (RM) from µCT-derived point cloud annotations, for integration into a COMSOL FEM model of the electrically stimulated human cochlea.

---

## Pipeline Evolution

### Phase 1 — Original approach: `Basilar_test.py` (9-step spine pipeline)

The initial pipeline modelled the BM as a **1D spine** extracted per cross-section, then connected spines into an open ribbon surface.

**Key steps:**

```
P1.1  Load BM CSV (metres → ×1000 → mm) + subtract Scalae_Slicer centroid
P1.2  DBSCAN clustering — separate cross-sectional slices (eps = 0.10 mm)
P1.3  PCA spine extraction — 25 spine points per cluster; modiolar-end anchor
P1.4  MST + DFS sorting — base → apex order (radial distance from p₀)
P1.5  Open-curve matching — forward/reverse diagonal cost comparison
P1.6  B-spline smoothing — 300-pt cubic spline per longitudinal track
P1.7  Restore coordinates — add back Scalae centroid
P1.8  Build open ribbon mesh — quad-strip triangulation (NOT closed)
P1.9  Export STL
```

**Why `r, c, γ` were used:** Pre-fitted spiral parameters from Katharina Schroll's 2025 thesis (fitted from 2.74 M-point Scalae data) provided the modiolar axis anchor `p₀` for radial sorting and the coordinate frame for centring. The BM is too thin (~80 µm) for reliable normal estimation, so Wimmer's algorithm was not re-run on BM data.

**Known limitations discovered:**

| Limitation | Impact |
|-----------|--------|
| Open ribbon surface (non-watertight) | Dragonfly cross-section view shows only intersection lines, not a filled region ("parallel wireframe" artifact) |
| ±40 µm manual normal offset for thickness | Numerically fragile on dense/open-boundary meshes |
| COMSOL import | Requires open-surface boundary conditions; complicates meshing |

---

### Phase 2 — Lukas's two-script watertight pipeline (current baseline)

After reviewing the limitations above, Lukas's pipeline was adopted as the new baseline.

**Core architectural difference:** instead of extracting a 1D spine per slice, Lukas's approach extracts a **closed 2D cross-sectional contour** per slice. Lofting closed contours + adding end caps naturally produces a watertight solid.

#### Two-script architecture

```
BM point cloud CSV  (raw Dragonfly export, metres)
    │
    ▼  group_points_for_slices.py
grouped CSV  (with slice_id, ordered base → apex)
    │
    ▼  interpolate_points_and_slices.py
BM_*_watertight.stl  (watertight = True, Euler = 2)
```

#### Script 1: `group_points_for_slices.py`

| Step | Operation |
|------|-----------|
| eps estimation | P95 of k=2 nearest-neighbour distances × 2.5 (fully automatic) |
| DBSCAN clustering | `min_samples=8`, `min_cluster_size=200` |
| Base→apex sorting | SVD on cluster centroids → radial distances → MST + DFS from highest-radius leaf |
| Output | CSV with `slice_id` column added |

Note: no `r, c, γ` parameters are used — the sorting is fully self-contained via SVD.

#### Script 2: `interpolate_points_and_slices.py` (14-step pipeline)

| Step | Operation |
|------|-----------|
| 1 | Load and validate grouped CSV |
| 2 | SVD plane fit → local 2D basis per slice |
| 3 | Estimate characteristic in-plane point spacing (cKDTree P80) |
| 4 | Alpha-shape boundary extraction in 2D (6 progressive radii; convex-hull fallback) |
| 5 | Periodic cubic spline resampling → 120 uniformly-spaced points per closed contour |
| 6 | Lift 2D contour back to 3D |
| 7 | Align neighbouring contours: exhaustive cyclic shift × 2 winding directions |
| 8 | Per-vertex arc-length cubic spline → 250 dense interpolated contours |
| 9 | Build side-face triangle mesh (quad-strip loft, cyclic around perimeter) |
| 10 | Cap first and last contours (ear-clipping) → closed solid |
| 11 | Global face winding consistency (adjacency graph + signed volume check) |
| 12 | trimesh cleanup (remove duplicates, degenerate faces, merge vertices) |
| 13 | Watertight validation (`is_watertight` / `is_winding_consistent` / Euler=2 / volume>0) |
| 14 | Export STL (+ optional surface-point CSV + PyVista preview) |

**Verified result (Lukas's reference run):**

```
watertight            = True
winding_consistent    = True
Euler number          = 2
volume                > 0
```

**Comparison with original approach:**

| Aspect | `Basilar_test.py` | Lukas's pipeline |
|--------|------------------|-----------------|
| Geometric primitive per slice | 1D spine (25 pts, open) | 2D closed contour (120 pts) |
| Output topology | Open ribbon | Watertight solid |
| Physical BM thickness | ±40 µm manual offset | Directly from contour shape |
| `r, c, γ` parameters | Required (sorting + anchor) | Not used |
| Slice alignment | Forward/reverse cost only | Exhaustive cyclic + winding |
| COMSOL import | Surface body (complex) | Solid body (direct) |

**Dragonfly validation (2026-04-07):** Import with unit set to `Meters` confirmed size and spatial position correctly match the µCT anatomy. Full slice-by-slice µCT overlay still pending.

---

### Phase 3 — COMSOL end-cap fix: `interpolate_points_and_slices_modified.py`

**Problem:** Importing `BM_Slices_24Apr_test7_grouped_watertight.stl` into COMSOL Multiphysics 6.3.0 produced:

```
Error: Intersecting edge elements. (Edges 140, 156)
Coordinates: -0.00382639, -0.00688259, 0.00849347
Failed to generate mesh for face.
Failure due to incomplete boundary mesh.
```

**Root cause:** The BM cross-section is a **concave (non-convex) polygon**. The original ear-clipping triangulation in Step 10 falls back to a fan triangulation for non-convex vertices, producing self-intersecting triangles. These are geometrically valid in trimesh's watertight check (which tests topology, not planar self-intersection) but are rejected by COMSOL's secondary meshing.

**Fix development history:**

| Attempt | Algorithm | trimesh | COMSOL |
|---------|-----------|---------|--------|
| Original | Ear-clipping, 120 pts | ✓ watertight | ✗ Edges 140/156 intersecting |
| Attempt 1 | Centroid fan, 40 pts (downsampled cap) | ✗ Euler = −78 (T-junction on side/cap seam) | — |
| Attempt 2 | Centroid fan, 120 pts | ✓ watertight | ✗ 7 intersecting edges |
| **Current** | **2D Delaunay + ray-casting, 120 pts** | ✓ (expected) | ⏳ validation pending |

**Current algorithm (`triangulate_cap_2d` in `interpolate_points_and_slices_modified.py`):**

```
Step 1  fit_plane(): project the 3D cap contour to its local SVD best-fit 2D plane
Step 2  scipy.spatial.Delaunay on all 120 boundary vertices in 2D
        (covers the convex hull — includes triangles inside and outside the polygon)
Step 3  _point_in_polygon() ray-casting filter:
        - for each Delaunay triangle, test whether its centroid is inside the polygon
        - ray-casting is correct for both convex and concave (non-convex) polygons
        - discard exterior triangles
Step 4  return filtered triangle indices (same vertex ordering as contour_3d)
```

**Why this is correct for concave polygons:**
- Delaunay triangles are non-intersecting in 2D by construction
- The 2D→3D lift is a linear map (via SVD basis), which preserves non-intersection
- Using all 120 original boundary vertices avoids T-junctions at the cap/side-face seam
- Ray-casting correctly handles concave re-entrant corners

**Scope of modification — only end-cap logic changed:**

| Region | Change |
|--------|--------|
| Side-face loft (Steps 1–9) | Unchanged |
| Basal end cap | Ear-clipping 120 pts → 2D Delaunay + ray-casting 120 pts |
| Apex end cap | Ear-clipping 120 pts → 2D Delaunay + ray-casting 120 pts |

**Also fixed:** NumPy 2.0 deprecation warning in `triangle_circumradius` — replaced `np.cross(b-a, c-a)` with explicit 2D scalar cross product.

---

## Current Project Status

| Phase | Description | Status |
|-------|-------------|--------|
| Phase 0 | Scalae spiral parameters (`r, c, γ`) pre-computed by K. Schroll | Complete |
| Phase 1 | `Basilar_test.py` — open ribbon BM reconstruction | Complete (reference; superseded) |
| Phase 2a | Lukas pipeline — watertight BM STL, local validation | Complete (`watertight=True`, Euler=2) |
| Phase 2b | Dragonfly import + size/position validation | Complete (unit fix: import as Meters) |
| Phase 2c | COMSOL end-cap fix (`interpolate_points_and_slices_modified.py`) | Implemented; **COMSOL validation pending** |
| Phase 2d | Dragonfly µCT slice-by-slice overlay (boundary accuracy) | Pending |
| Phase 3 | Reissner Membrane reconstruction | Not started |
| Final | COMSOL FEM integration + simulation | Not started |

---

## Repository Structure

```
Master_Thesis/
│
├── GitHub/
│   ├── Basilar_test.py            ← original 9-step spine pipeline (reference)
│   ├── Basilar_test_EN.md         ← original pipeline session log (English)
│   ├── Basilar_test_CN.md         ← original pipeline session log (Chinese)
│   ├── README_BM.md               ← previous README
│   └── README_modified.md         ← this file
│
├── Lukas/
│   └── BM/
│       ├── group_points_for_slices.py                  ← Script 1: DBSCAN slice grouping
│       ├── interpolate_points_and_slices.py             ← Script 2: original (Lukas's version)
│       ├── interpolate_points_and_slices_modified.py    ← Script 2: COMSOL cap fix (active)
│       ├── BM_Lukas_EN.md                              ← session log (English)
│       ├── BM_Lukas_CN.md                              ← session log (Chinese)
│       ├── BM_single_slice_watertight_pipeline.md      ← Lukas's pipeline documentation
│       └── BM_Dragonfly_final/
│           ├── BM_Slices_24Apr_test7.csv               ← raw BM point cloud (current)
│           ├── BM_Slices_24Apr_test7_grouped.csv       ← grouped (Script 1 output)
│           ├── BM_Slices_24Apr_test7_grouped_watertight.stl       ← original Script 2 output
│           └── BM_Slices_24Apr_test7_grouped_watertight_modified.stl  ← cap-fix output (active)
│
├── BM_Apex_modified/
│   └── BM_Apex_Sparse_verylong.csv    ← older BM input (Phase 1, re-annotated apex)
│
├── Scalae_Slicer.csv                  ← 2.74 M pts, mm — source of r, c, γ
└── Reissner_Membrane.stl              ← RM input (Phase 3, not yet processed)
```

---

## How to Run

**Prerequisites:**

```bash
cd /Users/huoyu/Documents/Visual_Studio_Code
source .venv/bin/activate
```

### Script 1 — Group point cloud into slices

```bash
python Master_Thesis/Lukas/BM/group_points_for_slices.py \
  --csv Master_Thesis/Lukas/BM/BM_Dragonfly_final/BM_Slices_24Apr_test7.csv
```

Output: `BM_Slices_24Apr_test7_grouped.csv`

### Script 2 — Watertight STL reconstruction (original)

```bash
python Master_Thesis/Lukas/BM/interpolate_points_and_slices.py \
  --csv Master_Thesis/Lukas/BM/BM_Dragonfly_final/BM_Slices_24Apr_test7_grouped.csv \
  --show
```

### Script 2 modified — COMSOL-compatible cap fix (active)

```bash
# Default (120-point caps, 2D Delaunay + ray-casting)
python Master_Thesis/Lukas/BM/interpolate_points_and_slices_modified.py --show

# Custom parameters
python Master_Thesis/Lukas/BM/interpolate_points_and_slices_modified.py \
  --contour-points 120 \
  --dense-slices 250 \
  --alpha-scale 8.0 \
  --show
```

Output: `BM_Slices_24Apr_test7_grouped_watertight_modified.stl`

---

## Key Parameters

| Parameter | Default | Script | Effect |
|-----------|---------|--------|--------|
| `eps` | auto | `group_points_for_slices.py` | DBSCAN neighbourhood radius (P95 × 2.5) |
| `min_cluster_size` | 200 | `group_points_for_slices.py` | Minimum points per detected slice |
| `contour_points` | 120 | `interpolate_*.py` | Circumferential resolution per contour |
| `dense_slices` | 250 | `interpolate_*.py` | Longitudinal resolution (interpolated slices) |
| `alpha_scale` | 8.0 | `interpolate_*.py` | Alpha-shape boundary tightness |

**Unit note:** All scripts preserve native Dragonfly **metre-scale coordinates**. When importing the output STL into Dragonfly, select `Meters` in the Mesh Import dialog (not the default `Millimeters`).

---

## Validation Checklist

### trimesh (programmatic, already passing)

```python
import trimesh
mesh = trimesh.load("BM_Slices_24Apr_test7_grouped_watertight_modified.stl")
assert mesh.is_watertight           # True
assert mesh.is_winding_consistent   # True
assert mesh.euler_number == 2       # True
assert mesh.volume > 0              # True
```

### COMSOL (pending for modified STL)

1. Import `BM_Slices_24Apr_test7_grouped_watertight_modified.stl` as a **Mesh** object
2. Under **Geometry**, run **Form Union**
3. Apply **Free Tetrahedral** mesh — target: no `Intersecting edge elements` error on Apex or Basal caps
4. Check mesh quality: minimum element quality > 0.1

### Dragonfly (slice-by-slice µCT overlay, pending)

1. Import STL with unit = `Meters`
2. Overlay with source µCT cross-sections
3. Verify BM boundary traces the µCT anatomy at Basal, Mid-turn, and Apex regions

---

## Open Questions for Discussion with Lukas

1. **COMSOL validation:** Is the modified STL (`_modified.stl`) now accepted by COMSOL Free Tetrahedral meshing without intersecting-edge errors?

2. **Dragonfly µCT overlay:** How should we quantify the alignment between the reconstructed BM and the µCT anatomy — are there target metrics (e.g., Hausdorff distance, per-slice boundary error)?

3. **Units — Option B:** Should the scripts apply an in-code ×1000 conversion (metres → mm) so that the default Dragonfly import unit (Millimeters) works without manual adjustment on every import?

4. **Phase 3 — Reissner Membrane:** Does the RM Dragonfly annotation have the same format as BM (`x, y, z` point cloud per cross-section)? Can the same two-script pipeline be applied directly?

5. **COMSOL integration path:** Should both BM and RM be imported as thin **solid bodies** (watertight STLs) into the existing COMSOL model, or is a new full µCT pipeline (path b) preferred?

---

## References

1. **W. Wimmer** (2019). *Robust Cochlear Modiolar Axis Detection in CT*. Frontiers in Neuroscience.
   Mathematical foundation for the screw motion parameterisation (`r, c, γ`). Distance error 0.13 mm, angular error 2.4° (n=23).

2. **Membrane electrical properties** (2024). *Applied Sciences* 14(22), 10408.
   Resistivity values for BM and RM needed for COMSOL boundary conditions.

3. **K. Schroll** (2025). *3D Parametrization of Cochlear Scalae* (Master Thesis, TUM).
   Source of the `r, c, γ` parameters used in `Basilar_test.py`.

4. **L. Driendl** (2026). *interpolate_points_and_slices.py* — internal GitLab repository.
   Two-script watertight BM reconstruction pipeline; adopted as current baseline.
