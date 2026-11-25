# Ideas

Notes from 2026-11-25 on where this could go. Nothing here is built unless it
says so.

## The one real edge

`fold` returns the combined value of every shape over a point in O(log n), no
matter how many shapes that is. With shapes listed back to front and
premultiplied colours, that is a finished alpha-blended pixel. Anything based on
bounding boxes has to visit every layer.

Measured (disks, Apple M1 Pro): at about 100 layers per point, `fold` took
0.7 µs for 1,000 disks and 2.2 µs for 10,000, against 3 µs and 7 µs to list the
layers.

The price is memory: one slab per crossing, 100 to 200 bytes each, and the
number of slabs grows with shapes times depth squared.

## The engine idea: push overdraw to the extreme

A renderer where transparent layers are free per pixel, so the art style can be
built on them: glass panes over everything, every blade of grass translucent,
crystal palaces, cel shading from many stacked flat bands.

Why it may be new: games solve the *ordering* of transparent layers (depth
peeling, per-pixel fragment lists, weighted blending) but still shade every
layer. Vector renderers composite layer by layer inside tiles. The closest
relative is PDF transparency flattening in print, which precomputes overlap
regions offline. Resolving overlap ahead of time and sampling it live is the
part that seems unexplored. This is from memory, not a literature search; do
one before committing months.

It loses to a GPU on ordinary scenes. It only wins at depths nobody designs for
today, so the engine has to be built around that style.

## Making the memory affordable

Fold results can themselves be combined, which allows two things:

- **Instancing.** One tree per asset (a grass clump, a crystal cluster, a
  window), placed many times with a transform. A plain BVH over the instances
  finds the few covering a pixel; each answers in O(log n); combine front to
  back. Memory is per unique asset. Moving, rotating and scaling an instance is
  just a transform on the query point.
- **Depth groups.** Crossings between shapes in different groups need not be
  stored. G depth bands cut crossings by about G and cost G lookups per pixel.
  Glass "over the top" is one extra band.

So: a conventional structure outside, this one inside wherever overlap is deep
and local.

## Motion

- Rigid motion of an instance is free (see above).
- Rebuilding per frame caps at about a thousand shapes: builds measured about
  7 ms for 1,000 sparse Bézier shapes and 0.8 s for 100,000.
- Precomputing all of time ("4D") does not work: a version is needed for every
  change in how curves cross, and three-curve coincidences alone can number
  around n³.
- A kinetic version is plausible and is research. Between topology changes only
  coordinates move, so arcs could be evaluated at (x, t) and each frame would
  only refresh slab boundaries. Two crossings swapping order in x is an easy
  O(log n) edit. A crossing pair appearing or vanishing at a tangency inserts or
  deletes slabs mid-history, which the current append-only persistence cannot
  do; it would need full persistence or a retroactive design.

## Gradients

`fold` is not limited to flat colours. With constant alpha per shape and colour
linear in (x, y), "over" keeps the result linear, so the coefficients combine
associatively. Radial gradients and alpha that varies across a shape do not fit.
So many flat bands do not beat one linear gradient; they beat things a gradient
cannot express.

## As a file format

Not a replacement for SVG: it is derived, not editable, and far larger. It could
be a baked runtime format shipped next to the source art, like a distance-field
atlas for fonts but for arbitrary overlapping translucent shapes:
resolution-independent, exact edges, any point sampled in O(log n) without
rasterising. Good fits: zoomable maps, vector art on 3D surfaces.

## What is unproven, riskiest first

1. **Sampling on a GPU.** Everything runs on a CPU at a few microseconds per
   sample, roughly 1 to 5 seconds for a 1080p frame on one thread. It needs a
   flattened layout readable from a fragment shader, and a dozen dependent reads
   per pixel may be too slow.
2. **Instancing and depth groups**, as above.
3. **Motion inside an asset**: bending, morphing.
4. **Antialiasing.** Queries are point samples. Supersampling works and
   multiplies the cost.

## First experiment

A shader that samples one baked tree full-screen, timed at depth 10, 100 and
1000. If it holds 60 fps the rest is engineering. If not, that is known in days.
