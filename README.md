# disk-tree

A disk stabbing structure in Zig: given n disks in the plane, report every disk
that contains a query point in O(log n + k) time, where k is the number of disks
reported.

```zig
const disk_tree = @import("disk_tree");

var tree = try disk_tree.DiskTree.init(gpa, &.{
    .{ .cx = 0, .cy = 0, .r = 2 },
    .{ .cx = 1, .cy = 0.5, .r = 1 },
});
defer tree.deinit(gpa);

const ids = try tree.query(gpa, 0.5, 0.2); // { 0, 1 }
defer gpa.free(ids);
```

`queryAppend` does the same into a list you already have.

## Folding without visiting

`Folded(T, combine)` attaches a value to every disk and returns the values of
all disks over a point already combined, in O(log n) no matter how many disks
that is. `combine` must be associative but need not be commutative: it sees the
disks in index order.

Counting the disks over a point:

```zig
fn add(a: u32, b: u32) u32 {
    return a + b;
}

var tree = try disk_tree.Folded(u32, add).init(gpa, circles, ones);
const depth = tree.fold(x, y) orelse 0;
```

With disks listed back to front, premultiplied colours and the "over" operator,
`fold` is the finished pixel. `zig build render -- out.ppm` draws 600
translucent disks that way, one lookup per pixel:

![600 translucent disks](disks.png)

This is what a bounding volume hierarchy cannot do: it has to visit every disk
over the point. The price is the memory below, so it pays off when the overlap
is deep, the values are constant per disk, and the queries are not coherent
enough to sweep directly.

## Bézier shapes

The repository also holds the stabbing code from `circle_renderer`, which
handles convex shapes bounded by cubic Béziers rather than disks:

| file | what it has |
|------|-------------|
| `src/math/cubic.zig` | cubic, quadratic and linear root finding |
| `src/geometry/bezier.zig` | cubic Bézier evaluation, extrema, splitting into x-monotonic pieces |
| `src/geometry/shapes.zig` | `ConvexBezierShape`: y-interval at an x, point containment |
| `src/stabbing/stabbable.zig` | `Stabbable`, the interface a shape offers to the tree |
| `src/stabbing/slab.zig` | one slab and its segment tree over y |
| `src/stabbing/tree.zig` | `StabbingTree`: slab decomposition over any `Stabbable` shapes, and the packing of it for a GPU |

It is a separate structure from `DiskTree`: it keeps a list of active shapes per
slab instead of persistent versions, and works in `f32`. Everything is exported
from `src/root.zig`.

Two things do not work as they stand. `StabbingTree.query` returns nothing,
because `build` records each slab's active shapes without building its segment
tree, so the `stabbing tree basic` test fails. `PackedTreeData` does not compile
on Zig 0.16, because it uses `std.Thread.Pool`.

## How the disk tree works

Sweep a vertical line across the plane. It cuts each circle it meets in two
points, one on the lower arc and one on the upper arc, and the arcs it cuts have
an order by y. Between two neighbouring arcs lies one face of the arrangement,
and every point of a face is covered by the same disks.

That order only changes at an event: a circle starts, a circle ends, or two
circles cross. So the x axis splits into slabs with a fixed order each, and
neighbouring slabs differ by a small edit:

- a start inserts two arcs,
- an end removes two arcs,
- a crossing swaps two neighbouring arcs, which changes the face between them
  by one disk.

The arcs are kept in a treap ordered by y, and each arc stores the set of disks
covering the face above it, with the combined value of each subtree kept in its
nodes. Both are persistent, so every slab keeps its own
version at O(log n) extra cost instead of a full copy. Treap nodes hold one
timestamped override before they are copied (Driscoll et al.), and the disk sets
are path-copied treaps.

A query binary searches the slabs by x, walks that slab's treap by y, and
enumerates the set stored on the arc just below the point.

With c circle-circle crossings, building takes O((n + c) log n) time and space
after finding the crossing pairs. c is at most n(n - 1), and that is reached
when every pair of circles crosses, so this structure suits inputs where disks
overlap a bounded number of others.

### Floating point

Crossing coordinates are computed in floating point, so the sweep does not trust
them to say what the order is. Events only mark arcs as possibly out of place.
The order itself comes from comparing those arcs with their neighbours at the
middle of the next slab, and every swap that comparison causes updates the face
sets. Several circles through one point, tangencies and duplicate circles need
no special cases.

A query point within rounding error of a circle's boundary may be reported
either way for that circle.

## Build

Needs Zig 0.16.

```sh
zig build test
zig build bench
```

## Benchmark

Uniformly random disks in a 1000 x 1000 square, 200,000 random queries, on an
Apple M1 Pro. "hits" is the average number of disks over a query point, "query"
lists them, "fold" counts them with `Folded`, and "scan" is a linear pass over
all disks.

| disks     | hits  | slabs     | build  | memory  | query   | fold    | scan    |
|-----------|-------|-----------|--------|---------|---------|---------|---------|
| 1,000     | 1.2   | 5,722     | 6 ms   | 0.6 MiB | 0.39 µs | 0.24 µs | 0.45 µs |
| 1,000     | 11.3  | 37,490    | 16 ms  | 5.0 MiB | 0.61 µs | 0.39 µs | 0.46 µs |
| 1,000     | 95.9  | 260,674   | 0.13 s | 50 MiB  | 3.3 µs  | 0.75 µs | 0.44 µs |
| 10,000    | 1.2   | 56,820    | 25 ms  | 5.7 MiB | 0.43 µs | 0.41 µs | 4.4 µs  |
| 10,000    | 11.6  | 390,748   | 0.18 s | 53 MiB  | 1.7 µs  | 1.0 µs  | 4.4 µs  |
| 10,000    | 108.0 | 3,345,438 | 2.1 s  | 672 MiB | 7.5 µs  | 2.4 µs  | 4.6 µs  |
| 100,000   | 1.2   | 578,216   | 0.36 s | 59 MiB  | 1.6 µs  | 1.5 µs  | 45 µs   |
| 100,000   | 11.8  | 4,041,312 | 2.4 s  | 549 MiB | 3.6 µs  | 2.6 µs  | 44 µs   |
| 1,000,000 | 1.2   | 5,800,312 | 4.5 s  | 591 MiB | 3.8 µs  | 3.6 µs  | 463 µs  |

The number of slabs grows with disks times hits squared, and memory is 100 to
200 bytes per slab, so deep overlap is expensive to store. Query time on large
inputs is dominated by cache misses on the walk down the treap; a thousand
disks are still quicker to scan.

## References

- Sarnak and Tarjan, "Planar point location using persistent search trees", 1986
- Driscoll, Sarnak, Sleator and Tarjan, "Making data structures persistent", 1989
