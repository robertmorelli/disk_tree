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

## How it works

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
covering the face above it. Both are persistent, so every slab keeps its own
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
Apple M1 Pro. "hits" is the average number of disks per query, and "scan" is a
linear pass over all disks.

| disks     | hits | slabs      | build   | memory   | query   | scan   |
|-----------|------|------------|---------|----------|---------|--------|
| 1,000     | 1.2  | 5,722      | 2 ms    | 0.5 MiB  | 0.27 µs | 0.45 µs |
| 1,000     | 11.3 | 37,490     | 15 ms   | 4.4 MiB  | 0.53 µs | 0.45 µs |
| 10,000    | 1.2  | 56,820     | 25 ms   | 5.4 MiB  | 0.43 µs | 4.4 µs |
| 10,000    | 11.6 | 390,748    | 0.18 s  | 47 MiB   | 1.8 µs  | 4.4 µs |
| 100,000   | 1.2  | 578,216    | 0.33 s  | 56 MiB   | 1.7 µs  | 45 µs  |
| 100,000   | 11.8 | 4,041,312  | 2.3 s   | 483 MiB  | 3.7 µs  | 44 µs  |
| 1,000,000 | 1.2  | 5,800,312  | 4.4 s   | 561 MiB  | 3.8 µs  | 481 µs |
| 1,000,000 | 11.9 | 40,793,948 | 33 s    | 4.9 GiB  | 7.6 µs  | 470 µs |

Memory is 100 to 125 bytes per slab, and the last row needs several GiB more
while building. Query time on large inputs is dominated by
cache misses on the walk down the treap.

## References

- Sarnak and Tarjan, "Planar point location using persistent search trees", 1986
- Driscoll, Sarnak, Sleator and Tarjan, "Making data structures persistent", 1989
