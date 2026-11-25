//! Point stabbing structures: which shapes contain a point?

const std = @import("std");

/// Disks, with a persistent sweep: O(log n + k) queries and O(log n) folds.
pub const disk_tree = @import("disk_tree.zig");
pub const Circle = disk_tree.Circle;
pub const DiskTree = disk_tree.DiskTree;
pub const Folded = disk_tree.Folded;

/// Convex shapes bounded by cubic Béziers, with a segment tree per slab.
pub const cubic = @import("math/cubic.zig");
pub const bezier = @import("geometry/bezier.zig");
pub const shapes = @import("geometry/shapes.zig");
pub const stabbable = @import("stabbing/stabbable.zig");
pub const slab = @import("stabbing/slab.zig");
pub const tree = @import("stabbing/tree.zig");

test {
    // tree.PackedTreeData is left out: it uses std.Thread.Pool, which Zig
    // 0.16 removed.
    inline for (.{
        @This(),                  cubic,                   bezier,
        bezier.CubicBezier,       bezier.MonotonicSegment, shapes,
        shapes.ConvexBezierShape, stabbable,               slab,
        slab.Slab,                tree,                    tree.StabbingTree,
    }) |namespace| std.testing.refAllDecls(namespace);
}
