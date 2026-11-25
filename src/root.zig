//! Point stabbing structures: which shapes contain a point?

const std = @import("std");

/// Disks, with a persistent sweep: O(log n + k) queries and O(log n) folds.
pub const disk_tree = @import("disk_tree.zig");
pub const Circle = disk_tree.Circle;
pub const DiskTree = disk_tree.DiskTree;
pub const Folded = disk_tree.Folded;

/// Shapes bounded by cubic Béziers, with the same persistent sweep.
pub const path_tree = @import("path_tree.zig");
pub const Path = path_tree.Path;
pub const PathTree = path_tree.PathTree;
pub const FoldedPaths = path_tree.FoldedPaths;

/// The sweep itself, for other kinds of shape.
pub const sweep = @import("sweep.zig");

/// Convex shapes bounded by cubic Béziers, with a segment tree per slab.
pub const cubic = @import("math/cubic.zig");
pub const bezier = @import("geometry/bezier.zig");
pub const shapes = @import("geometry/shapes.zig");
pub const stabbable = @import("stabbing/stabbable.zig");
pub const slab = @import("stabbing/slab.zig");
pub const tree = @import("stabbing/tree.zig");

test {
    inline for (.{
        @This(),                  cubic,                   bezier,
        bezier.CubicBezier,       bezier.MonotonicSegment, shapes,
        shapes.ConvexBezierShape, stabbable,               slab,
        slab.Slab,                tree,                    tree.StabbingTree,
        tree.PackedTreeData,      sweep,                   path_tree,
    }) |namespace| std.testing.refAllDecls(namespace);
}

fn randomShapes(gpa: std.mem.Allocator, random: std.Random, convex: []shapes.ConvexBezierShape, handles: []stabbable.Stabbable) !void {
    for (convex, handles) |*shape, *handle| {
        shape.* = try shapes.createRandomConvexShape(
            gpa,
            random,
            random.float(f32) * 100,
            random.float(f32) * 100,
            1 + random.float(f32) * 15,
            3 + random.uintLessThan(u32, 6),
        );
        handle.* = shape.asStabbable();
    }
}

test "StabbingTree agrees with a scan" {
    const gpa = std.testing.allocator;
    var prng: std.Random.DefaultPrng = .init(1);
    const random = prng.random();

    var convex: [120]shapes.ConvexBezierShape = undefined;
    var handles: [convex.len]stabbable.Stabbable = undefined;
    try randomShapes(gpa, random, &convex, &handles);
    defer for (&convex) |*shape| shape.deinit();

    var stabbing = tree.StabbingTree.init(gpa);
    defer stabbing.deinit();
    try stabbing.build(&handles);

    var got: std.ArrayList(u32) = .empty;
    defer got.deinit(gpa);
    var want: std.ArrayList(u32) = .empty;
    defer want.deinit(gpa);
    for (0..5000) |_| {
        const x = random.float(f32) * 110 - 5;
        const y = random.float(f32) * 110 - 5;
        want.clearRetainingCapacity();
        for (&convex, 0..) |*shape, i| {
            if (shape.containsPoint(x, y)) try want.append(gpa, @intCast(i));
        }
        got.clearRetainingCapacity();
        try stabbing.query(x, y, &got);
        std.mem.sortUnstable(u32, got.items, {}, std.sort.asc(u32));
        try std.testing.expectEqualSlices(u32, want.items, got.items);
    }
}

test "PackedTreeData packs with and without threads" {
    const gpa = std.testing.allocator;
    var prng: std.Random.DefaultPrng = .init(2);
    const random = prng.random();

    var convex: [120]shapes.ConvexBezierShape = undefined;
    var handles: [convex.len]stabbable.Stabbable = undefined;
    try randomShapes(gpa, random, &convex, &handles);
    defer for (&convex) |*shape| shape.deinit();

    var stabbing = tree.StabbingTree.init(gpa);
    defer stabbing.deinit();
    try stabbing.build(&handles);

    var packed_tree = tree.PackedTreeData.init(gpa);
    defer packed_tree.deinit();
    try packed_tree.pack(&stabbing, &handles, -20, 120);
    try std.testing.expect(packed_tree.segment_data.items.len > 0);
}
