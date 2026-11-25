const std = @import("std");
const disk_tree = @import("disk_tree.zig");
const path_tree = @import("path_tree.zig");
const shapes = @import("geometry/shapes.zig");
const stabbable = @import("stabbing/stabbable.zig");
const stabbing = @import("stabbing/tree.zig");
const Circle = disk_tree.Circle;

fn add(a: u32, b: u32) u32 {
    return a + b;
}

/// Folding a 1 per disk counts the disks over a point.
const Tree = disk_tree.Folded(u32, add);
const PathTree = path_tree.FoldedPaths(u32, add);

const side = 1000.0;
const queries = 200_000;

fn seconds(from: std.Io.Timestamp, io: std.Io) f64 {
    return @as(f64, @floatFromInt(from.untilNow(io, .awake).nanoseconds)) / std.time.ns_per_s;
}

fn runDisks(gpa: std.mem.Allocator, io: std.Io, n: usize, max_r: f64) !void {
    var prng: std.Random.DefaultPrng = .init(n);
    const random = prng.random();

    const circles = try gpa.alloc(Circle, n);
    defer gpa.free(circles);
    for (circles) |*c| c.* = .{
        .cx = random.float(f64) * side,
        .cy = random.float(f64) * side,
        .r = (0.1 + 0.9 * random.float(f64)) * max_r,
    };
    const points = try gpa.alloc([2]f64, queries);
    defer gpa.free(points);
    for (points) |*p| p.* = .{ random.float(f64) * side, random.float(f64) * side };

    const ones = try gpa.alloc(u32, n);
    defer gpa.free(ones);
    @memset(ones, 1);

    var start = std.Io.Timestamp.now(io, .awake);
    var tree = try Tree.init(gpa, circles, ones);
    defer tree.deinit(gpa);
    const build_s = seconds(start, io);

    var out: std.ArrayList(u32) = .empty;
    defer out.deinit(gpa);
    var reported: usize = 0;
    start = std.Io.Timestamp.now(io, .awake);
    for (points) |p| {
        out.clearRetainingCapacity();
        try tree.queryAppend(gpa, &out, p[0], p[1]);
        reported += out.items.len;
    }
    const query_ns = seconds(start, io) * 1e9 / queries;

    var counted: usize = 0;
    start = std.Io.Timestamp.now(io, .awake);
    for (points) |p| counted += tree.fold(p[0], p[1]) orelse 0;
    const fold_ns = seconds(start, io) * 1e9 / queries;

    // A linear scan over a sample, as the baseline and as a check.
    const sample = points[0..@min(queries, 200_000_000 / n)];
    var scanned: usize = 0;
    start = std.Io.Timestamp.now(io, .awake);
    for (sample) |p| {
        for (circles) |c| scanned += @intFromBool(c.contains(p[0], p[1]));
    }
    const scan_ns = seconds(start, io) * 1e9 / @as(f64, @floatFromInt(sample.len));

    var sampled: usize = 0;
    for (sample) |p| {
        out.clearRetainingCapacity();
        try tree.queryAppend(gpa, &out, p[0], p[1]);
        sampled += out.items.len;
    }

    std.debug.print("{d:>8} {d:>6.1} {d:>9} {d:>9.3} {d:>8.1} {d:>7.2} {d:>9.0} {d:>8.0} {d:>10.0} {s}\n", .{
        n,
        max_r,
        tree.xs.len,
        build_s,
        @as(f64, @floatFromInt(tree.bytes())) / (1 << 20),
        @as(f64, @floatFromInt(reported)) / queries,
        query_ns,
        fold_ns,
        scan_ns,
        if (sampled == scanned and counted == reported) "ok" else "MISMATCH",
    });
}

/// Convex Bézier shapes: the path tree against the per-slab segment trees of
/// `StabbingTree` and a scan.
fn runPaths(gpa: std.mem.Allocator, io: std.Io, n: usize, max_size: f32) !void {
    var prng: std.Random.DefaultPrng = .init(n);
    const random = prng.random();

    const convex = try gpa.alloc(shapes.ConvexBezierShape, n);
    defer gpa.free(convex);
    const paths = try gpa.alloc(path_tree.Path, n);
    defer gpa.free(paths);
    const handles = try gpa.alloc(stabbable.Stabbable, n);
    defer gpa.free(handles);
    for (convex, paths, handles) |*shape, *path, *handle| {
        shape.* = try shapes.createRandomConvexShape(
            gpa,
            random,
            random.float(f32) * side,
            random.float(f32) * side,
            (0.1 + 0.9 * random.float(f32)) * max_size,
            3 + random.uintLessThan(u32, 6),
        );
        path.* = .{ .curves = shape.curves.items };
        handle.* = shape.asStabbable();
    }
    defer for (convex) |*shape| shape.deinit();

    const ones = try gpa.alloc(u32, n);
    defer gpa.free(ones);
    @memset(ones, 1);
    const points = try gpa.alloc([2]f32, queries);
    defer gpa.free(points);
    for (points) |*p| p.* = .{ random.float(f32) * side, random.float(f32) * side };

    var start = std.Io.Timestamp.now(io, .awake);
    var tree = try PathTree.init(gpa, paths, ones);
    defer tree.deinit(gpa);
    const build_s = seconds(start, io);

    var out: std.ArrayList(u32) = .empty;
    defer out.deinit(gpa);
    var reported: usize = 0;
    start = std.Io.Timestamp.now(io, .awake);
    for (points) |p| {
        out.clearRetainingCapacity();
        try tree.queryAppend(gpa, &out, p[0], p[1]);
        reported += out.items.len;
    }
    const query_ns = seconds(start, io) * 1e9 / queries;

    var counted: usize = 0;
    start = std.Io.Timestamp.now(io, .awake);
    for (points) |p| counted += tree.fold(p[0], p[1]) orelse 0;
    const fold_ns = seconds(start, io) * 1e9 / queries;

    const sample = points[0..@min(queries, 20_000_000 / n)];
    var scanned: usize = 0;
    start = std.Io.Timestamp.now(io, .awake);
    for (sample) |p| {
        for (0..n) |shape| scanned += @intFromBool(tree.geometry.contains(shape, p[0], p[1]));
    }
    const scan_ns = seconds(start, io) * 1e9 / @as(f64, @floatFromInt(sample.len));
    var sampled: usize = 0;
    for (sample) |p| sampled += tree.fold(p[0], p[1]) orelse 0;

    // The slab tree stores every slab's active shapes, which is quadratic.
    var slab_build_s: f64 = std.math.nan(f64);
    var slab_query_ns: f64 = std.math.nan(f64);
    if (n <= 10_000) {
        start = std.Io.Timestamp.now(io, .awake);
        var slabs = stabbing.StabbingTree.init(gpa);
        defer slabs.deinit();
        try slabs.build(handles);
        slab_build_s = seconds(start, io);

        start = std.Io.Timestamp.now(io, .awake);
        for (points) |p| {
            out.clearRetainingCapacity();
            try slabs.query(p[0], p[1], &out);
        }
        slab_query_ns = seconds(start, io) * 1e9 / queries;
    }

    std.debug.print("{d:>8} {d:>6.1} {d:>9} {d:>9.3} {d:>8.1} {d:>7.2} {d:>9.0} {d:>8.0} {d:>10.0} {d:>10.3} {d:>9.0} {s}\n", .{
        n,
        max_size,
        tree.xs.len,
        build_s,
        @as(f64, @floatFromInt(tree.bytes())) / (1 << 20),
        @as(f64, @floatFromInt(reported)) / queries,
        query_ns,
        fold_ns,
        scan_ns,
        slab_build_s,
        slab_query_ns,
        if (sampled == scanned and counted == reported) "ok" else "MISMATCH",
    });
}

pub fn main(init: std.process.Init) !void {
    std.debug.print("{s:>8} {s:>6} {s:>9} {s:>9} {s:>8} {s:>7} {s:>9} {s:>8} {s:>10} {s}\n", .{
        "disks", "max r", "slabs", "build s", "MiB", "hits", "query ns", "fold ns", "scan ns", "vs scan",
    });
    for ([_]usize{ 1_000, 10_000, 100_000, 1_000_000 }) |n| {
        // Radii chosen so a point is covered by about 1, 10 and 100 disks.
        for ([_]f64{ 1.0, 3.2, 10.0 }) |depth| {
            // Slabs, and so memory, grow with n * depth^2; stay under a GiB.
            if (@as(f64, @floatFromInt(n)) * depth * depth > 5e6) continue;
            const max_r = depth * side / @sqrt(@as(f64, @floatFromInt(n)));
            try runDisks(init.gpa, init.io, n, max_r);
        }
    }

    std.debug.print("\n{s:>8} {s:>6} {s:>9} {s:>9} {s:>8} {s:>7} {s:>9} {s:>8} {s:>10} {s:>10} {s:>9} {s}\n", .{
        "shapes", "size", "slabs", "build s", "MiB", "hits", "query ns", "fold ns", "scan ns", "slab bld s", "slab q ns", "vs scan",
    });
    for ([_]usize{ 1_000, 10_000, 100_000 }) |n| {
        for ([_]f32{ 1.0, 3.2 }) |depth| {
            const max_size = depth * side / @sqrt(@as(f32, @floatFromInt(n)));
            try runPaths(init.gpa, init.io, n, max_size);
        }
    }
}
