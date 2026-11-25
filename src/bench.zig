const std = @import("std");
const disk_tree = @import("disk_tree.zig");
const Circle = disk_tree.Circle;
const DiskTree = disk_tree.DiskTree;

const side = 1000.0;
const queries = 200_000;

fn seconds(from: std.Io.Timestamp, io: std.Io) f64 {
    return @as(f64, @floatFromInt(from.untilNow(io, .awake).nanoseconds)) / std.time.ns_per_s;
}

fn run(gpa: std.mem.Allocator, io: std.Io, n: usize, max_r: f64) !void {
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

    var start = std.Io.Timestamp.now(io, .awake);
    var tree = try DiskTree.init(gpa, circles);
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

    std.debug.print("{d:>8} {d:>6.1} {d:>9} {d:>9.3} {d:>8.1} {d:>7.2} {d:>9.0} {d:>10.0} {s}\n", .{
        n,
        max_r,
        tree.xs.len,
        build_s,
        @as(f64, @floatFromInt(tree.bytes())) / (1 << 20),
        @as(f64, @floatFromInt(reported)) / queries,
        query_ns,
        scan_ns,
        if (sampled == scanned) "ok" else "MISMATCH",
    });
}

pub fn main(init: std.process.Init) !void {
    std.debug.print("{s:>8} {s:>6} {s:>9} {s:>9} {s:>8} {s:>7} {s:>9} {s:>10} {s}\n", .{
        "disks", "max r", "slabs", "build s", "MiB", "hits", "query ns", "scan ns", "vs scan",
    });
    for ([_]usize{ 1_000, 10_000, 100_000, 1_000_000 }) |n| {
        // Radii chosen so a point is covered by about 1 and about 10 disks.
        for ([_]f64{ 1.0, 3.2 }) |depth| {
            const max_r = depth * side / @sqrt(@as(f64, @floatFromInt(n)));
            try run(init.gpa, init.io, n, max_r);
        }
    }
}
