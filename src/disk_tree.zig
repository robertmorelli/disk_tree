//! Disk stabbing: preprocess a set of disks so that "which disks contain this
//! point?" is answered in O(log n + k) for k reported disks, and "what do the
//! values of those disks combine to?" in O(log n). See `sweep.zig` for how.

const std = @import("std");
const Allocator = std.mem.Allocator;
const sweep = @import("sweep.zig");
const Event = sweep.Event;

pub const Circle = struct {
    cx: f64,
    cy: f64,
    r: f64,

    pub fn contains(c: Circle, x: f64, y: f64) bool {
        const dx = x - c.cx;
        const dy = y - c.cy;
        return dx * dx + dy * dy <= c.r * c.r;
    }

    /// Circles that fail this are ignored: they are never reported.
    fn valid(c: Circle) bool {
        return std.math.isFinite(c.cx) and std.math.isFinite(c.cy) and
            std.math.isFinite(c.r * c.r) and c.cx - c.r < c.cx + c.r;
    }
};

/// A disk tree without per-disk values.
///
/// `init(gpa, circles)` copies the circles; ids reported by queries are indices
/// into them. Building takes O((n + c) log n) time and space for c
/// circle-circle crossings, after a pair search that is O(n log n) plus the
/// number of circles overlapping each other along the longer axis.
pub const DiskTree = Folded(void, struct {
    fn combine(_: void, _: void) void {}
}.combine);

/// A disk tree that also carries a value of type `T` per disk and can combine
/// the values of all disks containing a point without visiting them.
/// `init(gpa, circles, values)` takes one value per circle.
///
/// `combine(a, b)` must be associative. It is always called with the disks
/// behind `a` at lower indices than those behind `b`, so it need not be
/// commutative: with disks listed back to front, alpha compositing works.
pub fn Folded(comptime T: type, comptime combine: fn (T, T) T) type {
    return sweep.Tree(Circles, T, combine);
}

/// Each disk has two arcs: `2 * disk` is the lower half of its circle and
/// `2 * disk + 1` the upper half.
const Circles = struct {
    circles: []const Circle,

    pub const Input = []const Circle;

    pub fn init(gpa: Allocator, circles: []const Circle, events: *std.ArrayList(Event)) Allocator.Error!Circles {
        const owned = try gpa.dupe(Circle, circles);
        errdefer gpa.free(owned);
        try collectEvents(gpa, owned, events);
        return .{ .circles = owned };
    }

    pub fn deinit(self: *Circles, gpa: Allocator) void {
        gpa.free(self.circles);
    }

    pub fn shapeCount(self: *const Circles) usize {
        return self.circles.len;
    }

    pub fn arcCount(self: *const Circles) usize {
        return 2 * self.circles.len;
    }

    pub fn shapeOf(_: *const Circles, arc: u32) u32 {
        return arc >> 1;
    }

    pub fn arcY(self: *const Circles, arc: u32, x: f64) f64 {
        const c = self.circles[arc >> 1];
        const dx = x - c.cx;
        const h = @sqrt(@max(0.0, c.r * c.r - dx * dx));
        return if (arc & 1 == 0) c.cy - h else c.cy + h;
    }

    /// Crossing coordinates lose about half their digits near a tangency.
    pub fn slack(self: *const Circles, arc: u32) f64 {
        const c = self.circles[arc >> 1];
        return 1e-7 * @sqrt(c.r * @max(c.r, @abs(c.cx), @abs(c.cy)));
    }

    pub fn bytes(self: *const Circles) usize {
        return std.mem.sliceAsBytes(self.circles).len;
    }
};

/// Emits a start and an end event per circle, and an event per point where two
/// circles cross. Crossings are only hints telling the sweep which arcs to
/// re-examine, so near-tangent pairs are included rather than risk missing one.
fn collectEvents(gpa: Allocator, circles: []const Circle, events: *std.ArrayList(Event)) Allocator.Error!void {
    var boxes: std.ArrayList(sweep.Box) = .empty;
    defer boxes.deinit(gpa);
    try boxes.ensureTotalCapacity(gpa, circles.len);
    for (circles, 0..) |c, i| {
        if (!c.valid()) continue;
        const id: u32 = @intCast(i);
        try events.append(gpa, .{ .x = c.cx - c.r, .a = 2 * id, .b = @intCast(2 * id + 1), .kind = .start });
        try events.append(gpa, .{ .x = c.cx + c.r, .a = 2 * id, .b = @intCast(2 * id + 1), .kind = .end });
        boxes.appendAssumeCapacity(.{ .x = .{ c.cx - c.r, c.cx + c.r }, .y = .{ c.cy - c.r, c.cy + c.r }, .id = id });
    }
    try sweep.forEachOverlap(gpa, boxes.items, Crossings{ .gpa = gpa, .circles = circles, .events = events });
}

const Crossings = struct {
    gpa: Allocator,
    circles: []const Circle,
    events: *std.ArrayList(Event),

    /// Emits an event per point where circles `ia` and `ib` cross.
    pub fn pair(self: Crossings, ia: u32, ib: u32) Allocator.Error!void {
        const a = self.circles[ia];
        const b = self.circles[ib];
        const dx = b.cx - a.cx;
        const dy = b.cy - a.cy;
        const d2 = dx * dx + dy * dy;
        if (d2 == 0) return;

        const slack = 1.0 + 1e-9;
        const sum = a.r + b.r;
        const diff = a.r - b.r;
        if (d2 > sum * sum * slack or d2 * slack < diff * diff) return;

        // The crossing points are the foot of the radical line on the centre
        // line, offset by h along the perpendicular.
        const t = (a.r * a.r - b.r * b.r + d2) / (2 * d2);
        const h = @sqrt(@max(0.0, a.r * a.r / d2 - t * t));
        const mx = a.cx + t * dx;
        const my = a.cy + t * dy;

        inline for (.{ -h, h }) |s| {
            const px = mx + s * dy;
            const py = my - s * dx;
            try self.events.append(self.gpa, .{
                .x = px,
                .a = 2 * ia + @intFromBool(py > a.cy),
                .b = @intCast(2 * ib + @intFromBool(py > b.cy)),
                .kind = .cross,
            });
        }
    }
};

const testing = std.testing;

/// Checks `tree` against a linear scan at `queries` random points in `box`
/// (min x, min y, max x, max y), skipping points too close to a boundary to
/// have a well-defined answer.
fn expectMatchesScan(circles: []const Circle, box: [4]f64, queries: usize, random: std.Random) !void {
    const gpa = testing.allocator;
    var tree = try DiskTree.init(gpa, circles);
    defer tree.deinit(gpa);

    var got: std.ArrayList(u32) = .empty;
    defer got.deinit(gpa);
    var want: std.ArrayList(u32) = .empty;
    defer want.deinit(gpa);

    var done: usize = 0;
    points: while (done < queries) {
        const x = box[0] + random.float(f64) * (box[2] - box[0]);
        const y = box[1] + random.float(f64) * (box[3] - box[1]);

        want.clearRetainingCapacity();
        for (circles, 0..) |c, i| {
            if (!c.valid()) continue;
            const d = std.math.hypot(x - c.cx, y - c.cy);
            if (@abs(d - c.r) < 1e-9 * @max(1.0, c.r)) continue :points;
            if (d < c.r) try want.append(gpa, @intCast(i));
        }
        done += 1;

        got.clearRetainingCapacity();
        try tree.queryAppend(gpa, &got, x, y);
        try testing.expectEqualSlices(u32, want.items, got.items);
    }
}

fn expectQuery(tree: *const DiskTree, x: f64, y: f64, want: []const u32) !void {
    const got = try tree.query(testing.allocator, x, y);
    defer testing.allocator.free(got);
    try testing.expectEqualSlices(u32, want, got);
}

test "no circles" {
    var tree = try DiskTree.init(testing.allocator, &.{});
    defer tree.deinit(testing.allocator);
    try expectQuery(&tree, 0, 0, &.{});
}

test "overlapping, nested and disjoint" {
    const circles = [_]Circle{
        .{ .cx = 0, .cy = 0, .r = 2 },
        .{ .cx = 1, .cy = 0.5, .r = 1 },
        .{ .cx = -1, .cy = 0.5, .r = 1 },
        .{ .cx = 0, .cy = 0, .r = 5 },
        .{ .cx = 20, .cy = 0, .r = 1 },
    };
    var tree = try DiskTree.init(testing.allocator, &circles);
    defer tree.deinit(testing.allocator);

    try expectQuery(&tree, 0.5, 0.2, &.{ 0, 1, 3 });
    try expectQuery(&tree, -0.5, 0.3, &.{ 0, 2, 3 });
    try expectQuery(&tree, 0, 0.6, &.{ 0, 3 });
    try expectQuery(&tree, 0, -1.5, &.{ 0, 3 });
    try expectQuery(&tree, 3, 3, &.{3});
    try expectQuery(&tree, 20.5, 0.5, &.{4});
    try expectQuery(&tree, 10, 0, &.{});
    try expectQuery(&tree, -1000, 0, &.{});
    try expectQuery(&tree, 1000, 0, &.{});
    try expectQuery(&tree, 0, 1000, &.{});
    try expectQuery(&tree, 0, -1000, &.{});
    try expectQuery(&tree, std.math.nan(f64), 0, &.{});
    try expectQuery(&tree, 0, std.math.nan(f64), &.{});
}

test "invalid circles are ignored" {
    const circles = [_]Circle{
        .{ .cx = 0, .cy = 0, .r = 0 },
        .{ .cx = 0, .cy = 0, .r = -1 },
        .{ .cx = 0, .cy = 0, .r = std.math.nan(f64) },
        .{ .cx = 0, .cy = 0, .r = std.math.inf(f64) },
        .{ .cx = std.math.inf(f64), .cy = 0, .r = 1 },
        .{ .cx = 1e300, .cy = 0, .r = 1 },
        .{ .cx = 0, .cy = 0, .r = 1 },
    };
    var tree = try DiskTree.init(testing.allocator, &circles);
    defer tree.deinit(testing.allocator);
    try expectQuery(&tree, 0.1, 0.1, &.{6});
}

test "random circles" {
    var prng: std.Random.DefaultPrng = .init(1);
    const random = prng.random();
    var circles: [300]Circle = undefined;

    for (0..40) |round| {
        // Sparse to heavily overlapping.
        const max_r = 0.5 + @as(f64, @floatFromInt(round));
        const n = 1 + random.uintLessThan(usize, circles.len);
        for (circles[0..n]) |*c| c.* = .{
            .cx = random.float(f64) * 100,
            .cy = random.float(f64) * 100,
            .r = 0.01 + random.float(f64) * max_r,
        };
        try expectMatchesScan(circles[0..n], .{ -5, -5, 105, 105 }, 2000, random);
    }
}

test "tall and narrow" {
    var prng: std.Random.DefaultPrng = .init(2);
    const random = prng.random();
    var circles: [300]Circle = undefined;
    for (&circles) |*c| c.* = .{
        .cx = random.float(f64),
        .cy = random.float(f64) * 1000,
        .r = 1 + random.float(f64) * 20,
    };
    try expectMatchesScan(&circles, .{ -25, -25, 25, 1025 }, 5000, random);
}

test "mixed scales" {
    var prng: std.Random.DefaultPrng = .init(3);
    const random = prng.random();
    var circles: [300]Circle = undefined;
    for (&circles) |*c| c.* = .{
        .cx = 1e6 + random.float(f64) * 10,
        .cy = -1e6 + random.float(f64) * 10,
        .r = std.math.pow(f64, 10, random.float(f64) * 6 - 4),
    };
    try expectMatchesScan(&circles, .{ 1e6 - 1, -1e6 - 1, 1e6 + 11, -1e6 + 11 }, 5000, random);
}

test "integer grid" {
    // Many circles through common points, tangencies, and shared event xs.
    var prng: std.Random.DefaultPrng = .init(4);
    const random = prng.random();
    var circles: std.ArrayList(Circle) = .empty;
    defer circles.deinit(testing.allocator);

    for ([_]f64{ 0.5, 1, @sqrt(2.0), 2, 2.5 }) |r| {
        circles.clearRetainingCapacity();
        for (0..10) |i| for (0..10) |j| {
            try circles.append(testing.allocator, .{ .cx = @floatFromInt(i), .cy = @floatFromInt(j), .r = r });
        };
        try expectMatchesScan(circles.items, .{ -3, -3, 12, 12 }, 5000, random);
    }

    circles.clearRetainingCapacity();
    for (0..200) |_| try circles.append(testing.allocator, .{
        .cx = @floatFromInt(random.uintLessThan(u32, 12)),
        .cy = @floatFromInt(random.uintLessThan(u32, 12)),
        .r = @floatFromInt(1 + random.uintLessThan(u32, 4)),
    });
    try expectMatchesScan(circles.items, .{ -5, -5, 16, 16 }, 5000, random);
}

test "concentric and duplicate" {
    var prng: std.Random.DefaultPrng = .init(5);
    const random = prng.random();
    var circles: [200]Circle = undefined;
    for (&circles, 0..) |*c, i| c.* = .{ .cx = 3, .cy = -2, .r = @floatFromInt(1 + i / 2) };
    try expectMatchesScan(&circles, .{ -100, -105, 106, 101 }, 5000, random);
}

test "rings through one point" {
    // Every circle passes through the origin.
    var prng: std.Random.DefaultPrng = .init(6);
    const random = prng.random();
    var circles: [64]Circle = undefined;
    for (&circles, 0..) |*c, i| {
        const angle = std.math.tau * @as(f64, @floatFromInt(i)) / circles.len;
        c.* = .{ .cx = 5 * @cos(angle), .cy = 5 * @sin(angle), .r = 5 };
    }
    try expectMatchesScan(&circles, .{ -11, -11, 11, 11 }, 5000, random);
    try expectMatchesScan(&circles, .{ -0.01, -0.01, 0.01, 0.01 }, 5000, random);
}

test "fold matches a scan" {
    // Composition of affine maps mod a prime: associative, not commutative.
    const Affine = struct {
        a: u64,
        b: u64,

        const p = 1_000_003;

        fn then(f: @This(), g: @This()) @This() {
            return .{ .a = f.a * g.a % p, .b = (f.b * g.a + g.b) % p };
        }
    };
    const Tree = Folded(Affine, Affine.then);
    const gpa = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(7);
    const random = prng.random();

    var circles: [300]Circle = undefined;
    var values: [300]Affine = undefined;
    for ([_]f64{ 2, 10, 40 }) |max_r| {
        for (&circles, &values) |*c, *v| {
            c.* = .{ .cx = random.float(f64) * 100, .cy = random.float(f64) * 100, .r = 0.01 + random.float(f64) * max_r };
            v.* = .{ .a = 1 + random.uintLessThan(u64, Affine.p - 1), .b = random.uintLessThan(u64, Affine.p) };
        }
        var tree = try Tree.init(gpa, &circles, &values);
        defer tree.deinit(gpa);

        for (0..5000) |_| {
            const x = random.float(f64) * 100;
            const y = random.float(f64) * 100;
            const ids = try tree.query(gpa, x, y);
            defer gpa.free(ids);

            var want: ?Affine = null;
            for (ids) |id| want = if (want) |w| w.then(values[id]) else values[id];
            try testing.expectEqual(want, tree.fold(x, y));
        }
    }
}
