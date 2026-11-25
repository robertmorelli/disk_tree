//! Stabbing for shapes bounded by cubic Béziers: "which shapes contain this
//! point?" in O(log n + k), and "what do their values combine to?" in O(log n).
//!
//! A shape is one or more closed contours, filled by the even-odd rule, so it
//! may be concave, have holes, or cross itself. Straight edges are cubics with
//! their control points on the line. See `sweep.zig` for the structure.

const std = @import("std");
const Allocator = std.mem.Allocator;
const sweep = @import("sweep.zig");
const Event = sweep.Event;
const bezier = @import("geometry/bezier.zig");

pub const CubicBezier = bezier.CubicBezier;
pub const Vec2 = bezier.Vec2;

/// The outline of one shape. Consecutive curves that do not join start a new
/// contour, and a contour that does not end where it began is closed with a
/// straight edge.
pub const Path = struct {
    curves: []const CubicBezier,
    /// Added to every control point.
    offset: Vec2 = .{ .x = 0, .y = 0 },
};

/// A path tree without per-shape values.
///
/// `init(gpa, paths)` copies what it needs from the paths; ids reported by
/// queries are indices into them.
pub const PathTree = FoldedPaths(void, struct {
    fn combine(_: void, _: void) void {}
}.combine);

/// A path tree that also carries a value of type `T` per shape and can combine
/// the values of all shapes containing a point without visiting them.
/// `init(gpa, paths, values)` takes one value per path.
///
/// `combine(a, b)` must be associative. It is always called with the shapes
/// behind `a` at lower indices than those behind `b`, so it need not be
/// commutative: with shapes listed back to front, alpha compositing works.
pub fn FoldedPaths(comptime T: type, comptime combine: fn (T, T) T) type {
    return sweep.Tree(Paths, T, combine);
}

/// The polar form of a cubic: `blossom(p, t, t, t)` is the curve at t, and the
/// control points of the piece over [s, t] are the blossoms of (s, s, s),
/// (s, s, t), (s, t, t) and (t, t, t).
fn blossom(p: [4]f64, a: f64, b: f64, c: f64) f64 {
    const q0 = p[0] + a * (p[1] - p[0]);
    const q1 = p[1] + a * (p[2] - p[1]);
    const q2 = p[2] + a * (p[3] - p[2]);
    const r0 = q0 + b * (q1 - q0);
    const r1 = q1 + b * (q2 - q1);
    return r0 + c * (r1 - r0);
}

fn between(p: [4]f64, s: f64, t: f64) [4]f64 {
    return .{ blossom(p, s, s, s), blossom(p, s, s, t), blossom(p, s, t, t), blossom(p, t, t, t) };
}

/// Appends the roots of a t^2 + b t + c that lie strictly inside (0, 1).
fn appendUnitRoots(roots: *std.ArrayList(f64), a: f64, b: f64, c: f64) void {
    var found: [2]f64 = undefined;
    var n: usize = 0;
    if (a == 0) {
        if (b != 0) {
            found[0] = -c / b;
            n = 1;
        }
    } else {
        const disc = b * b - 4 * a * c;
        if (disc >= 0) {
            const q = -0.5 * (b + std.math.copysign(@sqrt(disc), b));
            found[0] = q / a;
            n = 1;
            if (q != 0) {
                found[1] = c / q;
                n = 2;
            }
        }
    }
    for (found[0..n]) |t| {
        if (t > 1e-9 and t < 1 - 1e-9) roots.appendAssumeCapacity(t);
    }
}

/// A stretch of a cubic that runs left to right and does not inflect, so that
/// as a function of x its height is single-valued and its slope is monotone.
const Piece = struct {
    /// Control points.
    x: [4]f64,
    y: [4]f64,
    /// The same cubics as polynomials in t, lowest power first.
    px: [4]f64,
    py: [4]f64,

    fn init(x: [4]f64, y: [4]f64) Piece {
        return .{ .x = x, .y = y, .px = power(x), .py = power(y) };
    }

    fn power(c: [4]f64) [4]f64 {
        return .{ c[0], 3 * (c[1] - c[0]), 3 * (c[2] - 2 * c[1] + c[0]), c[3] - 3 * c[2] + 3 * c[1] - c[0] };
    }

    fn eval(p: [4]f64, t: f64) f64 {
        return ((p[3] * t + p[2]) * t + p[1]) * t + p[0];
    }

    fn slopeOf(p: [4]f64, t: f64) f64 {
        return (3 * p[3] * t + 2 * p[2]) * t + p[1];
    }

    const Sample = struct {
        t: f64,
        y: f64,
        /// dy/dx, or NaN where the piece is too steep to tell.
        slope: f64,
    };

    /// The parameter at which the piece is at `x`, clamped to its extent.
    fn tAt(p: *const Piece, x: f64) f64 {
        if (x <= p.x[0]) return 0;
        if (x >= p.x[3]) return 1;
        var lo: f64 = 0;
        var hi: f64 = 1;
        var t = (x - p.x[0]) / (p.x[3] - p.x[0]);
        // Below this the residual is rounding error in evaluating the cubic.
        const noise = 8 * std.math.floatEps(f64) * @max(@abs(p.x[0]), @abs(p.x[3]));
        // Halley's method, falling back to bisection when a step leaves the
        // bracket.
        for (0..64) |_| {
            const f = eval(p.px, t) - x;
            if (@abs(f) <= noise) break;
            if (f > 0) hi = t else lo = t;
            const d1 = slopeOf(p.px, t);
            const d2 = 6 * p.px[3] * t + 2 * p.px[2];
            var next = t - 2 * f * d1 / (2 * d1 * d1 - f * d2);
            if (!(next > lo and next < hi)) next = lo + (hi - lo) / 2;
            const step = @abs(next - t);
            t = next;
            if (step <= 4 * std.math.floatEps(f64)) break;
        }
        return t;
    }

    fn yAt(p: *const Piece, x: f64) f64 {
        return eval(p.py, p.tAt(x));
    }

    fn sample(p: *const Piece, x: f64) Sample {
        const t = p.tAt(x);
        // Where the piece is nearly vertical, rounding decides the sign.
        const dx = slopeOf(p.px, t);
        const slope = if (dx > 1e-9 * (p.x[3] - p.x[0])) slopeOf(p.py, t) / dx else std.math.nan(f64);
        return .{ .t = t, .y = eval(p.py, t), .slope = slope };
    }

    /// Whether some line parallel to this piece's chord has the two pieces
    /// strictly on opposite sides. A piece lies within its control points.
    fn chordSeparates(p: *const Piece, other: *const Piece) bool {
        const nx = p.y[0] - p.y[3];
        const ny = p.x[3] - p.x[0];
        var min = std.math.inf(f64);
        var max = -std.math.inf(f64);
        var other_min = std.math.inf(f64);
        var other_max = -std.math.inf(f64);
        for (0..4) |i| {
            const d = nx * p.x[i] + ny * p.y[i];
            min = @min(min, d);
            max = @max(max, d);
            const e = nx * other.x[i] + ny * other.y[i];
            other_min = @min(other_min, e);
            other_max = @max(other_max, e);
        }
        return max < other_min or other_max < min;
    }

    /// Bounds on the height anywhere on the piece.
    fn yBounds(p: *const Piece) [2]f64 {
        return .{ @min(p.y[0], p.y[1], p.y[2], p.y[3]), @max(p.y[0], p.y[1], p.y[2], p.y[3]) };
    }

    /// Bounds on the height between two parameters.
    fn yRange(p: *const Piece, s: f64, t: f64) [2]f64 {
        const c = between(p.y, s, t);
        return .{ @min(c[0], c[1], c[2], c[3]), @max(c[0], c[1], c[2], c[3]) };
    }
};

/// A maximal run of pieces going the same way along a contour, stored left to
/// right. Every arc starts at a leftmost point of its contour, together with
/// the arc on the other side of that point, and ends at a rightmost point the
/// same way.
const Arc = struct {
    first: u32,
    len: u32,
    shape: u32,
};

const Paths = struct {
    pieces: []const Piece,
    arcs: []const Arc,
    /// The arcs of shape i are `shape_arcs[i]..shape_arcs[i + 1]`.
    shape_arcs: []const u32,
    /// Crossing hints are within this distance of the crossing in x.
    tolerance: f64,

    pub const Input = []const Path;

    pub fn init(gpa: Allocator, paths: []const Path, events: *std.ArrayList(Event)) Allocator.Error!Paths {
        var builder: Builder = .{ .gpa = gpa, .events = events };
        defer builder.deinit();

        var scale: f64 = 0;
        for (paths) |path| {
            for (path.curves) |c| {
                for ([_]Vec2{ c.p0, c.p1, c.p2, c.p3 }) |p| {
                    scale = @max(scale, @abs(@as(f64, p.x) + path.offset.x), @abs(@as(f64, p.y) + path.offset.y));
                }
            }
        }
        builder.tolerance = 1e-9 * scale;

        const shape_arcs = try gpa.alloc(u32, paths.len + 1);
        errdefer gpa.free(shape_arcs);
        for (paths, 0..) |path, shape| {
            shape_arcs[shape] = @intCast(builder.arcs.items.len);
            var begin: usize = 0;
            for (1..path.curves.len + 1) |end| {
                if (end < path.curves.len and std.meta.eql(path.curves[end].p0, path.curves[end - 1].p3)) continue;
                try builder.addContour(@intCast(shape), path.curves[begin..end], path.offset);
                begin = end;
            }
        }
        shape_arcs[paths.len] = @intCast(builder.arcs.items.len);

        try builder.findCrossings();

        const pieces = try builder.pieces.toOwnedSlice(gpa);
        errdefer gpa.free(pieces);
        const arcs = try builder.arcs.toOwnedSlice(gpa);
        return .{ .pieces = pieces, .arcs = arcs, .shape_arcs = shape_arcs, .tolerance = builder.tolerance };
    }

    pub fn deinit(self: *Paths, gpa: Allocator) void {
        gpa.free(self.pieces);
        gpa.free(self.arcs);
        gpa.free(self.shape_arcs);
    }

    pub fn shapeCount(self: *const Paths) usize {
        return self.shape_arcs.len - 1;
    }

    pub fn arcCount(self: *const Paths) usize {
        return self.arcs.len;
    }

    pub fn shapeOf(self: *const Paths, arc: u32) u32 {
        return self.arcs[arc].shape;
    }

    pub fn arcY(self: *const Paths, arc: u32, x: f64) f64 {
        return arcPiece(self.pieces, self.arcs[arc], x).yAt(x);
    }

    pub fn arcBounds(self: *const Paths, arc: u32, x: f64) [2]f64 {
        return arcPiece(self.pieces, self.arcs[arc], x).yBounds();
    }

    pub fn slack(self: *const Paths, _: u32) f64 {
        return self.tolerance;
    }

    pub fn bytes(self: *const Paths) usize {
        return std.mem.sliceAsBytes(self.pieces).len + std.mem.sliceAsBytes(self.arcs).len +
            std.mem.sliceAsBytes(self.shape_arcs).len;
    }

    /// Whether (x, y) is inside `shape`, by looking at every arc of the shape.
    pub fn contains(self: *const Paths, shape: usize, x: f64, y: f64) bool {
        var inside = false;
        for (self.arcs[self.shape_arcs[shape]..self.shape_arcs[shape + 1]]) |arc| {
            const first = self.pieces[arc.first];
            const last = self.pieces[arc.first + arc.len - 1];
            if (x < first.x[0] or x >= last.x[3]) continue;
            if (arcPiece(self.pieces, arc, x).yAt(x) <= y) inside = !inside;
        }
        return inside;
    }
};

fn arcPiece(pieces: []const Piece, arc: Arc, x: f64) *const Piece {
    var i = arc.first;
    const last = arc.first + arc.len - 1;
    while (i < last and x > pieces[i].x[3]) i += 1;
    return &pieces[i];
}

const Builder = struct {
    gpa: Allocator,
    events: *std.ArrayList(Event),
    tolerance: f64 = 0,
    pieces: std.ArrayList(Piece) = .empty,
    arcs: std.ArrayList(Arc) = .empty,
    /// For each arc, the arcs it meets at its left and right ends.
    partners: std.ArrayList([2]u32) = .empty,
    /// Pieces of the contour being added, in drawing order.
    drawn: std.ArrayList(Drawn) = .empty,

    const Drawn = struct {
        piece: Piece,
        rightward: bool,
    };

    fn deinit(self: *Builder) void {
        self.pieces.deinit(self.gpa);
        self.arcs.deinit(self.gpa);
        self.partners.deinit(self.gpa);
        self.drawn.deinit(self.gpa);
    }

    fn addContour(self: *Builder, shape: u32, curves: []const CubicBezier, offset: Vec2) Allocator.Error!void {
        self.drawn.clearRetainingCapacity();
        for (curves) |c| {
            try self.draw(
                .{ c.p0.x, c.p1.x, c.p2.x, c.p3.x },
                .{ c.p0.y, c.p1.y, c.p2.y, c.p3.y },
                offset,
            );
        }
        const from = curves[curves.len - 1].p3;
        const to = curves[0].p0;
        if (!std.meta.eql(from, to)) {
            const a = from.lerp(to, 1.0 / 3.0);
            const b = from.lerp(to, 2.0 / 3.0);
            try self.draw(.{ from.x, a.x, b.x, to.x }, .{ from.y, a.y, b.y, to.y }, offset);
        }

        // Rotate so that the contour is a whole number of runs.
        const drawn = self.drawn.items;
        const n = drawn.len;
        var begin: usize = 0;
        while (begin < n and drawn[begin].rightward == drawn[(begin + n - 1) % n].rightward) begin += 1;
        if (begin >= n) return;

        const first_arc: u32 = @intCast(self.arcs.items.len);
        var i: usize = 0;
        while (i < n) {
            const rightward = drawn[(begin + i) % n].rightward;
            var len: usize = 1;
            while (i + len < n and drawn[(begin + i + len) % n].rightward == rightward) len += 1;

            try self.arcs.append(self.gpa, .{
                .first = @intCast(self.pieces.items.len),
                .len = @intCast(len),
                .shape = shape,
            });
            for (0..len) |k| {
                const at = if (rightward) i + k else i + len - 1 - k;
                try self.pieces.append(self.gpa, drawn[(begin + at) % n].piece);
            }
            i += len;
        }

        // Runs alternate direction, so neighbours meet at an extreme point.
        const arc_count: u32 = @intCast(self.arcs.items.len - first_arc);
        try self.partners.appendNTimes(self.gpa, .{ sweep.nil, sweep.nil }, arc_count);
        for (0..arc_count) |k| {
            const a = first_arc + @as(u32, @intCast(k));
            const b = first_arc + @as(u32, @intCast((k + 1) % arc_count));
            const a_rightward = drawn[begin].rightward == (k % 2 == 0);
            const side: usize = @intFromBool(a_rightward);
            self.partners.items[a][side] = b;
            self.partners.items[b][side] = a;
            try self.events.append(self.gpa, .{
                .x = if (a_rightward) @max(self.arcEnd(a), self.arcEnd(b)) else @min(self.arcStart(a), self.arcStart(b)),
                .a = a,
                .b = @intCast(b),
                .kind = if (a_rightward) .end else .start,
            });
        }
    }

    /// Cuts a cubic where it turns back in x and where it inflects, and keeps
    /// the pieces that are not vertical.
    fn draw(self: *Builder, xs: [4]f32, ys: [4]f32, offset: Vec2) Allocator.Error!void {
        var x: [4]f64 = undefined;
        var y: [4]f64 = undefined;
        for (0..4) |i| {
            x[i] = @as(f64, xs[i]) + offset.x;
            y[i] = @as(f64, ys[i]) + offset.y;
        }

        // The derivative is 3 (a t^2 + b t + c) in each coordinate.
        const ax = (x[1] - x[0]) - 2 * (x[2] - x[1]) + (x[3] - x[2]);
        const bx = 2 * ((x[2] - x[1]) - (x[1] - x[0]));
        const cx = x[1] - x[0];
        const ay = (y[1] - y[0]) - 2 * (y[2] - y[1]) + (y[3] - y[2]);
        const by = 2 * ((y[2] - y[1]) - (y[1] - y[0]));
        const cy = y[1] - y[0];

        var buffer: [6]f64 = undefined;
        var cuts: std.ArrayList(f64) = .initBuffer(&buffer);
        cuts.appendAssumeCapacity(0);
        appendUnitRoots(&cuts, ax, bx, cx);
        // Inflections: where the first and second derivatives are parallel.
        appendUnitRoots(&cuts, bx * ay - ax * by, 2 * (cx * ay - cy * ax), cx * by - cy * bx);
        cuts.appendAssumeCapacity(1);
        std.mem.sortUnstable(f64, cuts.items, {}, std.sort.asc(f64));

        for (cuts.items[0 .. cuts.items.len - 1], cuts.items[1..]) |s, t| {
            if (s == t) continue;
            var px = between(x, s, t);
            var py = between(y, s, t);
            if (px[0] == px[3]) continue;
            const rightward = px[0] < px[3];
            if (!rightward) {
                std.mem.reverse(f64, &px);
                std.mem.reverse(f64, &py);
            }
            try self.drawn.append(self.gpa, .{ .piece = .init(px, py), .rightward = rightward });
        }
    }

    fn arcStart(self: *const Builder, arc: u32) f64 {
        return self.pieces.items[self.arcs.items[arc].first].x[0];
    }

    fn arcEnd(self: *const Builder, arc: u32) f64 {
        const a = self.arcs.items[arc];
        return self.pieces.items[a.first + a.len - 1].x[3];
    }

    /// Emits a cross event near every point where two arcs cross.
    fn findCrossings(self: *Builder) Allocator.Error!void {
        const boxes = try self.gpa.alloc(sweep.Box, self.arcs.items.len);
        defer self.gpa.free(boxes);
        for (boxes, self.arcs.items, 0..) |*box, arc, i| {
            var bottom = std.math.inf(f64);
            var top = -std.math.inf(f64);
            for (self.pieces.items[arc.first..][0..arc.len]) |*piece| {
                const bounds = piece.yBounds();
                bottom = @min(bottom, bounds[0]);
                top = @max(top, bounds[1]);
            }
            const id: u32 = @intCast(i);
            box.* = .{ .x = .{ self.arcStart(id), self.arcEnd(id) }, .y = .{ bottom, top }, .id = id };
        }
        try sweep.forEachOverlap(self.gpa, boxes, self);
    }

    /// Emits a cross event near every point where arcs `a` and `b` cross.
    pub fn pair(self: *Builder, a: u32, b: u32) Allocator.Error!void {
        var lo = @max(self.arcStart(a), self.arcStart(b));
        var hi = @min(self.arcEnd(a), self.arcEnd(b));
        // Arcs that meet at an end are handled there by their start or end
        // event, and would otherwise look like a crossing.
        if (self.partners.items[a][0] == b) lo += self.tolerance;
        if (self.partners.items[a][1] == b) hi -= self.tolerance;
        if (!(lo <= hi)) return;

        const arc_a = self.arcs.items[a];
        const arc_b = self.arcs.items[b];
        try self.crossJumps(a, b);
        try self.crossJumps(b, a);
        for (self.pieces.items[arc_a.first..][0..arc_a.len]) |*pa| {
            for (self.pieces.items[arc_b.first..][0..arc_b.len]) |*pb| {
                const u = @max(lo, pa.x[0], pb.x[0]);
                const v = @min(hi, pa.x[3], pb.x[3]);
                if (!(u <= v)) continue;
                const ra = pa.yBounds();
                const rb = pb.yBounds();
                if (ra[1] < rb[0] or rb[1] < ra[0]) continue;
                if (pa.chordSeparates(pb) or pb.chordSeparates(pa)) continue;
                var search: Search = .{ .builder = self, .a = a, .b = b, .pa = pa, .pb = pb };
                if (u == v or std.meta.eql(pa.*, pb.*)) {
                    try search.hint(u);
                    try search.hint(v);
                } else {
                    try search.run(u, v, pa.sample(u), pa.sample(v), pb.sample(u), pb.sample(v));
                }
            }
        }
    }

    /// A vertical edge inside arc `a` makes it jump, passing whatever lies in
    /// between without meeting it.
    fn crossJumps(self: *Builder, a: u32, b: u32) Allocator.Error!void {
        const arc = self.arcs.items[a];
        const pieces = self.pieces.items[arc.first..][0..arc.len];
        for (pieces[0 .. pieces.len - 1], pieces[1..]) |before, after| {
            if (before.y[3] == after.y[0]) continue;
            const x = before.x[3];
            if (x < self.arcStart(b) or x > self.arcEnd(b)) continue;
            try self.events.append(self.gpa, .{ .x = x, .a = a, .b = @intCast(b), .kind = .cross });
        }
    }

    /// Looks for the crossings of two pieces over a range of x by subdivision.
    const Search = struct {
        builder: *Builder,
        a: u32,
        b: u32,
        pa: *const Piece,
        pb: *const Piece,
        /// Subdivisions left before giving up on telling near-coincident
        /// pieces apart.
        budget: u32 = 256,

        fn hint(self: *Search, x: f64) Allocator.Error!void {
            try self.builder.events.append(self.builder.gpa, .{ .x = x, .a = self.a, .b = @intCast(self.b), .kind = .cross });
        }

        fn run(self: *Search, u: f64, v: f64, au: Piece.Sample, av: Piece.Sample, bu: Piece.Sample, bv: Piece.Sample) Allocator.Error!void {
            const gu = au.y - bu.y;
            const gv = av.y - bv.y;

            // Pieces that leave a shared point along the same tangent cannot
            // be told apart until they have separated, which may be well past
            // the point. Their range is halved all the way down, with a hint
            // at every cut, so that they are re-examined once they can be.
            const touching = 1e-2 * self.builder.tolerance;
            const kiss_u = @abs(gu) <= touching and !(@abs(au.slope - bu.slope) > 1e-3);
            const kiss_v = @abs(gv) <= touching and !(@abs(av.slope - bv.slope) > 1e-3);
            const kiss = kiss_u or kiss_v;
            if (@abs(gu) <= touching) try self.hint(u);
            if (@abs(gv) <= touching) try self.hint(v);

            if (!kiss) {
                const ra = self.pa.yRange(au.t, av.t);
                const rb = self.pb.yRange(bu.t, bv.t);
                if (ra[1] < rb[0] or rb[1] < ra[0]) return;

                // Neither piece inflects, so its slope over the range lies
                // between its slopes at the ends. That bounds the slope of
                // the gap.
                const steepest = @max(au.slope, av.slope) - @min(bu.slope, bv.slope);
                const shallowest = @min(au.slope, av.slope) - @max(bu.slope, bv.slope);
                const known = !std.math.isNan(au.slope + av.slope + bu.slope + bv.slope) and
                    !std.math.isNan(steepest) and !std.math.isNan(shallowest);
                if (known) {
                    if (shallowest > 0 or steepest < 0) {
                        // The gap is monotone: one crossing or none.
                        if (gu != 0 and gv != 0 and (gu < 0) != (gv < 0)) try self.hint(self.root(u, v, gu));
                        return;
                    }
                    // The gap cannot close faster than its steepest slope.
                    const reach = @max(@abs(steepest), @abs(shallowest)) * (v - u);
                    if ((gu < 0) == (gv < 0) and @abs(gu + gv) > reach) return;
                }
            }

            if (v - u <= self.builder.tolerance or self.budget == 0) {
                try self.hint(u);
                try self.hint(v);
                return;
            }
            self.budget -= 1;
            const mid = u + (v - u) / 2;
            const am = self.pa.sample(mid);
            const bm = self.pb.sample(mid);
            if (kiss) {
                // Still together at both ends and in the middle: nothing in
                // this range can be told apart.
                if (@abs(gu) <= touching and @abs(gv) <= touching and @abs(am.y - bm.y) <= touching) return;
                try self.hint(mid);
            }
            try self.run(u, mid, au, am, bu, bm);
            try self.run(mid, v, am, av, bm, bv);
        }

        /// The crossing inside a range whose ends have gaps of opposite sign,
        /// by Newton's method, bisecting when a step leaves the bracket.
        fn root(self: *Search, u: f64, v: f64, gu: f64) f64 {
            var lo = u;
            var hi = v;
            var x = u + (v - u) / 2;
            for (0..64) |_| {
                const sa = self.pa.sample(x);
                const sb = self.pb.sample(x);
                const g = sa.y - sb.y;
                if (g == 0) break;
                if ((g < 0) == (gu < 0)) lo = x else hi = x;
                var next = x - g / (sa.slope - sb.slope);
                if (!(next > lo and next < hi)) next = lo + (hi - lo) / 2;
                const step = @abs(next - x);
                x = next;
                if (step <= 0.1 * self.builder.tolerance) break;
            }
            return x;
        }
    };
};

const testing = std.testing;
const shapes = @import("geometry/shapes.zig");

fn line(ax: f32, ay: f32, bx: f32, by: f32) CubicBezier {
    const a = Vec2.init(ax, ay);
    const b = Vec2.init(bx, by);
    return CubicBezier.init(a, a.lerp(b, 1.0 / 3.0), a.lerp(b, 2.0 / 3.0), b);
}

/// A closed polygon through the given points.
fn polygon(gpa: Allocator, curves: *std.ArrayList(CubicBezier), points: []const [2]f32) !void {
    for (points, 0..) |p, i| {
        const q = points[(i + 1) % points.len];
        try curves.append(gpa, line(p[0], p[1], q[0], q[1]));
    }
}

/// Checks `paths` against a scan of every arc at random points in `box` (min x,
/// min y, max x, max y), skipping points within `margin` of a boundary.
fn expectMatchesScan(paths: []const Path, box: [4]f64, queries: usize, margin: f64, random: std.Random) !void {
    const gpa = testing.allocator;
    var tree = try PathTree.init(gpa, paths);
    defer tree.deinit(gpa);
    const geometry = &tree.geometry;

    var got: std.ArrayList(u32) = .empty;
    defer got.deinit(gpa);
    var want: std.ArrayList(u32) = .empty;
    defer want.deinit(gpa);

    var done: usize = 0;
    points: while (done < queries) {
        const x = box[0] + random.float(f64) * (box[2] - box[0]);
        const y = box[1] + random.float(f64) * (box[3] - box[1]);

        for (geometry.arcs, 0..) |arc, i| {
            const first = geometry.pieces[arc.first];
            const last = geometry.pieces[arc.first + arc.len - 1];
            if (x < first.x[0] - margin or x > last.x[3] + margin) continue;
            if (@abs(x - first.x[0]) < margin or @abs(x - last.x[3]) < margin) continue :points;
            if (@abs(geometry.arcY(@intCast(i), x) - y) < margin) continue :points;
        }
        want.clearRetainingCapacity();
        for (0..paths.len) |shape| {
            if (geometry.contains(shape, x, y)) try want.append(gpa, @intCast(shape));
        }
        done += 1;

        got.clearRetainingCapacity();
        try tree.queryAppend(gpa, &got, x, y);
        try testing.expectEqualSlices(u32, want.items, got.items);
    }
}

test "no paths" {
    var tree = try PathTree.init(testing.allocator, &.{});
    defer tree.deinit(testing.allocator);
    const got = try tree.query(testing.allocator, 0, 0);
    defer testing.allocator.free(got);
    try testing.expectEqual(0, got.len);
}

test "square with a hole, by the even-odd rule" {
    const gpa = testing.allocator;
    var curves: std.ArrayList(CubicBezier) = .empty;
    defer curves.deinit(gpa);
    try polygon(gpa, &curves, &.{ .{ 0, 0 }, .{ 10, 0 }, .{ 10, 10 }, .{ 0, 10 } });
    try polygon(gpa, &curves, &.{ .{ 3, 3 }, .{ 7, 3 }, .{ 7, 7 }, .{ 3, 7 } });

    var tree = try PathTree.init(gpa, &.{.{ .curves = curves.items }});
    defer tree.deinit(gpa);

    for ([_]struct { f64, f64, usize }{
        .{ 1, 1, 1 },   .{ 5, 5, 0 },  .{ 5, 8.5, 1 }, .{ 5, 1.5, 1 },
        .{ 8.5, 5, 1 }, .{ 11, 5, 0 }, .{ 5, 11, 0 },  .{ -1, 5, 0 },
    }) |case| {
        const got = try tree.query(gpa, case[0], case[1]);
        defer gpa.free(got);
        try testing.expectEqual(case[2], got.len);
    }
}

test "agrees with ConvexBezierShape" {
    const gpa = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(1);
    const random = prng.random();

    var convex: [60]shapes.ConvexBezierShape = undefined;
    var paths: [convex.len]Path = undefined;
    for (&convex, &paths) |*shape, *path| {
        shape.* = try shapes.createRandomConvexShape(
            gpa,
            random,
            random.float(f32) * 100,
            random.float(f32) * 100,
            2 + random.float(f32) * 18,
            3 + random.uintLessThan(u32, 6),
        );
        shape.translate(random.float(f32) * 4, random.float(f32) * 4);
        path.* = .{ .curves = shape.curves.items, .offset = shape.offset };
    }
    defer for (&convex) |*shape| shape.deinit();

    var tree = try PathTree.init(gpa, &paths);
    defer tree.deinit(gpa);

    var checked: usize = 0;
    points: while (checked < 5000) {
        const x = random.float(f32) * 110 - 5;
        const y = random.float(f32) * 110 - 5;
        // The shape's own test is in single precision with a tolerance, so
        // only points clearly on one side of every boundary are compared.
        var want: usize = 0;
        for (&convex) |*shape| {
            const interval = shape.yIntervalAtX(x) orelse {
                if (@abs(x - shape.bbox.x_min) < 0.05 or @abs(x - shape.bbox.x_max) < 0.05) continue :points;
                continue;
            };
            if (@abs(y - interval.y_min) < 0.05 or @abs(y - interval.y_max) < 0.05) continue :points;
            if (@abs(x - shape.bbox.x_min) < 0.05 or @abs(x - shape.bbox.x_max) < 0.05) continue :points;
            want += @intFromBool(y > interval.y_min and y < interval.y_max);
        }
        checked += 1;

        const got = try tree.query(gpa, x, y);
        defer gpa.free(got);
        try testing.expectEqual(want, got.len);
    }
}

test "random convex shapes" {
    const gpa = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(2);
    const random = prng.random();

    for ([_]f32{ 2, 8, 30 }) |size| {
        var convex: [150]shapes.ConvexBezierShape = undefined;
        var paths: [convex.len]Path = undefined;
        for (&convex, &paths) |*shape, *path| {
            shape.* = try shapes.createRandomConvexShape(
                gpa,
                random,
                random.float(f32) * 100,
                random.float(f32) * 100,
                size * (0.2 + random.float(f32)),
                3 + random.uintLessThan(u32, 8),
            );
            path.* = .{ .curves = shape.curves.items };
        }
        defer for (&convex) |*shape| shape.deinit();
        try expectMatchesScan(&paths, .{ -40, -40, 140, 140 }, 4000, 1e-7, random);
    }
}

test "random curved loops that cross themselves" {
    const gpa = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(3);
    const random = prng.random();

    var curves: std.ArrayList(CubicBezier) = .empty;
    defer curves.deinit(gpa);
    var lens: [40]usize = undefined;
    for (&lens) |*len| {
        // A closed chain of arbitrary cubics.
        const n = 2 + random.uintLessThan(usize, 5);
        const cx = random.float(f32) * 100;
        const cy = random.float(f32) * 100;
        var first: Vec2 = undefined;
        var from: Vec2 = undefined;
        for (0..n) |i| {
            const point = struct {
                fn at(r: std.Random, x: f32, y: f32) Vec2 {
                    return Vec2.init(x + (r.float(f32) - 0.5) * 50, y + (r.float(f32) - 0.5) * 50);
                }
            }.at;
            if (i == 0) {
                first = point(random, cx, cy);
                from = first;
            }
            const to = if (i == n - 1) first else point(random, cx, cy);
            try curves.append(gpa, CubicBezier.init(from, point(random, cx, cy), point(random, cx, cy), to));
            from = to;
        }
        len.* = n;
    }

    var paths: [lens.len]Path = undefined;
    var begin: usize = 0;
    for (&paths, lens) |*path, len| {
        path.* = .{ .curves = curves.items[begin..][0..len] };
        begin += len;
    }
    try expectMatchesScan(&paths, .{ -30, -30, 130, 130 }, 20000, 1e-7, random);
}

test "polygons sharing edges and corners" {
    const gpa = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(4);
    const random = prng.random();

    var curves: std.ArrayList(CubicBezier) = .empty;
    defer curves.deinit(gpa);
    // A checkerboard of squares, then triangles and a star over its corners.
    for (0..8) |i| for (0..8) |j| {
        const x: f32 = @floatFromInt(i);
        const y: f32 = @floatFromInt(j);
        try polygon(gpa, &curves, &.{ .{ x, y }, .{ x + 1, y }, .{ x + 1, y + 1 }, .{ x, y + 1 } });
    };
    for (0..7) |i| {
        const x: f32 = @floatFromInt(i);
        try polygon(gpa, &curves, &.{ .{ x, x }, .{ x + 2, x }, .{ x, x + 2 } });
    }
    try polygon(gpa, &curves, &.{ .{ 4, 8 }, .{ 6, 0 }, .{ 0, 5 }, .{ 8, 5 }, .{ 2, 0 } });

    var paths: std.ArrayList(Path) = .empty;
    defer paths.deinit(gpa);
    for (0..64) |i| try paths.append(gpa, .{ .curves = curves.items[4 * i ..][0..4] });
    for (0..7) |i| try paths.append(gpa, .{ .curves = curves.items[256 + 3 * i ..][0..3] });
    try paths.append(gpa, .{ .curves = curves.items[277..] });
    try expectMatchesScan(paths.items, .{ -1, -1, 9, 9 }, 20000, 1e-7, random);
}

test "fold over paths" {
    const gpa = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(5);
    const random = prng.random();

    var convex: [80]shapes.ConvexBezierShape = undefined;
    var paths: [convex.len]Path = undefined;
    var ones: [convex.len]u32 = undefined;
    for (&convex, &paths, &ones) |*shape, *path, *one| {
        shape.* = try shapes.createRandomConvexShape(gpa, random, random.float(f32) * 50, random.float(f32) * 50, 12, 5);
        path.* = .{ .curves = shape.curves.items };
        one.* = 1;
    }
    defer for (&convex) |*shape| shape.deinit();

    const Counting = FoldedPaths(u32, struct {
        fn add(a: u32, b: u32) u32 {
            return a + b;
        }
    }.add);
    var tree = try Counting.init(gpa, &paths, &ones);
    defer tree.deinit(gpa);

    for (0..5000) |_| {
        const x = random.float(f64) * 50;
        const y = random.float(f64) * 50;
        const ids = try tree.query(gpa, x, y);
        defer gpa.free(ids);
        try testing.expectEqual(ids.len, tree.fold(x, y) orelse 0);
    }
}

test "vertical edge inside an arc" {
    const gpa = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(6);
    var curves: std.ArrayList(CubicBezier) = .empty;
    defer curves.deinit(gpa);
    // A staircase whose risers are vertical, and a bar passing through them.
    try polygon(gpa, &curves, &.{ .{ 0, 0 }, .{ 9, 0 }, .{ 9, 9 }, .{ 6, 9 }, .{ 6, 6 }, .{ 3, 6 }, .{ 3, 3 }, .{ 0, 3 } });
    try polygon(gpa, &curves, &.{ .{ -1, 4 }, .{ 10, 5 }, .{ 10, 8 }, .{ -1, 5 } });
    try expectMatchesScan(
        &.{ .{ .curves = curves.items[0..8] }, .{ .curves = curves.items[8..] } },
        .{ -2, -1, 11, 10 },
        20000,
        1e-7,
        prng.random(),
    );
}

test "curves leaving a point along the same tangent" {
    var prng: std.Random.DefaultPrng = .init(7);
    const P = Vec2.init;
    const C = CubicBezier.init;
    // Both shapes pass through (1, 1) with slope 5, and cross again later.
    const a = [_]CubicBezier{
        C(P(1, 11), P(1, 7), P(4, 2), P(11, 10)), C(P(11, 10), P(8, 1), P(3, 0), P(1, 1)),
        C(P(1, 1), P(3, 11), P(8, 8), P(2, 11)),  C(P(2, 11), P(0, 11), P(11, 7), P(7, 0)),
        C(P(7, 0), P(7, 8), P(0, 10), P(1, 11)),
    };
    const b = [_]CubicBezier{
        C(P(6, 4), P(5, 3), P(3, 11), P(1, 1)),   C(P(1, 1), P(4, 8), P(0, 1), P(4, 6)),
        C(P(4, 6), P(5, 4), P(7, 7), P(11, 0)),   C(P(11, 0), P(6, 0), P(10, 8), P(11, 7)),
        C(P(11, 7), P(7, 4), P(10, 10), P(6, 4)),
    };
    try expectMatchesScan(&.{ .{ .curves = &a }, .{ .curves = &b } }, .{ -1, -1, 12, 12 }, 20000, 1e-7, prng.random());
}

test "curves with control points on a lattice" {
    // Cusps, tangencies, vertical curves and shared end points are all common.
    const gpa = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(8);
    const random = prng.random();
    const lattice = struct {
        fn point(r: std.Random) Vec2 {
            return Vec2.init(@floatFromInt(r.uintLessThan(u32, 12)), @floatFromInt(r.uintLessThan(u32, 12)));
        }
    };

    var curves: std.ArrayList(CubicBezier) = .empty;
    defer curves.deinit(gpa);
    var paths: std.ArrayList(Path) = .empty;
    defer paths.deinit(gpa);
    for (0..30) |_| {
        curves.clearRetainingCapacity();
        paths.clearRetainingCapacity();
        var lens: [40]usize = undefined;
        for (&lens) |*len| {
            len.* = 2 + random.uintLessThan(usize, 4);
            const first = lattice.point(random);
            var from = first;
            for (0..len.*) |i| {
                const to = if (i == len.* - 1) first else lattice.point(random);
                try curves.append(gpa, CubicBezier.init(from, lattice.point(random), lattice.point(random), to));
                from = to;
            }
        }
        var begin: usize = 0;
        for (lens) |len| {
            try paths.append(gpa, .{ .curves = curves.items[begin..][0..len] });
            begin += len;
        }
        try expectMatchesScan(paths.items, .{ -1, -1, 12, 12 }, 3000, 1e-7, random);
    }
}
