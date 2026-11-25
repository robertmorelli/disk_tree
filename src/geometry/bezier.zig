//! ALGORITHMS FIRST: slab-decomposed segment-tree renderer.
//! Do not replace the algorithm with brute-force point checks.
//! Preserve algorithmic structure; micro-optimizations come second.

const std = @import("std");
const math = std.math;
const cubic = @import("../math/cubic.zig");
const Allocator = std.mem.Allocator;

/// A 2D point/vector
pub const Vec2 = struct {
    x: f32,
    y: f32,

    pub fn init(x: f32, y: f32) Vec2 {
        return .{ .x = x, .y = y };
    }

    pub fn add(self: Vec2, other: Vec2) Vec2 {
        return .{ .x = self.x + other.x, .y = self.y + other.y };
    }

    pub fn sub(self: Vec2, other: Vec2) Vec2 {
        return .{ .x = self.x - other.x, .y = self.y - other.y };
    }

    pub fn scale(self: Vec2, s: f32) Vec2 {
        return .{ .x = self.x * s, .y = self.y * s };
    }

    pub fn lerp(self: Vec2, other: Vec2, t: f32) Vec2 {
        return self.scale(1 - t).add(other.scale(t));
    }

    pub fn length(self: Vec2) f32 {
        return @sqrt(self.x * self.x + self.y * self.y);
    }
};

/// Axis-aligned bounding box
pub const BoundingBox = struct {
    x_min: f32,
    y_min: f32,
    x_max: f32,
    y_max: f32,

    pub fn init(x_min: f32, y_min: f32, x_max: f32, y_max: f32) BoundingBox {
        return .{ .x_min = x_min, .y_min = y_min, .x_max = x_max, .y_max = y_max };
    }

    pub fn empty() BoundingBox {
        return .{
            .x_min = math.floatMax(f32),
            .y_min = math.floatMax(f32),
            .x_max = -math.floatMax(f32),
            .y_max = -math.floatMax(f32),
        };
    }

    pub fn expand(self: *BoundingBox, x: f32, y: f32) void {
        self.x_min = @min(self.x_min, x);
        self.y_min = @min(self.y_min, y);
        self.x_max = @max(self.x_max, x);
        self.y_max = @max(self.y_max, y);
    }

    pub fn expandBox(self: *BoundingBox, other: BoundingBox) void {
        self.x_min = @min(self.x_min, other.x_min);
        self.y_min = @min(self.y_min, other.y_min);
        self.x_max = @max(self.x_max, other.x_max);
        self.y_max = @max(self.y_max, other.y_max);
    }

    pub fn overlaps(self: BoundingBox, other: BoundingBox) bool {
        return !(self.x_max < other.x_min or self.x_min > other.x_max or
            self.y_max < other.y_min or self.y_min > other.y_max);
    }

    pub fn containsPoint(self: BoundingBox, x: f32, y: f32) bool {
        return x >= self.x_min and x <= self.x_max and
            y >= self.y_min and y <= self.y_max;
    }
};

/// Polynomial coefficients for x(t) or y(t) = at^3 + bt^2 + ct + d
pub const Coefficients = struct {
    a: f32,
    b: f32,
    c: f32,
    d: f32,
};

/// A cubic Bezier curve segment: P(t) = (1-t)^3*P0 + 3(1-t)^2*t*P1 + 3(1-t)*t^2*P2 + t^3*P3
pub const CubicBezier = struct {
    p0: Vec2, // Start point
    p1: Vec2, // Control point 1
    p2: Vec2, // Control point 2
    p3: Vec2, // End point

    pub fn init(p0: Vec2, p1: Vec2, p2: Vec2, p3: Vec2) CubicBezier {
        return .{ .p0 = p0, .p1 = p1, .p2 = p2, .p3 = p3 };
    }

    /// Evaluate the curve at parameter t in [0, 1]
    pub fn evaluate(self: CubicBezier, t: f32) Vec2 {
        const t2 = t * t;
        const t3 = t2 * t;
        const mt = 1 - t;
        const mt2 = mt * mt;
        const mt3 = mt2 * mt;

        return Vec2.init(
            mt3 * self.p0.x + 3 * mt2 * t * self.p1.x + 3 * mt * t2 * self.p2.x + t3 * self.p3.x,
            mt3 * self.p0.y + 3 * mt2 * t * self.p1.y + 3 * mt * t2 * self.p2.y + t3 * self.p3.y,
        );
    }

    /// Evaluate the derivative at parameter t
    pub fn derivative(self: CubicBezier, t: f32) Vec2 {
        const t2 = t * t;
        const mt = 1 - t;
        const mt2 = mt * mt;

        const d01 = self.p1.sub(self.p0);
        const d12 = self.p2.sub(self.p1);
        const d23 = self.p3.sub(self.p2);

        return Vec2.init(
            3 * (mt2 * d01.x + 2 * mt * t * d12.x + t2 * d23.x),
            3 * (mt2 * d01.y + 2 * mt * t * d12.y + t2 * d23.y),
        );
    }

    /// Get coefficients for x(t) = at^3 + bt^2 + ct + d
    pub fn xCoefficients(self: CubicBezier) Coefficients {
        const x0 = self.p0.x;
        const x1 = self.p1.x;
        const x2 = self.p2.x;
        const x3 = self.p3.x;
        return .{
            .a = -x0 + 3 * x1 - 3 * x2 + x3,
            .b = 3 * x0 - 6 * x1 + 3 * x2,
            .c = -3 * x0 + 3 * x1,
            .d = x0,
        };
    }

    /// Get coefficients for y(t) = at^3 + bt^2 + ct + d
    pub fn yCoefficients(self: CubicBezier) Coefficients {
        const y0 = self.p0.y;
        const y1 = self.p1.y;
        const y2 = self.p2.y;
        const y3 = self.p3.y;
        return .{
            .a = -y0 + 3 * y1 - 3 * y2 + y3,
            .b = 3 * y0 - 6 * y1 + 3 * y2,
            .c = -3 * y0 + 3 * y1,
            .d = y0,
        };
    }

    pub const ExtremaResult = struct { values: [2]f32, count: u8 };

    /// Find t values where dx/dt = 0 (x-extrema), in (0, 1)
    pub fn findXExtrema(self: CubicBezier) ExtremaResult {
        const coeffs = self.xCoefficients();
        // dx/dt = 3at^2 + 2bt + c = 0
        const roots = cubic.solveCubic(0, 3 * coeffs.a, 2 * coeffs.b, coeffs.c);

        var result = ExtremaResult{ .values = .{ 0, 0 }, .count = 0 };
        for (roots.slice()) |t| {
            if (t > 0.0 and t < 1.0) {
                if (result.count < 2) {
                    result.values[result.count] = t;
                    result.count += 1;
                }
            }
        }
        return result;
    }

    /// Find t values where dy/dt = 0 (y-extrema), in (0, 1)
    pub fn findYExtrema(self: CubicBezier) ExtremaResult {
        const coeffs = self.yCoefficients();
        const roots = cubic.solveCubic(0, 3 * coeffs.a, 2 * coeffs.b, coeffs.c);

        var result = ExtremaResult{ .values = .{ 0, 0 }, .count = 0 };
        for (roots.slice()) |t| {
            if (t > 0.0 and t < 1.0) {
                if (result.count < 2) {
                    result.values[result.count] = t;
                    result.count += 1;
                }
            }
        }
        return result;
    }

    /// Split the curve at parameter t using de Casteljau's algorithm
    pub fn splitAt(self: CubicBezier, t: f32) struct { left: CubicBezier, right: CubicBezier } {
        const p01 = self.p0.lerp(self.p1, t);
        const p12 = self.p1.lerp(self.p2, t);
        const p23 = self.p2.lerp(self.p3, t);
        const p012 = p01.lerp(p12, t);
        const p123 = p12.lerp(p23, t);
        const p0123 = p012.lerp(p123, t);

        return .{
            .left = CubicBezier.init(self.p0, p01, p012, p0123),
            .right = CubicBezier.init(p0123, p123, p23, self.p3),
        };
    }

    /// Compute tight bounding box
    pub fn boundingBox(self: CubicBezier) BoundingBox {
        var bbox = BoundingBox.empty();

        // Start with endpoints
        bbox.expand(self.p0.x, self.p0.y);
        bbox.expand(self.p3.x, self.p3.y);

        // Add x-extrema
        const x_ext = self.findXExtrema();
        for (x_ext.values[0..x_ext.count]) |t| {
            const pt = self.evaluate(t);
            bbox.expand(pt.x, pt.y);
        }

        // Add y-extrema
        const y_ext = self.findYExtrema();
        for (y_ext.values[0..y_ext.count]) |t| {
            const pt = self.evaluate(t);
            bbox.expand(pt.x, pt.y);
        }

        return bbox;
    }

    /// Translate all control points
    pub fn translate(self: *CubicBezier, dx: f32, dy: f32) void {
        self.p0.x += dx;
        self.p0.y += dy;
        self.p1.x += dx;
        self.p1.y += dy;
        self.p2.x += dx;
        self.p2.y += dy;
        self.p3.x += dx;
        self.p3.y += dy;
    }
};

/// A Bezier segment that is x-monotonic (x always increasing or decreasing).
/// This allows efficient y-at-x queries.
pub const MonotonicSegment = struct {
    bezier: CubicBezier,
    x_coeffs: Coefficients,
    y_coeffs: Coefficients,
    x_min: f32,
    x_max: f32,
    y_min: f32,
    y_max: f32,
    x_increasing: bool,

    pub fn init(bezier: CubicBezier) MonotonicSegment {
        const bbox = bezier.boundingBox();
        return .{
            .bezier = bezier,
            .x_coeffs = bezier.xCoefficients(),
            .y_coeffs = bezier.yCoefficients(),
            .x_min = bbox.x_min,
            .x_max = bbox.x_max,
            .y_min = bbox.y_min,
            .y_max = bbox.y_max,
            .x_increasing = bezier.p3.x >= bezier.p0.x,
        };
    }

    /// Find y value at given x. Returns null if x is outside range.
    pub fn yAtX(self: MonotonicSegment, x: f32) ?f32 {
        const eps: f32 = 1e-9;
        if (x < self.x_min - eps or x > self.x_max + eps) {
            return null;
        }

        // Solve x(t) = x for t
        const t = cubic.solveCubicInUnit(
            self.x_coeffs.a,
            self.x_coeffs.b,
            self.x_coeffs.c,
            self.x_coeffs.d,
            x,
        );

        // Evaluate y(t)
        const t2 = t * t;
        const t3 = t2 * t;
        return self.y_coeffs.a * t3 + self.y_coeffs.b * t2 + self.y_coeffs.c * t + self.y_coeffs.d;
    }

    /// Get bounding box
    pub fn boundingBox(self: MonotonicSegment) BoundingBox {
        return BoundingBox.init(self.x_min, self.y_min, self.x_max, self.y_max);
    }
};

/// Split a Bezier curve into x-monotonic segments
pub fn splitToMonotonic(allocator: Allocator, bezier: CubicBezier) !std.ArrayList(MonotonicSegment) {
    var segments: std.ArrayList(MonotonicSegment) = .empty;
    errdefer segments.deinit(allocator);

    const extrema = bezier.findXExtrema();

    if (extrema.count == 0) {
        // Already monotonic
        try segments.append(allocator, MonotonicSegment.init(bezier));
        return segments;
    }

    // Sort extrema values
    var sorted_extrema: [2]f32 = extrema.values;
    if (extrema.count == 2 and sorted_extrema[0] > sorted_extrema[1]) {
        const tmp = sorted_extrema[0];
        sorted_extrema[0] = sorted_extrema[1];
        sorted_extrema[1] = tmp;
    }

    var current = bezier;
    var offset: f32 = 0.0;

    var i: u8 = 0;
    while (i < extrema.count) : (i += 1) {
        const t = sorted_extrema[i];
        // Adjust t for already-split portions
        const t_adj = (t - offset) / (1.0 - offset);
        if (t_adj <= 0 or t_adj >= 1) {
            continue;
        }

        const split = current.splitAt(t_adj);
        try segments.append(allocator, MonotonicSegment.init(split.left));
        current = split.right;
        offset = t;
    }

    // Add final segment
    try segments.append(allocator, MonotonicSegment.init(current));

    return segments;
}

// Tests
test "bezier evaluate" {
    const b = CubicBezier.init(
        Vec2.init(0, 0),
        Vec2.init(1, 2),
        Vec2.init(2, 2),
        Vec2.init(3, 0),
    );

    const p0 = b.evaluate(0);
    try std.testing.expectApproxEqAbs(@as(f32, 0), p0.x, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f32, 0), p0.y, 1e-6);

    const p1 = b.evaluate(1);
    try std.testing.expectApproxEqAbs(@as(f32, 3), p1.x, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f32, 0), p1.y, 1e-6);

    const mid = b.evaluate(0.5);
    try std.testing.expectApproxEqAbs(@as(f32, 1.5), mid.x, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f32, 1.5), mid.y, 1e-6);
}

test "bezier split" {
    const b = CubicBezier.init(
        Vec2.init(0, 0),
        Vec2.init(1, 2),
        Vec2.init(2, 2),
        Vec2.init(3, 0),
    );

    const split = b.splitAt(0.5);

    // Both halves should meet at the same point
    try std.testing.expectApproxEqAbs(split.left.p3.x, split.right.p0.x, 1e-6);
    try std.testing.expectApproxEqAbs(split.left.p3.y, split.right.p0.y, 1e-6);

    // Start and end should match original
    try std.testing.expectApproxEqAbs(b.p0.x, split.left.p0.x, 1e-6);
    try std.testing.expectApproxEqAbs(b.p3.x, split.right.p3.x, 1e-6);
}

test "monotonic y at x" {
    const b = CubicBezier.init(
        Vec2.init(0, 0),
        Vec2.init(1, 1),
        Vec2.init(2, 1),
        Vec2.init(3, 0),
    );

    const seg = MonotonicSegment.init(b);

    // At x=0, y should be 0
    const y0 = seg.yAtX(0);
    try std.testing.expect(y0 != null);
    try std.testing.expectApproxEqAbs(@as(f32, 0), y0.?, 1e-5);

    // At x=3, y should be 0
    const y3 = seg.yAtX(3);
    try std.testing.expect(y3 != null);
    try std.testing.expectApproxEqAbs(@as(f32, 0), y3.?, 1e-5);

    // At x=1.5 (middle), y should be around 0.75
    const ymid = seg.yAtX(1.5);
    try std.testing.expect(ymid != null);
    try std.testing.expect(ymid.? > 0.5 and ymid.? < 1.0);
}
