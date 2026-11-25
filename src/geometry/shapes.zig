//! ALGORITHMS FIRST: slab-decomposed segment-tree renderer.
//! Do not replace the algorithm with brute-force point checks.
//! Preserve algorithmic structure; micro-optimizations come second.

const std = @import("std");
const Allocator = std.mem.Allocator;
const bezier = @import("bezier.zig");
const stabbable = @import("../stabbing/stabbable.zig");

pub const Vec2 = bezier.Vec2;
pub const CubicBezier = bezier.CubicBezier;
pub const MonotonicSegment = bezier.MonotonicSegment;
pub const BoundingBox = bezier.BoundingBox;
pub const YInterval = stabbable.YInterval;
pub const Stabbable = stabbable.Stabbable;

/// RGBA color
pub const Color = struct {
    r: f32,
    g: f32,
    b: f32,
    a: f32,

    pub fn init(r: f32, g: f32, b: f32, a: f32) Color {
        return .{ .r = r, .g = g, .b = b, .a = a };
    }

    pub fn toArray(self: Color) [4]f32 {
        return .{ self.r, self.g, self.b, self.a };
    }

    /// Create color from HSL values
    pub fn fromHSL(h: f32, s: f32, l: f32, a: f32) Color {
        const c = (1.0 - @abs(2.0 * l - 1.0)) * s;
        const hue6 = h * 6.0;
        const x = c * (1.0 - @abs(@mod(hue6, 2.0) - 1.0));
        const m = l - c / 2.0;

        var r: f32 = 0;
        var g: f32 = 0;
        var b: f32 = 0;

        if (hue6 < 1.0) {
            r = c;
            g = x;
            b = 0;
        } else if (hue6 < 2.0) {
            r = x;
            g = c;
            b = 0;
        } else if (hue6 < 3.0) {
            r = 0;
            g = c;
            b = x;
        } else if (hue6 < 4.0) {
            r = 0;
            g = x;
            b = c;
        } else if (hue6 < 5.0) {
            r = x;
            g = 0;
            b = c;
        } else {
            r = c;
            g = 0;
            b = x;
        }

        return .{ .r = r + m, .g = g + m, .b = b + m, .a = a };
    }
};

/// A convex shape defined by cubic Bezier curves forming a closed boundary.
/// Curves should be in counterclockwise order.
pub const ConvexBezierShape = struct {
    allocator: Allocator,
    curves: std.ArrayList(CubicBezier),
    monotonic_segments: std.ArrayList(MonotonicSegment),
    color: Color,
    bbox: BoundingBox,
    offset: Vec2,
    z_index: i32,

    pub fn init(allocator: Allocator) ConvexBezierShape {
        return .{
            .allocator = allocator,
            .curves = .empty,
            .monotonic_segments = .empty,
            .color = Color.init(1, 1, 1, 1),
            .bbox = BoundingBox.empty(),
            .offset = Vec2.init(0, 0),
            .z_index = 0,
        };
    }

    pub fn deinit(self: *ConvexBezierShape) void {
        self.curves.deinit(self.allocator);
        self.monotonic_segments.deinit(self.allocator);
    }

    /// Add a curve to the shape
    pub fn addCurve(self: *ConvexBezierShape, curve: CubicBezier) !void {
        try self.curves.append(self.allocator, curve);
    }

    /// Set shape color
    pub fn setColor(self: *ConvexBezierShape, color: Color) void {
        self.color = color;
    }

    pub fn setZIndex(self: *ConvexBezierShape, z_index: i32) void {
        self.z_index = z_index;
    }

    /// Precompute monotonic segments and bounding box (call after adding all curves)
    pub fn precompute(self: *ConvexBezierShape) !void {
        self.monotonic_segments.clearRetainingCapacity();
        self.bbox = BoundingBox.empty();

        for (self.curves.items) |curve| {
            var segments = try bezier.splitToMonotonic(self.allocator, curve);
            defer segments.deinit(self.allocator);

            for (segments.items) |seg| {
                try self.monotonic_segments.append(self.allocator, seg);
                self.bbox.expandBox(seg.boundingBox());
            }
        }
    }

    /// Translate all control points
    pub fn translate(self: *ConvexBezierShape, dx: f32, dy: f32) void {
        self.translateFast(dx, dy);
    }

    /// Fast translate that keeps precomputed segments and bbox in sync.
    /// Safe when only translation occurs (no topology changes).
    pub fn translateFast(self: *ConvexBezierShape, dx: f32, dy: f32) void {
        if (dx == 0 and dy == 0) return;

        self.offset.x += dx;
        self.offset.y += dy;

        if (self.bbox.x_min <= self.bbox.x_max) {
            self.bbox.x_min += dx;
            self.bbox.x_max += dx;
            self.bbox.y_min += dy;
            self.bbox.y_max += dy;
        }
    }

    // Stabbable interface implementation

    pub fn boundingBox(self: *ConvexBezierShape) BoundingBox {
        return self.bbox;
    }

    /// Find the y-interval [y_min, y_max] of this shape at the given x
    pub fn yIntervalAtX(self: *ConvexBezierShape, x: f32) ?YInterval {
        const eps: f32 = 1e-9;
        if (x < self.bbox.x_min - eps or x > self.bbox.x_max + eps) {
            return null;
        }

        const local_x = x - self.offset.x;
        var y_min: f32 = std.math.floatMax(f32);
        var y_max: f32 = -std.math.floatMax(f32);
        var count: u32 = 0;

        for (self.monotonic_segments.items) |seg| {
            if (seg.yAtX(local_x)) |y| {
                const yw = y + self.offset.y;
                y_min = @min(y_min, yw);
                y_max = @max(y_max, yw);
                count += 1;
            }
        }

        if (count < 2) {
            return null;
        }

        return YInterval.init(y_min, y_max);
    }

    /// Check if point (x, y) is inside this convex shape
    pub fn containsPoint(self: *ConvexBezierShape, x: f32, y: f32) bool {
        const interval = self.yIntervalAtX(x);
        if (interval) |ival| {
            return ival.contains(y);
        }
        return false;
    }

    pub fn segmentCount(self: *ConvexBezierShape) usize {
        return self.monotonic_segments.items.len;
    }

    /// Pack segment data into buffer (12 floats per segment)
    pub fn packSegments(self: *ConvexBezierShape, buffer: []f32) void {
        for (self.monotonic_segments.items, 0..) |seg, i| {
            const base = i * 12;
            if (base + 12 > buffer.len) break;

            buffer[base + 0] = seg.x_coeffs.a;
            buffer[base + 1] = seg.x_coeffs.b;
            buffer[base + 2] = seg.x_coeffs.c;
            buffer[base + 3] = seg.x_coeffs.d;
            buffer[base + 4] = seg.y_coeffs.a;
            buffer[base + 5] = seg.y_coeffs.b;
            buffer[base + 6] = seg.y_coeffs.c;
            buffer[base + 7] = seg.y_coeffs.d;
            buffer[base + 8] = seg.x_min;
            buffer[base + 9] = seg.x_max;
            buffer[base + 10] = 0.5; // t_guess
            buffer[base + 11] = 0; // padding
        }
    }

    pub fn getColor(self: *ConvexBezierShape) [4]f32 {
        return self.color.toArray();
    }

    pub fn getZIndex(self: *ConvexBezierShape) i32 {
        return self.z_index;
    }

    pub fn getOffset(self: *ConvexBezierShape) Vec2 {
        return self.offset;
    }

    /// Create a Stabbable interface for this shape
    pub fn asStabbable(self: *ConvexBezierShape) Stabbable {
        return stabbable.makeStabbable(ConvexBezierShape, self);
    }
};

/// Create a random convex shape using cubic Bezier curves
pub fn createRandomConvexShape(
    allocator: Allocator,
    rng: std.Random,
    cx: f32,
    cy: f32,
    size: f32,
    num_sides: u32,
) !ConvexBezierShape {
    var shape = ConvexBezierShape.init(allocator);
    errdefer shape.deinit();

    // Generate convex hull points with evenly-ish spaced angles
    var points: std.ArrayList(Vec2) = .empty;
    defer points.deinit(allocator);

    const n = @max(3, num_sides);
    const angle_step = 2.0 * std.math.pi / @as(f32, @floatFromInt(n));

    var i: u32 = 0;
    while (i < n) : (i += 1) {
        const base_angle = @as(f32, @floatFromInt(i)) * angle_step;
        // Add small random perturbation
        const angle = base_angle + (rng.float(f32) - 0.5) * 0.6;
        // Random radius (not too varied to ensure convexity)
        const r = size * (0.7 + rng.float(f32) * 0.3);

        try points.append(allocator, Vec2.init(
            cx + r * @cos(angle),
            cy + r * @sin(angle),
        ));
    }

    // Create Bezier curves between consecutive points
    i = 0;
    while (i < n) : (i += 1) {
        const p0 = points.items[i];
        const p3 = points.items[(i + 1) % n];

        const chord = p3.sub(p0);
        const chord_len = chord.length();

        var p1: Vec2 = undefined;
        var p2: Vec2 = undefined;

        if (chord_len > 1e-6) {
            // Perpendicular direction (outward from center)
            const mid = p0.lerp(p3, 0.5);
            var outward = mid.sub(Vec2.init(cx, cy));
            const outward_len = outward.length();
            if (outward_len > 1e-6) {
                outward = outward.scale(1.0 / outward_len);
            } else {
                outward = Vec2.init(0, 0);
            }

            // Control points with slight outward bulge for roundness
            const bulge = chord_len * 0.15;
            p1 = p0.add(chord.scale(1.0 / 3.0)).add(outward.scale(bulge));
            p2 = p0.add(chord.scale(2.0 / 3.0)).add(outward.scale(bulge));
        } else {
            p1 = p0;
            p2 = p3;
        }

        try shape.addCurve(CubicBezier.init(p0, p1, p2, p3));
    }

    // Random pastel color
    const hue = rng.float(f32);
    const sat = 0.3 + rng.float(f32) * 0.3;
    const light = 0.7 + rng.float(f32) * 0.2;
    // Slightly more transparent to improve overlap readability.
    const alpha = 0.4 + rng.float(f32) * 0.3;

    shape.setColor(Color.fromHSL(hue, sat, light, alpha));

    try shape.precompute();

    return shape;
}

// Tests
test "convex shape creation" {
    const allocator = std.testing.allocator;

    var shape = ConvexBezierShape.init(allocator);
    defer shape.deinit();

    // Add a simple diamond shape
    try shape.addCurve(CubicBezier.init(
        Vec2.init(1, 0),
        Vec2.init(1, 0.5),
        Vec2.init(0.5, 1),
        Vec2.init(0, 1),
    ));
    try shape.addCurve(CubicBezier.init(
        Vec2.init(0, 1),
        Vec2.init(-0.5, 1),
        Vec2.init(-1, 0.5),
        Vec2.init(-1, 0),
    ));
    try shape.addCurve(CubicBezier.init(
        Vec2.init(-1, 0),
        Vec2.init(-1, -0.5),
        Vec2.init(-0.5, -1),
        Vec2.init(0, -1),
    ));
    try shape.addCurve(CubicBezier.init(
        Vec2.init(0, -1),
        Vec2.init(0.5, -1),
        Vec2.init(1, -0.5),
        Vec2.init(1, 0),
    ));

    try shape.precompute();

    // Check bounding box
    const bbox = shape.boundingBox();
    try std.testing.expect(bbox.x_min <= -1 and bbox.x_max >= 1);
    try std.testing.expect(bbox.y_min <= -1 and bbox.y_max >= 1);

    // Check containment
    try std.testing.expect(shape.containsPoint(0, 0)); // Center should be inside
}

test "random shape creation" {
    const allocator = std.testing.allocator;
    var prng = std.Random.DefaultPrng.init(12345);

    var shape = try createRandomConvexShape(allocator, prng.random(), 5.0, 5.0, 1.0, 6);
    defer shape.deinit();

    // Shape should have been created successfully
    try std.testing.expect(shape.curves.items.len == 6);
    try std.testing.expect(shape.monotonic_segments.items.len >= 6);
}
