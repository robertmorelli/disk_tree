//! ALGORITHMS FIRST: slab-decomposed segment-tree renderer.
//! Do not replace the algorithm with brute-force point checks.
//! Preserve algorithmic structure; micro-optimizations come second.

const std = @import("std");
const math = std.math;

/// Result of solving a cubic equation - up to 3 real roots
pub const CubicRoots = struct {
    roots: [3]f32,
    count: u8,

    pub fn init() CubicRoots {
        return .{ .roots = .{ 0, 0, 0 }, .count = 0 };
    }

    pub fn add(self: *CubicRoots, root: f32) void {
        if (self.count < 3) {
            self.roots[self.count] = root;
            self.count += 1;
        }
    }

    pub fn slice(self: *const CubicRoots) []const f32 {
        return self.roots[0..self.count];
    }
};

/// Solve cubic equation ax^3 + bx^2 + cx + d = 0 using Cardano's formula.
/// Returns all real roots.
pub fn solveCubic(a: f32, b: f32, c: f32, d: f32) CubicRoots {
    var result = CubicRoots.init();
    const eps: f32 = 1e-10;

    // Degenerate case: not actually cubic
    if (@abs(a) < eps) {
        // Quadratic: bx^2 + cx + d = 0
        if (@abs(b) < eps) {
            // Linear: cx + d = 0
            if (@abs(c) < eps) {
                return result; // No roots (or infinite)
            }
            result.add(-d / c);
            return result;
        }

        const disc = c * c - 4 * b * d;
        if (disc < 0) {
            return result; // No real roots
        } else if (disc < eps) {
            result.add(-c / (2 * b));
            return result;
        } else {
            const sqrt_disc = @sqrt(disc);
            result.add((-c + sqrt_disc) / (2 * b));
            result.add((-c - sqrt_disc) / (2 * b));
            return result;
        }
    }

    // Normalize: x^3 + px^2 + qx + r = 0
    const p = b / a;
    const q = c / a;
    const r = d / a;

    // Substitute x = t - p/3 to get depressed cubic: t^3 + At + B = 0
    const A = q - p * p / 3.0;
    const B = r - p * q / 3.0 + 2.0 * p * p * p / 27.0;

    // Discriminant
    const disc = B * B / 4.0 + A * A * A / 27.0;

    if (disc > eps) {
        // One real root
        const sqrt_disc = @sqrt(disc);
        const u_val = -B / 2.0 + sqrt_disc;
        const v_val = -B / 2.0 - sqrt_disc;

        const u = math.copysign(math.pow(f32, @abs(u_val), 1.0 / 3.0), u_val);
        const v = math.copysign(math.pow(f32, @abs(v_val), 1.0 / 3.0), v_val);

        result.add(u + v - p / 3.0);
    } else if (disc < -eps) {
        // Three real roots (casus irreducibilis)
        const sqrt_neg_A3_27 = @sqrt(-A * A * A / 27.0);
        const arg = std.math.clamp(-B / 2.0 / sqrt_neg_A3_27, -1.0, 1.0);
        const theta = math.acos(arg);
        const sqrt_A = 2.0 * @sqrt(-A / 3.0);

        var k: u8 = 0;
        while (k < 3) : (k += 1) {
            const t = sqrt_A * @cos((theta + 2.0 * math.pi * @as(f32, @floatFromInt(k))) / 3.0);
            result.add(t - p / 3.0);
        }
    } else {
        // Double or triple root
        if (@abs(A) < eps) {
            result.add(-p / 3.0);
        } else {
            const u = math.copysign(math.pow(f32, @abs(B / 2.0), 1.0 / 3.0), -B);
            result.add(2.0 * u - p / 3.0);
            result.add(-u - p / 3.0);
        }
    }

    return result;
}

/// Find the t in [0, 1] where a*t^3 + b*t^2 + c*t + d equals x_target, for a
/// cubic that is monotonic over [0, 1]. Returns the nearer end if the target is
/// outside the cubic's range there.
pub fn solveCubicInUnit(a: f32, b: f32, c: f32, d: f32, x_target: f32) f32 {
    const f0 = d - x_target;
    const f1 = a + b + c + d - x_target;
    if (f0 == 0) return 0;
    if (f1 == 0) return 1;
    if ((f0 < 0) == (f1 < 0)) return if (@abs(f0) < @abs(f1)) 0 else 1;

    // Newton-Raphson, bisecting whenever a step would leave the bracket.
    var lo: f32 = 0;
    var hi: f32 = 1;
    var t = f0 / (f0 - f1);
    var iter: u32 = 0;
    while (iter < 48) : (iter += 1) {
        const f = ((a * t + b) * t + c) * t + f0;
        if (f == 0) break;
        if ((f < 0) == (f0 < 0)) lo = t else hi = t;

        const fp = (3.0 * a * t + 2.0 * b) * t + c;
        var next = t - f / fp;
        if (!(next > lo and next < hi)) next = lo + (hi - lo) / 2.0;
        if (next == t) break;
        t = next;
    }

    return t;
}

// Tests
test "solve linear" {
    // 2x + 4 = 0 -> x = -2
    const roots = solveCubic(0, 0, 2, 4);
    try std.testing.expectEqual(@as(u8, 1), roots.count);
    try std.testing.expectApproxEqAbs(@as(f32, -2.0), roots.roots[0], 1e-6);
}

test "solve quadratic" {
    // x^2 - 5x + 6 = 0 -> x = 2, 3
    const roots = solveCubic(0, 1, -5, 6);
    try std.testing.expectEqual(@as(u8, 2), roots.count);

    var found2 = false;
    var found3 = false;
    for (roots.slice()) |r| {
        if (@abs(r - 2.0) < 1e-6) found2 = true;
        if (@abs(r - 3.0) < 1e-6) found3 = true;
    }
    try std.testing.expect(found2 and found3);
}

test "solve cubic one root" {
    // x^3 - 1 = 0 -> x = 1
    const roots = solveCubic(1, 0, 0, -1);
    try std.testing.expect(roots.count >= 1);

    var found1 = false;
    for (roots.slice()) |r| {
        if (@abs(r - 1.0) < 1e-5) found1 = true;
    }
    try std.testing.expect(found1);
}

test "solve cubic three roots" {
    // x^3 - 6x^2 + 11x - 6 = 0 -> x = 1, 2, 3
    const roots = solveCubic(1, -6, 11, -6);
    try std.testing.expectEqual(@as(u8, 3), roots.count);

    var found1 = false;
    var found2 = false;
    var found3 = false;
    for (roots.slice()) |r| {
        if (@abs(r - 1.0) < 1e-5) found1 = true;
        if (@abs(r - 2.0) < 1e-5) found2 = true;
        if (@abs(r - 3.0) < 1e-5) found3 = true;
    }
    try std.testing.expect(found1 and found2 and found3);
}
