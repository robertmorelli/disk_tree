//! A crystal castle: thousands of translucent Bézier shapes stacked deep, drawn
//! with one fold per sample. Writes a PPM image to the path given as the first
//! argument. The second argument is how many particles to float over it all,
//! and the third their largest radius in pixels.

const std = @import("std");
const path_tree = @import("path_tree.zig");

const V = path_tree.Vec2;
const Cubic = path_tree.CubicBezier;

const width = 1400;
const height = 1000;
/// Samples per pixel along each axis.
const supersample = 2;

const Rgba = @Vector(4, f32);

/// A premultiplied colour, and how many shapes went into it.
const Value = struct {
    rgba: Rgba,
    layers: u32 = 1,
};

/// `front` composited over `back`.
fn over(back: Value, front: Value) Value {
    const keep: Rgba = @splat(1 - front.rgba[3]);
    return .{ .rgba = front.rgba + keep * back.rgba, .layers = back.layers + front.layers };
}

const Tree = path_tree.FoldedPaths(Value, over);

/// Hue in degrees; the result is premultiplied by `alpha`.
fn tint(hue: f32, saturation: f32, lightness: f32, alpha: f32) Value {
    const h = @mod(hue, 360) / 60;
    const chroma = (1 - @abs(2 * lightness - 1)) * saturation;
    const x = chroma * (1 - @abs(@mod(h, 2) - 1));
    const rgb: [3]f32 = switch (@as(u32, @intFromFloat(h))) {
        0 => .{ chroma, x, 0 },
        1 => .{ x, chroma, 0 },
        2 => .{ 0, chroma, x },
        3 => .{ 0, x, chroma },
        4 => .{ x, 0, chroma },
        else => .{ chroma, 0, x },
    };
    const m = lightness - chroma / 2;
    return .{ .rgba = .{ alpha * (rgb[0] + m), alpha * (rgb[1] + m), alpha * (rgb[2] + m), alpha } };
}

const Scene = struct {
    gpa: std.mem.Allocator,
    random: std.Random,
    curves: std.ArrayList(Cubic) = .empty,
    /// Start and length in `curves` of each finished shape, back to front.
    shapes: std.ArrayList([2]usize) = .empty,
    values: std.ArrayList(Value) = .empty,
    /// Where the shape being drawn starts.
    open: usize = 0,

    /// A point, nudged so that no two edges are exactly aligned: hand-cut
    /// crystal, and no degenerate geometry.
    fn at(self: *Scene, x: f32, y: f32) V {
        return V.init(x + (self.random.float(f32) - 0.5) * 0.4, y + (self.random.float(f32) - 0.5) * 0.4);
    }

    fn range(self: *Scene, lo: f32, hi: f32) f32 {
        return lo + self.random.float(f32) * (hi - lo);
    }

    fn line(self: *Scene, a: V, b: V) !void {
        try self.curves.append(self.gpa, Cubic.init(a, a.lerp(b, 1.0 / 3.0), a.lerp(b, 2.0 / 3.0), b));
    }

    fn curve(self: *Scene, p0: V, p1: V, p2: V, p3: V) !void {
        try self.curves.append(self.gpa, Cubic.init(p0, p1, p2, p3));
    }

    /// Finishes the shape drawn since the last fill.
    fn fill(self: *Scene, value: Value) !void {
        try self.shapes.append(self.gpa, .{ self.open, self.curves.items.len - self.open });
        try self.values.append(self.gpa, value);
        self.open = self.curves.items.len;
    }

    fn outline(self: *Scene, points: []const V) !void {
        for (points, 0..) |p, i| try self.line(p, points[(i + 1) % points.len]);
    }

    fn polygon(self: *Scene, points: []const V, value: Value) !void {
        try self.outline(points);
        try self.fill(value);
    }

    /// A circular contour. Two of them in one shape make a ring.
    fn circle(self: *Scene, c: V, r: f32) !void {
        const k = 0.5522848 * r;
        const e = V.init(c.x + r, c.y);
        const s = V.init(c.x, c.y + r);
        const w = V.init(c.x - r, c.y);
        const n = V.init(c.x, c.y - r);
        try self.curve(e, V.init(e.x, e.y + k), V.init(s.x + k, s.y), s);
        try self.curve(s, V.init(s.x - k, s.y), V.init(w.x, w.y + k), w);
        try self.curve(w, V.init(w.x, w.y - k), V.init(n.x - k, n.y), n);
        try self.curve(n, V.init(n.x + k, n.y), V.init(e.x, e.y - k), e);
    }

    fn ring(self: *Scene, c: V, r: f32, thickness: f32, value: Value) !void {
        try self.circle(c, r);
        try self.circle(V.init(c.x + 0.01, c.y), r - thickness);
        try self.fill(value);
    }

    fn bead(self: *Scene, c: V, r: f32, value: Value) !void {
        try self.circle(c, r);
        try self.fill(value);
    }

    /// A drop hanging from `top`.
    fn teardrop(self: *Scene, top: V, w: f32, h: f32, value: Value) !void {
        const bottom = V.init(top.x, top.y + h);
        try self.curve(top, V.init(top.x + w, top.y + h * 0.45), V.init(top.x + w, bottom.y), bottom);
        try self.curve(bottom, V.init(top.x - w, bottom.y), V.init(top.x - w, top.y + h * 0.45), top);
        try self.fill(value);
    }

    /// A cut gem: a diamond with a darker half and a bright table.
    fn gem(self: *Scene, c: V, w: f32, h: f32, hue: f32, alpha: f32) !void {
        const top = self.at(c.x, c.y - h);
        const bottom = self.at(c.x, c.y + h);
        const left = self.at(c.x - w, c.y);
        const right = self.at(c.x + w, c.y);
        try self.polygon(&.{ top, right, bottom, left }, tint(hue, 0.7, 0.6, alpha));
        try self.polygon(&.{ top, right, bottom }, tint(hue + 20, 0.8, 0.3, alpha * 0.7));
        try self.polygon(&.{ top, c.lerp(right, 0.45), c.lerp(left, 0.45) }, tint(hue - 20, 0.5, 0.9, alpha * 0.8));
    }

    /// A tapering band along a cubic.
    fn ribbon(self: *Scene, p0: V, p1: V, p2: V, p3: V, w: f32, value: Value) !void {
        const up = V.init(0, -w);
        const down = V.init(0, w);
        try self.curve(p0.add(up), p1.add(up), p2.add(up.scale(0.5)), p3);
        try self.curve(p3, p2.add(down.scale(0.5)), p1.add(down), p0.add(down));
        try self.fill(value);
    }

    /// A four-pointed glint.
    fn sparkle(self: *Scene, c: V, r: f32, alpha: f32) !void {
        const t = r * 0.12;
        try self.polygon(&.{
            self.at(c.x, c.y - r), self.at(c.x + t, c.y - t),
            self.at(c.x + r, c.y), self.at(c.x + t, c.y + t),
            self.at(c.x, c.y + r), self.at(c.x - t, c.y + t),
            self.at(c.x - r, c.y), self.at(c.x - t, c.y - t),
        }, tint(320, 0.6, 0.95, alpha));
    }

    /// A pointed crystal growing from `base`. `angle` is from straight up.
    fn shard(self: *Scene, base: V, angle: f32, len: f32, w: f32, hue: f32, alpha: f32) !void {
        const dir = V.init(@sin(angle), -@cos(angle));
        const side = V.init(@cos(angle), @sin(angle));
        const waist = base.add(dir.scale(len * 0.3));
        const tip = self.at(base.x + dir.x * len, base.y + dir.y * len);
        const left = waist.sub(side.scale(w));
        const right = waist.add(side.scale(w));
        const foot = self.at(base.x, base.y);
        try self.polygon(&.{ foot, left, tip }, tint(hue, 0.65, 0.58, alpha));
        try self.polygon(&.{ foot, right, tip }, tint(hue + 18, 0.75, 0.32, alpha));
        try self.polygon(&.{ foot, waist.sub(side.scale(w * 0.3)), tip, waist.add(side.scale(w * 0.1)) }, tint(hue - 15, 0.5, 0.88, alpha * 0.55));
    }

    fn cluster(self: *Scene, x: f32, y: f32, count: usize, len: f32, hue: f32) !void {
        for (0..count) |_| {
            const angle = self.range(-1.25, 1.25);
            const l = len * self.range(0.35, 1.0) * (1.15 - 0.5 * @abs(angle));
            try self.shard(V.init(x + self.range(-len, len) * 0.25, y), angle, l, l * self.range(0.07, 0.16), hue + self.range(-25, 25), self.range(0.28, 0.5));
        }
    }

    /// A pointed window, lit from inside.
    fn arch(self: *Scene, x: f32, bottom: f32, half: f32, h: f32, value: Value) !void {
        const a = self.at(x - half, bottom);
        const b = self.at(x - half, bottom - h * 0.5);
        const top = self.at(x, bottom - h);
        const c = self.at(x + half, bottom - h * 0.5);
        const d = self.at(x + half, bottom);
        try self.line(a, b);
        try self.curve(b, V.init(b.x, bottom - h * 0.8), V.init(x - half * 0.4, top.y), top);
        try self.curve(top, V.init(x + half * 0.4, top.y), V.init(c.x, bottom - h * 0.8), c);
        try self.line(c, d);
        try self.fill(value);
    }

    /// A faceted tower with a spire, bands and windows.
    fn tower(self: *Scene, x: f32, base: f32, w: f32, h: f32, hue: f32, alpha: f32) !void {
        const facets = 5;
        const lean = self.range(-0.02, 0.02) * h;
        const apex = self.at(x + lean * 1.6, base - h - w * self.range(1.3, 2.2));
        for (0..facets) |i| {
            const f: f32 = @floatFromInt(i);
            const x0 = x - w / 2 + w * f / facets - w * 0.05;
            const x1 = x - w / 2 + w * (f + 1) / facets + w * 0.05;
            const light = 0.22 + 0.42 * @abs(@sin(f * 1.9 + hue));
            const top0 = self.at(x0 + lean, base - h);
            const top1 = self.at(x1 + lean, base - h);
            try self.polygon(&.{ self.at(x0, base), self.at(x1, base), top1, top0 }, tint(hue + f * 7, 0.65, light, alpha));
            try self.polygon(&.{ top0, top1, apex }, tint(hue + f * 7 + 12, 0.7, light + 0.12, alpha * 1.1));
        }

        const bands: usize = @intFromFloat(@max(2, h / 55));
        for (1..bands) |i| {
            const y = base - h * @as(f32, @floatFromInt(i)) / @as(f32, @floatFromInt(bands));
            const shift = lean * @as(f32, @floatFromInt(i)) / @as(f32, @floatFromInt(bands));
            try self.polygon(&.{
                self.at(x - w * 0.58 + shift, y),     self.at(x + w * 0.58 + shift, y),
                self.at(x + w * 0.58 + shift, y - 5), self.at(x - w * 0.58 + shift, y - 5),
            }, tint(hue + 35, 0.5, 0.82, 0.22));
            if (w > 26) {
                try self.arch(x + shift, y + h / @as(f32, @floatFromInt(bands)) * 0.75, w * 0.16, w * 0.62, tint(335, 0.9, 0.86, 0.5));
                try self.arch(x + shift, y + h / @as(f32, @floatFromInt(bands)) * 0.72, w * 0.09, w * 0.42, tint(40, 0.9, 0.93, 0.45));
            }
        }
        try self.sparkle(apex, w * 0.5, 0.75);
    }

    /// An arched crystal bridge between two points.
    fn bridge(self: *Scene, a: V, b: V, rise: f32, hue: f32) !void {
        const p1 = a.lerp(b, 0.3).add(V.init(0, -rise));
        const p2 = a.lerp(b, 0.7).add(V.init(0, -rise));
        try self.ribbon(a, p1, p2, b, 9, tint(hue, 0.6, 0.6, 0.35));
        try self.ribbon(a.add(V.init(0, 14)), p1.add(V.init(0, 30)), p2.add(V.init(0, 30)), b.add(V.init(0, 14)), 4, tint(hue + 30, 0.6, 0.85, 0.3));
        const span = Cubic.init(a, p1, p2, b);
        for (1..12) |i| {
            const p = span.evaluate(@as(f32, @floatFromInt(i)) / 12);
            try self.gem(V.init(p.x, p.y - 14), 4, 9, hue + 40, 0.5);
        }
    }

    /// A string of gems hanging from `top`, ending in a drop.
    fn strand(self: *Scene, top: V, count: usize, size: f32, hue: f32) !void {
        var y = top.y;
        for (0..count) |i| {
            const s = size * (1 - 0.08 * @as(f32, @floatFromInt(i)));
            y += s * 1.5;
            try self.gem(V.init(top.x, y), s * 0.55, s, hue + @as(f32, @floatFromInt(i)) * 9, 0.55);
            y += s * 1.2;
        }
        try self.teardrop(V.init(top.x, y), size * 0.9, size * 2.6, tint(hue + 25, 0.8, 0.7, 0.6));
        try self.teardrop(V.init(top.x - size * 0.2, y + size * 0.6), size * 0.3, size * 1.1, tint(hue, 0.4, 0.95, 0.6));
    }

    fn chandelier(self: *Scene, cx: f32, top: f32, size: f32, hue: f32) !void {
        // Chain and crown.
        try self.polygon(&.{ self.at(cx - 1.5, -5), self.at(cx + 1.5, -5), self.at(cx + 1.5, top), self.at(cx - 1.5, top) }, tint(45, 0.6, 0.75, 0.7));
        var link: f32 = 8;
        while (link < top) : (link += 14) try self.ring(V.init(cx, link), 5, 1.4, tint(45, 0.7, 0.8, 0.7));
        try self.ring(V.init(cx, top), size * 0.16, size * 0.03, tint(45, 0.8, 0.8, 0.8));
        try self.gem(V.init(cx, top), size * 0.07, size * 0.12, hue, 0.7);

        const tiers = 3;
        for (0..tiers) |tier| {
            const k: f32 = @floatFromInt(tier);
            const y = top + size * (0.25 + 0.5 * k);
            const reach = size * (0.55 + 0.4 * k);
            const arms: usize = 4 + 2 * tier;

            try self.ring(V.init(cx, y), size * 0.1, size * 0.025, tint(45, 0.8, 0.78, 0.8));
            var tips: [2][8]V = undefined;
            for (0..2) |side| {
                const sign: f32 = if (side == 0) -1 else 1;
                for (0..arms) |arm| {
                    const f = (@as(f32, @floatFromInt(arm)) + 1) / @as(f32, @floatFromInt(arms));
                    const hub = V.init(cx, y);
                    const tip = self.at(cx + sign * reach * f, y - reach * 0.16);
                    // An S-curved arm, a finer scroll under it, and a curl at the hub.
                    try self.ribbon(hub, V.init(cx + sign * reach * f * 0.45, y + reach * 0.42), V.init(tip.x, y + reach * 0.3), tip, size * 0.02, tint(48, 0.75, 0.7, 0.75));
                    try self.ribbon(hub, V.init(cx + sign * reach * f * 0.2, y - reach * 0.3), V.init(cx + sign * reach * f * 0.7, y - reach * 0.34), V.init(cx + sign * reach * f * 0.62, y + reach * 0.12), size * 0.012, tint(hue + 20, 0.6, 0.8, 0.5));

                    // Cup, candle and flame.
                    try self.ring(tip, size * 0.035, size * 0.012, tint(45, 0.8, 0.8, 0.8));
                    try self.polygon(&.{ self.at(tip.x - size * 0.012, tip.y), self.at(tip.x + size * 0.012, tip.y), self.at(tip.x + size * 0.012, tip.y - size * 0.11), self.at(tip.x - size * 0.012, tip.y - size * 0.11) }, tint(hue, 0.3, 0.92, 0.7));
                    try self.teardrop(V.init(tip.x, tip.y - size * 0.2), size * 0.03, size * 0.09, tint(35, 1.0, 0.8, 0.85));
                    try self.bead(V.init(tip.x, tip.y - size * 0.15), size * 0.11, tint(330, 1.0, 0.85, 0.1));

                    try self.strand(V.init(tip.x, tip.y + size * 0.03), 2 + tier + arm % 3, size * 0.035, hue + 10 * @as(f32, @floatFromInt(arm)));
                    tips[side][arm] = tip;
                }
                // Swags of beads between neighbouring arms.
                for (0..arms - 1) |arm| {
                    const a = tips[side][arm];
                    const b = tips[side][arm + 1];
                    const sag = V.init(0, size * 0.22);
                    const swag = Cubic.init(a, a.lerp(b, 0.25).add(sag), a.lerp(b, 0.75).add(sag), b);
                    try self.ribbon(swag.p0, swag.p1, swag.p2, swag.p3, 0.8, tint(hue, 0.4, 0.9, 0.5));
                    for (1..9) |i| {
                        const p = swag.evaluate(@as(f32, @floatFromInt(i)) / 9);
                        try self.bead(p, size * 0.014, tint(hue + 30, 0.8, 0.75, 0.6));
                    }
                }
            }
        }

        // The finial.
        const bottom = top + size * 1.55;
        try self.strand(V.init(cx, bottom - size * 0.3), 4, size * 0.06, hue);
        try self.ring(V.init(cx, bottom - size * 0.3), size * 0.2, size * 0.02, tint(45, 0.8, 0.8, 0.7));
        try self.sparkle(V.init(cx, bottom + size * 0.4), size * 0.25, 0.8);
    }

    /// Motes of light drifting in front of everything.
    fn particles(self: *Scene, count: usize, max_radius: f32) !void {
        for (0..count) |_| {
            // Mostly tiny, a few large.
            const size = self.random.float(f32);
            const r = 0.35 + (max_radius - 0.35) * size * size * size;
            const hue: f32 = switch (self.random.uintLessThan(u32, 4)) {
                0 => 45,
                1 => 330,
                2 => 285,
                else => 190,
            };
            try self.bead(V.init(self.range(0, width), self.range(0, height)), r, tint(hue, 0.9, 0.85, self.range(0.15, 0.6)));
        }
    }

    fn build(self: *Scene) !void {
        const ground: f32 = 905;
        const centre: f32 = width / 2;

        // Far spires, huge and faint.
        for (0..70) |_| {
            try self.shard(V.init(self.range(-50, width + 50), ground + 60), self.range(-0.35, 0.35), self.range(350, 900), self.range(18, 60), self.range(255, 300), self.range(0.06, 0.14));
        }
        // Moons.
        try self.bead(V.init(1010, 330), 95, tint(310, 0.7, 0.85, 0.3));
        try self.ring(V.init(1010, 330), 125, 6, tint(290, 0.7, 0.8, 0.25));
        try self.bead(V.init(330, 300), 40, tint(280, 0.6, 0.9, 0.3));

        // The back row of the castle, then walls, then the front row.
        const back = [_][4]f32{
            .{ -430, 46, 360, 285 }, .{ -330, 38, 440, 300 }, .{ -210, 54, 520, 275 }, .{ -95, 44, 620, 290 },
            .{ 95, 44, 620, 290 },   .{ 210, 54, 520, 275 },  .{ 330, 38, 440, 300 },  .{ 430, 46, 360, 285 },
            .{ -540, 34, 250, 265 }, .{ 540, 34, 250, 265 },
        };
        for (back) |t| try self.tower(centre + t[0], ground, t[1], t[2], t[3], 0.26);

        for (0..6) |i| {
            const f: f32 = @floatFromInt(i);
            const half = 520 - f * 62;
            const top = ground - 150 - f * 48;
            try self.polygon(&.{ self.at(centre - half, ground), self.at(centre + half, ground), self.at(centre + half * 0.93, top), self.at(centre - half * 0.93, top) }, tint(275 + f * 9, 0.6, 0.28 + f * 0.04, 0.3));
            // Battlements along each wall.
            var x = centre - half * 0.93;
            while (x < centre + half * 0.93 - 14) : (x += 28) {
                try self.polygon(&.{ self.at(x, top), self.at(x + 14, top), self.at(x + 14, top - 16), self.at(x, top - 16) }, tint(285 + f * 9, 0.6, 0.5, 0.35));
            }
        }

        try self.bridge(V.init(centre - 330, 520), V.init(centre - 95, 420), 90, 300);
        try self.bridge(V.init(centre + 95, 420), V.init(centre + 330, 520), 90, 300);
        try self.bridge(V.init(centre - 480, 640), V.init(centre - 270, 560), 70, 320);
        try self.bridge(V.init(centre + 270, 560), V.init(centre + 480, 640), 70, 320);

        const front = [_][4]f32{
            .{ -480, 62, 300, 320 }, .{ -270, 70, 420, 305 }, .{ 0, 96, 640, 290 },
            .{ 270, 70, 420, 305 },  .{ 480, 62, 300, 320 },  .{ -150, 40, 330, 330 },
            .{ 150, 40, 330, 330 },  .{ -380, 30, 210, 340 }, .{ 380, 30, 210, 340 },
        };
        for (front) |t| try self.tower(centre + t[0], ground, t[1], t[2], t[3], 0.34);

        // The great gate and its rose window.
        try self.arch(centre, ground, 62, 210, tint(330, 0.9, 0.8, 0.5));
        try self.arch(centre, ground, 44, 170, tint(30, 1.0, 0.9, 0.5));
        for (0..12) |i| {
            const a = std.math.tau * @as(f32, @floatFromInt(i)) / 12;
            const c = V.init(centre, 470);
            const p = V.init(c.x + 46 * @cos(a), c.y + 46 * @sin(a));
            try self.teardrop(V.init(p.x, p.y - 12), 9, 24, tint(300 + 12 * @as(f32, @floatFromInt(i % 5)), 0.9, 0.7, 0.5));
            try self.ring(c, 14 + 5 * @as(f32, @floatFromInt(i % 6)), 2, tint(45, 0.8, 0.8, 0.5));
        }
        try self.sparkle(V.init(centre, 470), 70, 0.7);

        // Amethyst growing everywhere along the ground and up the walls.
        var x: f32 = -20;
        while (x < width + 20) : (x += self.range(22, 46)) {
            try self.cluster(x, ground + self.range(0, 70), 9, self.range(70, 230), self.range(265, 325));
        }
        for (0..40) |_| {
            try self.cluster(centre + self.range(-560, 560), ground - self.range(60, 420), 6, self.range(30, 80), self.range(270, 330));
        }

        // Chandeliers, hung in front of it all.
        try self.chandelier(centre, 60, 135, 300);
        try self.chandelier(centre - 470, 55, 88, 285);
        try self.chandelier(centre + 470, 55, 88, 320);
        try self.chandelier(centre - 245, 30, 52, 330);
        try self.chandelier(centre + 245, 30, 52, 275);

        // Panes of glass over the top, and light catching everything.
        for (0..9) |i| {
            const f: f32 = @floatFromInt(i);
            const left = -200 + f * 190;
            try self.polygon(&.{
                self.at(left, -20),                     self.at(left + 260, -20),
                self.at(left + 260 - 380, height + 20), self.at(left - 380, height + 20),
            }, tint(190 + f * 14, 0.7, 0.85, 0.07));
        }
        for (0..150) |_| {
            try self.sparkle(V.init(self.range(0, width), self.range(0, height)), self.range(4, 22), self.range(0.4, 0.9));
        }
    }
};

/// The sky behind everything: dusk, with a glow behind the castle.
fn sky(x: f32, y: f32) [3]f32 {
    const t = y / height;
    const dx = (x - width / 2.0) / 620;
    const dy = (y - 560) / 420;
    const glow = 0.55 * @exp(-(dx * dx + dy * dy));
    return .{
        0.05 + 0.26 * t + glow * 0.9,
        0.02 + 0.05 * t + glow * 0.3,
        0.13 + 0.2 * t + glow * 0.75,
    };
}

fn seconds(from: std.Io.Timestamp, io: std.Io) f64 {
    return @as(f64, @floatFromInt(from.untilNow(io, .awake).nanoseconds)) / std.time.ns_per_s;
}

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;
    const io = init.io;

    var args = init.minimal.args.iterate();
    _ = args.skip();
    const path = args.next() orelse "castle.ppm";
    const particle_count = if (args.next()) |arg| try std.fmt.parseInt(usize, arg, 10) else 0;
    const particle_radius = if (args.next()) |arg| try std.fmt.parseFloat(f32, arg) else 2.5;

    var prng: std.Random.DefaultPrng = .init(7);
    var scene: Scene = .{ .gpa = gpa, .random = prng.random() };
    try scene.build();
    try scene.particles(particle_count, particle_radius);

    const paths = try gpa.alloc(path_tree.Path, scene.shapes.items.len);
    defer gpa.free(paths);
    for (paths, scene.shapes.items) |*p, shape| p.* = .{ .curves = scene.curves.items[shape[0]..][0..shape[1]] };

    var start = std.Io.Timestamp.now(io, .awake);
    var tree = try Tree.init(gpa, paths, scene.values.items);
    defer tree.deinit(gpa);
    const build_s = seconds(start, io);

    const pixels = try gpa.alloc([3]u8, width * height);
    defer gpa.free(pixels);
    var deepest: u32 = 0;
    var layers: u64 = 0;
    start = std.Io.Timestamp.now(io, .awake);
    for (pixels, 0..) |*pixel, i| {
        var sum: @Vector(3, f32) = @splat(0);
        for (0..supersample * supersample) |s| {
            const x = @as(f32, @floatFromInt(i % width)) + (@as(f32, @floatFromInt(s % supersample)) + 0.5) / supersample;
            const y = @as(f32, @floatFromInt(i / width)) + (@as(f32, @floatFromInt(s / supersample)) + 0.5) / supersample;
            const behind: @Vector(3, f32) = sky(x, y);
            if (tree.fold(x, y)) |v| {
                deepest = @max(deepest, v.layers);
                layers += v.layers;
                const keep: @Vector(3, f32) = @splat(1 - v.rgba[3]);
                sum += @Vector(3, f32){ v.rgba[0], v.rgba[1], v.rgba[2] } + keep * behind;
            } else {
                sum += behind;
            }
        }
        const colour: [3]f32 = sum / @as(@Vector(3, f32), @splat(supersample * supersample));
        for (pixel, colour) |*channel, c| channel.* = @intFromFloat(@round(255 * std.math.clamp(c, 0, 1)));
    }
    const render_s = seconds(start, io);
    const samples = width * height * supersample * supersample;

    std.debug.print(
        \\{d} shapes, {d} curves
        \\tree: {d} slabs, {d:.0} MiB, built in {d:.2} s
        \\render: {d} samples in {d:.2} s, {d:.2} us each
        \\layers under a sample: {d:.1} on average, {d} at most
        \\
    , .{
        paths.len,                                                          scene.curves.items.len,
        tree.xs.len,                                                        @as(f64, @floatFromInt(tree.bytes())) / (1 << 20),
        build_s,                                                            samples,
        render_s,                                                           render_s * 1e6 / @as(f64, @floatFromInt(samples)),
        @as(f64, @floatFromInt(layers)) / @as(f64, @floatFromInt(samples)), deepest,
    });

    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);
    var buffer: [1 << 16]u8 = undefined;
    var writer = file.writer(io, &buffer);
    try writer.interface.print("P6\n{d} {d}\n255\n", .{ width, height });
    try writer.interface.writeAll(std.mem.sliceAsBytes(pixels));
    try writer.interface.flush();
}
