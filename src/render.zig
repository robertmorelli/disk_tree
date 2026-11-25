//! Renders translucent shapes by asking a tree for the blended colour under
//! each pixel. The cost per pixel does not depend on how many shapes overlap.
//!
//!     render disks out.ppm
//!     render paths out.ppm

const std = @import("std");
const disk_tree = @import("disk_tree.zig");
const path_tree = @import("path_tree.zig");
const shapes = @import("geometry/shapes.zig");

/// Premultiplied RGBA.
const Colour = [4]f32;

/// `front` composited over `back`.
fn over(back: Colour, front: Colour) Colour {
    const keep: @Vector(4, f32) = @splat(1 - front[3]);
    return @as(@Vector(4, f32), front) + keep * @as(@Vector(4, f32), back);
}

const DiskTree = disk_tree.Folded(Colour, over);
const PathTree = path_tree.FoldedPaths(Colour, over);

const size = 1024;

fn randomColour(random: std.Random) Colour {
    const alpha = 0.15 + 0.35 * random.float(f32);
    return .{ alpha * random.float(f32), alpha * random.float(f32), alpha * random.float(f32), alpha };
}

fn diskScene(gpa: std.mem.Allocator, random: std.Random) !DiskTree {
    var circles: [600]disk_tree.Circle = undefined;
    var colours: [circles.len]Colour = undefined;
    for (&circles, &colours) |*c, *colour| {
        c.* = .{
            .cx = random.float(f64) * size,
            .cy = random.float(f64) * size,
            .r = 20 + random.float(f64) * 140,
        };
        colour.* = randomColour(random);
    }
    return DiskTree.init(gpa, &circles, &colours);
}

/// Rounded blobs and five-pointed stars. A star is drawn as one contour that
/// crosses itself, so the even-odd rule leaves its centre empty.
fn pathScene(gpa: std.mem.Allocator, random: std.Random) !PathTree {
    var curves: std.ArrayList(path_tree.CubicBezier) = .empty;
    defer curves.deinit(gpa);
    var lens: [400]usize = undefined;
    var colours: [lens.len]Colour = undefined;

    for (&lens, &colours, 0..) |*len, *colour, i| {
        const cx = random.float(f32) * size;
        const cy = random.float(f32) * size;
        const radius = 30 + random.float(f32) * 130;
        const before = curves.items.len;
        if (i % 4 == 0) {
            const turn = random.float(f32) * std.math.tau;
            var points: [5]path_tree.Vec2 = undefined;
            for (&points, 0..) |*p, k| {
                const angle = turn + std.math.tau * 2.0 * @as(f32, @floatFromInt(k)) / 5.0;
                p.* = .init(cx + radius * @cos(angle), cy + radius * @sin(angle));
            }
            for (points, 0..) |a, k| {
                const b = points[(k + 1) % points.len];
                try curves.append(gpa, .init(a, a.lerp(b, 1.0 / 3.0), a.lerp(b, 2.0 / 3.0), b));
            }
        } else {
            var blob = try shapes.createRandomConvexShape(gpa, random, cx, cy, radius, 3 + random.uintLessThan(u32, 6));
            defer blob.deinit();
            try curves.appendSlice(gpa, blob.curves.items);
        }
        len.* = curves.items.len - before;
        colour.* = randomColour(random);
    }

    var paths: [lens.len]path_tree.Path = undefined;
    var begin: usize = 0;
    for (&paths, lens) |*path, len| {
        path.* = .{ .curves = curves.items[begin..][0..len] };
        begin += len;
    }
    return PathTree.init(gpa, &paths, &colours);
}

fn shade(tree: anytype, pixels: [][3]u8) void {
    for (pixels, 0..) |*pixel, i| {
        const x: f64 = @floatFromInt(i % size);
        const y: f64 = @floatFromInt(i / size);
        const colour = tree.fold(x + 0.5, y + 0.5) orelse Colour{ 0, 0, 0, 0 };
        // Over a white background.
        for (pixel, colour[0..3]) |*channel, c| {
            channel.* = @intFromFloat(@round(255 * std.math.clamp(c + 1 - colour[3], 0, 1)));
        }
    }
}

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;
    const io = init.io;

    var args = init.minimal.args.iterate();
    _ = args.skip();
    const scene = args.next() orelse "disks";
    const path = args.next() orelse "out.ppm";

    var prng: std.Random.DefaultPrng = .init(11);
    const random = prng.random();
    const pixels = try gpa.alloc([3]u8, size * size);
    defer gpa.free(pixels);

    if (std.mem.eql(u8, scene, "paths")) {
        var tree = try pathScene(gpa, random);
        defer tree.deinit(gpa);
        shade(&tree, pixels);
    } else {
        var tree = try diskScene(gpa, random);
        defer tree.deinit(gpa);
        shade(&tree, pixels);
    }

    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);
    var buffer: [1 << 16]u8 = undefined;
    var writer = file.writer(io, &buffer);
    try writer.interface.print("P6\n{d} {d}\n255\n", .{ size, size });
    try writer.interface.writeAll(std.mem.sliceAsBytes(pixels));
    try writer.interface.flush();
}
