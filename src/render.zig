//! Renders translucent disks by asking the tree for the blended colour under
//! each pixel. The cost per pixel does not depend on how many disks overlap.
//! Writes a PPM image to the path given as the first argument.

const std = @import("std");
const disk_tree = @import("disk_tree.zig");

/// Premultiplied RGBA.
const Colour = [4]f32;

/// `front` composited over `back`.
fn over(back: Colour, front: Colour) Colour {
    const keep: @Vector(4, f32) = @splat(1 - front[3]);
    return @as(@Vector(4, f32), front) + keep * @as(@Vector(4, f32), back);
}

const Tree = disk_tree.Folded(Colour, over);

const size = 1024;
const disks = 600;

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;
    const io = init.io;

    var args = init.minimal.args.iterate();
    _ = args.skip();
    const path = args.next() orelse "disks.ppm";

    var prng: std.Random.DefaultPrng = .init(11);
    const random = prng.random();
    var circles: [disks]disk_tree.Circle = undefined;
    var colours: [disks]Colour = undefined;
    for (&circles, &colours) |*c, *colour| {
        c.* = .{
            .cx = random.float(f64) * size,
            .cy = random.float(f64) * size,
            .r = 20 + random.float(f64) * 140,
        };
        const alpha = 0.15 + 0.35 * random.float(f32);
        colour.* = .{ alpha * random.float(f32), alpha * random.float(f32), alpha * random.float(f32), alpha };
    }

    var tree = try Tree.init(gpa, &circles, &colours);
    defer tree.deinit(gpa);

    const pixels = try gpa.alloc([3]u8, size * size);
    defer gpa.free(pixels);
    for (pixels, 0..) |*pixel, i| {
        const x: f64 = @floatFromInt(i % size);
        const y: f64 = @floatFromInt(i / size);
        const colour = tree.fold(x + 0.5, y + 0.5) orelse Colour{ 0, 0, 0, 0 };
        // Over a white background.
        for (pixel, colour[0..3]) |*channel, c| {
            channel.* = @intFromFloat(@round(255 * std.math.clamp(c + 1 - colour[3], 0, 1)));
        }
    }

    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);
    var buffer: [1 << 16]u8 = undefined;
    var writer = file.writer(io, &buffer);
    try writer.interface.print("P6\n{d} {d}\n255\n", .{ size, size });
    try writer.interface.writeAll(std.mem.sliceAsBytes(pixels));
    try writer.interface.flush();
}
