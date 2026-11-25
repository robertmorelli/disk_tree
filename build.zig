const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const mod = b.addModule("disk_tree", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    const tests = b.addTest(.{ .root_module = mod });
    b.step("test", "Run the tests").dependOn(&b.addRunArtifact(tests).step);

    const bench = b.addExecutable(.{
        .name = "bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/bench.zig"),
            .target = target,
            .optimize = .ReleaseFast,
        }),
    });
    b.step("bench", "Run the benchmark").dependOn(&b.addRunArtifact(bench).step);

    const render = b.addExecutable(.{
        .name = "render",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/render.zig"),
            .target = target,
            .optimize = .ReleaseFast,
        }),
    });
    const run_render = b.addRunArtifact(render);
    if (b.args) |args| run_render.addArgs(args);
    b.step("render", "Render translucent disks or paths to a PPM image").dependOn(&run_render.step);
}
