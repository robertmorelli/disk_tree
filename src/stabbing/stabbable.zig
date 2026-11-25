//! ALGORITHMS FIRST: slab-decomposed segment-tree renderer.
//! Do not replace the algorithm with brute-force point checks.
//! Preserve algorithmic structure; micro-optimizations come second.

const std = @import("std");
const bezier = @import("../geometry/bezier.zig");

pub const BoundingBox = bezier.BoundingBox;
pub const Vec2 = bezier.Vec2;

/// Y-interval at a given x coordinate
pub const YInterval = struct {
    y_min: f32,
    y_max: f32,

    pub fn init(y_min: f32, y_max: f32) YInterval {
        return .{ .y_min = y_min, .y_max = y_max };
    }

    pub fn contains(self: YInterval, y: f32) bool {
        const eps: f32 = 0.001;
        return y >= self.y_min - eps and y <= self.y_max + eps;
    }
};

/// Generic interface for objects that can be queried via stabbing.
/// Uses runtime polymorphism via vtable for flexibility.
pub const Stabbable = struct {
    ptr: *anyopaque,
    vtable: *const VTable,

    pub const VTable = struct {
        /// Get the bounding box of this shape
        boundingBox: *const fn (ptr: *anyopaque) BoundingBox,

        /// Get the Y interval at a given X coordinate.
        /// Returns null if x is outside the shape's x-range.
        yIntervalAtX: *const fn (ptr: *anyopaque, x: f32) ?YInterval,

        /// Check if a point is contained within this shape
        containsPoint: *const fn (ptr: *anyopaque, x: f32, y: f32) bool,

        /// Get the number of segments (for GPU data packing)
        segmentCount: *const fn (ptr: *anyopaque) usize,

        /// Pack segment data into a float buffer.
        /// Buffer format: [x_a, x_b, x_c, x_d, y_a, y_b, y_c, y_d, x_min, x_max, _, _]
        /// 12 floats per segment
        packSegments: *const fn (ptr: *anyopaque, buffer: []f32) void,

        /// Get shape color as RGBA
        getColor: *const fn (ptr: *anyopaque) [4]f32,

        /// Get shape z-index for stable ordering
        getZIndex: *const fn (ptr: *anyopaque) i32,

        /// Get shape translation offset
        getOffset: *const fn (ptr: *anyopaque) Vec2,
    };

    /// Get bounding box
    pub fn boundingBox(self: Stabbable) BoundingBox {
        return self.vtable.boundingBox(self.ptr);
    }

    /// Get Y interval at X
    pub fn yIntervalAtX(self: Stabbable, x: f32) ?YInterval {
        return self.vtable.yIntervalAtX(self.ptr, x);
    }

    /// Check if point is inside shape
    pub fn containsPoint(self: Stabbable, x: f32, y: f32) bool {
        return self.vtable.containsPoint(self.ptr, x, y);
    }

    /// Get segment count for GPU packing
    pub fn segmentCount(self: Stabbable) usize {
        return self.vtable.segmentCount(self.ptr);
    }

    /// Pack segments to buffer
    pub fn packSegments(self: Stabbable, buffer: []f32) void {
        self.vtable.packSegments(self.ptr, buffer);
    }

    /// Get RGBA color
    pub fn getColor(self: Stabbable) [4]f32 {
        return self.vtable.getColor(self.ptr);
    }

    /// Get z-index
    pub fn getZIndex(self: Stabbable) i32 {
        return self.vtable.getZIndex(self.ptr);
    }

    /// Get translation offset
    pub fn getOffset(self: Stabbable) Vec2 {
        return self.vtable.getOffset(self.ptr);
    }
};

/// Helper to create a Stabbable from a concrete type
pub fn makeStabbable(comptime T: type, ptr: *T) Stabbable {
    const gen = struct {
        fn boundingBox(p: *anyopaque) BoundingBox {
            const self: *T = @ptrCast(@alignCast(p));
            return self.boundingBox();
        }

        fn yIntervalAtX(p: *anyopaque, x: f32) ?YInterval {
            const self: *T = @ptrCast(@alignCast(p));
            return self.yIntervalAtX(x);
        }

        fn containsPoint(p: *anyopaque, x: f32, y: f32) bool {
            const self: *T = @ptrCast(@alignCast(p));
            return self.containsPoint(x, y);
        }

        fn segmentCount(p: *anyopaque) usize {
            const self: *T = @ptrCast(@alignCast(p));
            return self.segmentCount();
        }

        fn packSegments(p: *anyopaque, buffer: []f32) void {
            const self: *T = @ptrCast(@alignCast(p));
            self.packSegments(buffer);
        }

        fn getColor(p: *anyopaque) [4]f32 {
            const self: *T = @ptrCast(@alignCast(p));
            return self.getColor();
        }

        fn getZIndex(p: *anyopaque) i32 {
            const self: *T = @ptrCast(@alignCast(p));
            return self.getZIndex();
        }

        fn getOffset(p: *anyopaque) Vec2 {
            const self: *T = @ptrCast(@alignCast(p));
            return self.getOffset();
        }

        const vtable = Stabbable.VTable{
            .boundingBox = boundingBox,
            .yIntervalAtX = yIntervalAtX,
            .containsPoint = containsPoint,
            .segmentCount = segmentCount,
            .packSegments = packSegments,
            .getColor = getColor,
            .getZIndex = getZIndex,
            .getOffset = getOffset,
        };
    };

    return .{
        .ptr = ptr,
        .vtable = &gen.vtable,
    };
}
