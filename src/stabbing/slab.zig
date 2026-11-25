//! ALGORITHMS FIRST: slab-decomposed segment-tree renderer.
//! Do not replace the algorithm with brute-force point checks.
//! Preserve algorithmic structure; micro-optimizations come second.

const std = @import("std");
const Allocator = std.mem.Allocator;
const stabbable = @import("stabbable.zig");

pub const BoundingBox = stabbable.BoundingBox;
pub const YInterval = stabbable.YInterval;
pub const Stabbable = stabbable.Stabbable;

/// Segment tree node for y-interval queries
pub const TreeNode = struct {
    lo_idx: i32,
    hi_idx: i32,
    split_y: f32,
    left_child: i32,
    right_child: i32,
    bucket_start: u32,
    bucket_count: u32,

    pub fn leaf(idx: i32, y_val: f32, bucket_start: u32) TreeNode {
        return .{
            .lo_idx = idx,
            .hi_idx = idx,
            .split_y = y_val,
            .left_child = -1,
            .right_child = -1,
            .bucket_start = bucket_start,
            .bucket_count = 0,
        };
    }

    pub fn internal(lo: i32, hi: i32, split_y: f32, left: i32, right: i32, bucket_start: u32) TreeNode {
        return .{
            .lo_idx = lo,
            .hi_idx = hi,
            .split_y = split_y,
            .left_child = left,
            .right_child = right,
            .bucket_start = bucket_start,
            .bucket_count = 0,
        };
    }

    /// Pack node into 7 floats for GPU
    pub fn pack(self: TreeNode, buffer: []f32) void {
        buffer[0] = @floatFromInt(self.lo_idx);
        buffer[1] = @floatFromInt(self.hi_idx);
        buffer[2] = self.split_y;
        buffer[3] = @floatFromInt(self.left_child);
        buffer[4] = @floatFromInt(self.right_child);
        buffer[5] = @floatFromInt(self.bucket_start);
        buffer[6] = @floatFromInt(self.bucket_count);
    }
};

/// Endpoint for segment tree construction
const Endpoint = struct {
    y: f32,
    shape_local_idx: u32,
    is_upper: bool, // false = lower, true = upper
};

fn endpointLessThan(_: void, a: Endpoint, b: Endpoint) bool {
    return a.y < b.y;
}

/// A vertical slab with cached y-intervals for active shapes
pub const Slab = struct {
    allocator: Allocator,
    x_left: f32,
    x_right: f32,
    x_mid: f32,

    /// Active shape indices (global indices into the shapes array)
    active_shapes: std.ArrayList(u32),

    /// Y-intervals for each active shape at x_mid
    y_intervals: std.ArrayList(YInterval),

    /// Segment tree nodes (flattened)
    tree_nodes: std.ArrayList(TreeNode),

    /// Shape indices stored in tree buckets
    buckets: std.ArrayList(u32),

    pub fn init(allocator: Allocator, x_left: f32, x_right: f32) Slab {
        return .{
            .allocator = allocator,
            .x_left = x_left,
            .x_right = x_right,
            .x_mid = (x_left + x_right) / 2.0,
            .active_shapes = .empty,
            .y_intervals = .empty,
            .tree_nodes = .empty,
            .buckets = .empty,
        };
    }

    pub fn deinit(self: *Slab) void {
        self.active_shapes.deinit(self.allocator);
        self.y_intervals.deinit(self.allocator);
        self.tree_nodes.deinit(self.allocator);
        self.buckets.deinit(self.allocator);
    }

    /// Reset slab contents while retaining backing capacity for reuse.
    pub fn reset(self: *Slab, x_left: f32, x_right: f32) void {
        self.x_left = x_left;
        self.x_right = x_right;
        self.x_mid = (x_left + x_right) / 2.0;
        self.active_shapes.clearRetainingCapacity();
        self.y_intervals.clearRetainingCapacity();
        self.tree_nodes.clearRetainingCapacity();
        self.buckets.clearRetainingCapacity();
    }

    /// Build the slab's segment tree from the active shapes
    pub fn buildTree(self: *Slab) !void {
        if (self.active_shapes.items.len == 0) {
            // Empty slab - add a dummy node
            try self.tree_nodes.append(self.allocator, TreeNode.leaf(0, 0, 0));
            try self.buckets.append(self.allocator, 0);
            return;
        }

        // Create endpoints list sorted by y
        var endpoints: std.ArrayList(Endpoint) = .empty;
        defer endpoints.deinit(self.allocator);

        for (self.y_intervals.items, 0..) |interval, i| {
            try endpoints.append(self.allocator, .{
                .y = interval.y_min,
                .shape_local_idx = @intCast(i),
                .is_upper = false,
            });
            try endpoints.append(self.allocator, .{
                .y = interval.y_max,
                .shape_local_idx = @intCast(i),
                .is_upper = true,
            });
        }

        std.mem.sort(Endpoint, endpoints.items, {}, endpointLessThan);

        if (endpoints.items.len == 0) {
            try self.tree_nodes.append(self.allocator, TreeNode.leaf(0, 0, 0));
            try self.buckets.append(self.allocator, 0);
            return;
        }

        // Build tree recursively
        _ = try self.buildNode(endpoints.items, 0, @as(i32, @intCast(endpoints.items.len)) - 1);

        // Find ranks for each interval's endpoints
        const ranks = try self.allocator.alloc([2]i32, self.active_shapes.items.len);
        defer self.allocator.free(ranks);
        for (endpoints.items, 0..) |ep, rank| {
            ranks[ep.shape_local_idx][@intFromBool(ep.is_upper)] = @intCast(rank);
        }

        // Insert intervals into tree. Each node's bucket is one contiguous
        // range of `buckets`, so count every bucket before filling any.
        for (self.active_shapes.items, ranks) |shape_idx, rank| {
            self.insertInterval(0, @min(rank[0], rank[1]), @max(rank[0], rank[1]), shape_idx, .count);
        }
        var bucket_total: u32 = 0;
        for (self.tree_nodes.items) |*node| {
            node.bucket_start = bucket_total;
            bucket_total += node.bucket_count;
            node.bucket_count = 0;
        }
        try self.buckets.resize(self.allocator, bucket_total);
        for (self.active_shapes.items, ranks) |shape_idx, rank| {
            self.insertInterval(0, @min(rank[0], rank[1]), @max(rank[0], rank[1]), shape_idx, .fill);
        }

        // Add dummy bucket if empty
        if (self.buckets.items.len == 0) {
            try self.buckets.append(self.allocator, 0);
        }
    }

    fn buildNode(self: *Slab, endpoints: []Endpoint, lo: i32, hi: i32) !i32 {
        const node_idx: i32 = @intCast(self.tree_nodes.items.len);

        if (lo == hi) {
            const y_val = endpoints[@intCast(lo)].y;
            try self.tree_nodes.append(self.allocator, TreeNode.leaf(lo, y_val, @intCast(self.buckets.items.len)));
            return node_idx;
        }

        // Reserve space for this node
        try self.tree_nodes.append(self.allocator, undefined);

        const mid = @divFloor(lo + hi, 2);
        const split_y = endpoints[@intCast(mid)].y;

        const left_child = try self.buildNode(endpoints, lo, mid);
        const right_child = try self.buildNode(endpoints, mid + 1, hi);

        self.tree_nodes.items[@intCast(node_idx)] = TreeNode.internal(
            lo,
            hi,
            split_y,
            left_child,
            right_child,
            @intCast(self.buckets.items.len),
        );

        return node_idx;
    }

    const InsertPass = enum { count, fill };

    fn insertInterval(self: *Slab, node_idx: i32, lo_rank: i32, hi_rank: i32, shape_idx: u32, pass: InsertPass) void {
        if (node_idx < 0 or node_idx >= @as(i32, @intCast(self.tree_nodes.items.len))) {
            return;
        }

        var node = &self.tree_nodes.items[@intCast(node_idx)];
        const node_lo = node.lo_idx;
        const node_hi = node.hi_idx;

        if (lo_rank <= node_lo and node_hi <= hi_rank) {
            // Fully covered - add to bucket
            if (pass == .fill) {
                self.buckets.items[node.bucket_start + node.bucket_count] = shape_idx;
            }
            node.bucket_count += 1;
            return;
        }

        const mid = @divFloor(node_lo + node_hi, 2);
        const left_child = node.left_child;
        const right_child = node.right_child;

        if (lo_rank <= mid and left_child >= 0) {
            self.insertInterval(left_child, lo_rank, hi_rank, shape_idx, pass);
        }
        if (hi_rank > mid and right_child >= 0) {
            self.insertInterval(right_child, lo_rank, hi_rank, shape_idx, pass);
        }
    }

    /// Query the slab for shapes potentially containing point with y-coordinate y
    pub fn queryCandidates(self: *const Slab, y: f32, allocator: Allocator, candidates: *std.ArrayList(u32)) !void {
        if (self.tree_nodes.items.len == 0) return;

        // Stack-based traversal
        var stack: [32]i32 = undefined;
        var stack_ptr: usize = 0;
        stack[stack_ptr] = 0;
        stack_ptr += 1;

        while (stack_ptr > 0) {
            stack_ptr -= 1;
            const node_idx = stack[stack_ptr];
            if (node_idx < 0) continue;

            const node = self.tree_nodes.items[@intCast(node_idx)];

            // Add bucket contents
            var i: u32 = 0;
            while (i < node.bucket_count) : (i += 1) {
                const bucket_idx = node.bucket_start + i;
                if (bucket_idx < self.buckets.items.len) {
                    try candidates.append(allocator, self.buckets.items[bucket_idx]);
                }
            }

            // Descend
            if (node.left_child < 0 and node.right_child < 0) continue;

            if (y <= node.split_y) {
                if (node.left_child >= 0 and stack_ptr < stack.len) {
                    stack[stack_ptr] = node.left_child;
                    stack_ptr += 1;
                }
            } else {
                if (node.right_child >= 0 and stack_ptr < stack.len) {
                    stack[stack_ptr] = node.right_child;
                    stack_ptr += 1;
                }
            }
        }
    }
};

// Tests
test "slab basic" {
    const allocator = std.testing.allocator;

    var slab = Slab.init(allocator, 0.0, 1.0);
    defer slab.deinit();

    // Add some active shapes
    try slab.active_shapes.append(allocator, 0);
    try slab.active_shapes.append(allocator, 1);

    try slab.y_intervals.append(allocator, YInterval.init(0.0, 1.0));
    try slab.y_intervals.append(allocator, YInterval.init(0.5, 1.5));

    try slab.buildTree();

    // Query
    var candidates: std.ArrayList(u32) = .empty;
    defer candidates.deinit(allocator);

    try slab.queryCandidates(0.75, allocator, &candidates);
    try std.testing.expect(candidates.items.len > 0);
}
