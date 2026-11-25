//! ALGORITHMS FIRST: slab-decomposed segment-tree renderer.
//! Do not replace the algorithm with brute-force point checks.
//! Preserve algorithmic structure; micro-optimizations come second.

const std = @import("std");
const Allocator = std.mem.Allocator;
const stabbable = @import("stabbable.zig");
const slab_mod = @import("slab.zig");

pub const BoundingBox = stabbable.BoundingBox;
pub const YInterval = stabbable.YInterval;
pub const Stabbable = stabbable.Stabbable;
pub const Slab = slab_mod.Slab;
pub const TreeNode = slab_mod.TreeNode;

fn nanoTimestamp() i128 {
    return std.Io.Timestamp.now(std.Io.Threaded.global_single_threaded.io(), .real).nanoseconds;
}

pub const BinMode = enum {
    lists,
    range_filter,
    premerge,
};

pub const TreeConfig = struct {
    slab_target: u32 = 1024,
    num_bins: u32 = 8,
    wide_threshold: u32 = 0,
    bin_mode: BinMode = .premerge,
    pack_retain: bool = true,
    pack_reuse_bins: bool = true,
    pack_reuse_shape_bins: bool = true,
    zsort_check: bool = true,
};

fn normalizeConfig(cfg: *TreeConfig) void {
    if (cfg.num_bins == 0) {
        cfg.wide_threshold = 0;
    } else {
        if (cfg.num_bins > 64) cfg.num_bins = 64;
        if (cfg.wide_threshold == 0 or cfg.wide_threshold > cfg.num_bins) {
            cfg.wide_threshold = cfg.num_bins;
        }
    }
}

/// Stabbing tree for efficient point-in-shape queries
pub const StabbingTree = struct {
    pub const BuildProfile = struct {
        total_ms: f64 = 0.0,
        bbox_ms: f64 = 0.0,
        xlist_ms: f64 = 0.0,
        sortx_ms: f64 = 0.0,
        unique_ms: f64 = 0.0,
        events_ms: f64 = 0.0,
        sort_events_ms: f64 = 0.0,
        zorder_ms: f64 = 0.0,
        slabs_ms: f64 = 0.0,
        coalesce_ms: f64 = 0.0,
        slabs_built: u32 = 0,
    };

    pub const PackProfile = struct {
        total_ms: f64 = 0.0,
        segments_ms: f64 = 0.0,
        colors_ms: f64 = 0.0,
        bboxes_ms: f64 = 0.0,
        offsets_ms: f64 = 0.0,
        slabs_ms: f64 = 0.0,
        slabs_count_ms: f64 = 0.0,
        slabs_prefix_ms: f64 = 0.0,
        slabs_global_ms: f64 = 0.0,
        slabs_bin_ms: f64 = 0.0,
    };

    allocator: Allocator,
    slabs: std.ArrayList(Slab),
    shapes: []Stabbable,
    slab_target: u32,
    num_bins: u32,
    wide_threshold: u32,
    bin_mode: BinMode,
    pack_retain: bool,
    pack_reuse_bins: bool,
    pack_reuse_shape_bins: bool,
    zsort_check: bool,
    profile_enabled: bool,
    last_build_profile: BuildProfile,
    last_pack_profile: PackProfile,

    pub fn init(allocator: Allocator) StabbingTree {
        return initWithConfig(allocator, TreeConfig{});
    }

    pub fn initWithConfig(allocator: Allocator, config: TreeConfig) StabbingTree {
        var cfg = config;
        normalizeConfig(&cfg);
        return .{
            .allocator = allocator,
            .slabs = .empty,
            .shapes = &[_]Stabbable{},
            .slab_target = cfg.slab_target,
            .num_bins = cfg.num_bins,
            .wide_threshold = cfg.wide_threshold,
            .bin_mode = cfg.bin_mode,
            .pack_retain = cfg.pack_retain,
            .pack_reuse_bins = cfg.pack_reuse_bins,
            .pack_reuse_shape_bins = cfg.pack_reuse_shape_bins,
            .zsort_check = cfg.zsort_check,
            .profile_enabled = false,
            .last_build_profile = .{},
            .last_pack_profile = .{},
        };
    }

    pub fn deinit(self: *StabbingTree) void {
        for (self.slabs.items) |*slab| {
            slab.deinit();
        }
        self.slabs.deinit(self.allocator);
    }

    /// Build the stabbing tree from a list of shapes
    pub fn build(self: *StabbingTree, shapes: []Stabbable) !void {
        var t_start: i128 = 0;
        var t_last: i128 = 0;
        if (self.profile_enabled) {
            t_start = nanoTimestamp();
            t_last = t_start;
            self.last_build_profile = .{};
        }

        self.shapes = shapes;

        if (shapes.len == 0) return;

        // Collect all x-boundaries (optionally quantized)
        var x_list: std.ArrayList(f32) = .empty;
        defer x_list.deinit(self.allocator);

        var global_min: f32 = std.math.floatMax(f32);
        var global_max: f32 = -std.math.floatMax(f32);
        for (shapes) |shape| {
            const bbox = shape.boundingBox();
            global_min = @min(global_min, bbox.x_min);
            global_max = @max(global_max, bbox.x_max);
        }
        if (self.profile_enabled) {
            const now = nanoTimestamp();
            self.last_build_profile.bbox_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            t_last = now;
        }

        const range = global_max - global_min;
        var do_quantize = self.slab_target > 0 and range > 1e-6;
        var step: f32 = 0.0;
        if (do_quantize) {
            step = range / @as(f32, @floatFromInt(self.slab_target));
            if (step <= 1e-6) {
                do_quantize = false;
            }
        }

        for (shapes) |shape| {
            const bbox = shape.boundingBox();
            var x_min = bbox.x_min;
            var x_max = bbox.x_max;

            if (do_quantize) {
                const tmin = (x_min - global_min) / step;
                const tmax = (x_max - global_min) / step;
                x_min = global_min + @floor(tmin) * step;
                x_max = global_min + @ceil(tmax) * step;

                if (x_max <= x_min) {
                    x_max = x_min + step;
                    if (x_max > global_max) {
                        x_min = x_max - step;
                    }
                }
                x_min = std.math.clamp(x_min, global_min, global_max);
                x_max = std.math.clamp(x_max, global_min, global_max);
            }

            try x_list.append(self.allocator, x_min);
            try x_list.append(self.allocator, x_max);
        }
        if (self.profile_enabled) {
            const now = nanoTimestamp();
            self.last_build_profile.xlist_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            t_last = now;
        }

        // Note: We skip expensive O(n²) intersection finding.
        // Bbox boundaries are sufficient - coalescing merges redundant slabs,
        // and shape_contains_point handles actual geometry.

        // Sort and deduplicate x coordinates
        std.mem.sort(f32, x_list.items, {}, std.sort.asc(f32));
        if (self.profile_enabled) {
            const now = nanoTimestamp();
            self.last_build_profile.sortx_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            t_last = now;
        }

        var unique_xs: std.ArrayList(f32) = .empty;
        defer unique_xs.deinit(self.allocator);

        var prev: f32 = -std.math.floatMax(f32);
        for (x_list.items) |x| {
            if (x - prev > 1e-9) {
                try unique_xs.append(self.allocator, x);
                prev = x;
            }
        }
        if (self.profile_enabled) {
            const now = nanoTimestamp();
            self.last_build_profile.unique_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            t_last = now;
        }

        // Build sweep-line events: (x, shape_idx, is_start)
        const Event = struct {
            x: f32,
            shape_idx: u32,
            is_start: bool,

            fn lessThan(_: void, a: @This(), b: @This()) bool {
                if (a.x != b.x) return a.x < b.x;
                // Starts before ends at same x
                if (a.is_start != b.is_start) return a.is_start;
                return a.shape_idx < b.shape_idx;
            }
        };

        var events: std.ArrayList(Event) = .empty;
        defer events.deinit(self.allocator);

        for (shapes, 0..) |shape, j| {
            const bbox = shape.boundingBox();
            var x_min = bbox.x_min;
            var x_max = bbox.x_max;
            if (do_quantize) {
                const tmin = (x_min - global_min) / step;
                const tmax = (x_max - global_min) / step;
                x_min = global_min + @floor(tmin) * step;
                x_max = global_min + @ceil(tmax) * step;
                if (x_max <= x_min) {
                    x_max = x_min + step;
                    if (x_max > global_max) {
                        x_min = x_max - step;
                    }
                }
                x_min = std.math.clamp(x_min, global_min, global_max);
                x_max = std.math.clamp(x_max, global_min, global_max);
            }

            try events.append(self.allocator, .{ .x = x_min, .shape_idx = @intCast(j), .is_start = true });
            try events.append(self.allocator, .{ .x = x_max, .shape_idx = @intCast(j), .is_start = false });
        }
        if (self.profile_enabled) {
            const now = nanoTimestamp();
            self.last_build_profile.events_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            t_last = now;
        }

        std.mem.sort(Event, events.items, {}, Event.lessThan);
        if (self.profile_enabled) {
            const now = nanoTimestamp();
            self.last_build_profile.sort_events_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            t_last = now;
        }

        const n = shapes.len;

        // Build z-order mapping once per build
        const z_order_to_shape = try self.allocator.alloc(u32, n);
        defer self.allocator.free(z_order_to_shape);
        for (z_order_to_shape, 0..) |*slot, idx| slot.* = @intCast(idx);

        const SortCtx = struct {
            shapes: []Stabbable,
        };
        const sort_ctx = SortCtx{ .shapes = shapes };
        const lessThan = struct {
            fn lt(ctx: SortCtx, a: u32, b: u32) bool {
                const za = ctx.shapes[a].getZIndex();
                const zb = ctx.shapes[b].getZIndex();
                if (za != zb) return za > zb; // front-to-back
                return a > b;
            }
        }.lt;
        var already_sorted = false;
        if (self.zsort_check and n > 1) {
            already_sorted = true;
            var prev_idx: u32 = z_order_to_shape[0];
            var prev_z = sort_ctx.shapes[prev_idx].getZIndex();
            var i: usize = 1;
            while (i < n) : (i += 1) {
                const cur_idx: u32 = z_order_to_shape[i];
                const cur_z = sort_ctx.shapes[cur_idx].getZIndex();
                if (cur_z > prev_z or (cur_z == prev_z and cur_idx > prev_idx)) {
                    already_sorted = false;
                    break;
                }
                prev_z = cur_z;
                prev_idx = cur_idx;
            }
        }
        if (!already_sorted) {
            std.mem.sort(u32, z_order_to_shape, sort_ctx, lessThan);
        }

        var shape_to_zorder = try self.allocator.alloc(u32, n);
        defer self.allocator.free(shape_to_zorder);
        for (z_order_to_shape, 0..) |shape_idx, order| {
            shape_to_zorder[shape_idx] = @intCast(order);
        }
        if (self.profile_enabled) {
            const now = nanoTimestamp();
            self.last_build_profile.zorder_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            t_last = now;
        }

        // Active set as bitset in z-order space
        const word_count = (n + 63) / 64;
        var active_bits = try self.allocator.alloc(u64, word_count);
        defer self.allocator.free(active_bits);
        @memset(active_bits, 0);

        var event_idx: usize = 0;

        var i: usize = 0;
        while (i + 1 < unique_xs.items.len) : (i += 1) {
            const x_left = unique_xs.items[i];
            const x_right = unique_xs.items[i + 1];

            if (x_right <= x_left + 1e-12) continue;

            // Process events up to x_left (inclusive for starts)
            while (event_idx < events.items.len and events.items[event_idx].x <= x_left) {
                const ev = events.items[event_idx];
                if (ev.is_start) {
                    const z_idx = shape_to_zorder[ev.shape_idx];
                    const word_idx = z_idx / 64;
                    const bit: u6 = @intCast(z_idx & 63);
                    active_bits[word_idx] |= (@as(u64, 1) << bit);
                } else {
                    const z_idx = shape_to_zorder[ev.shape_idx];
                    const word_idx = z_idx / 64;
                    const bit: u6 = @intCast(z_idx & 63);
                    active_bits[word_idx] &= ~(@as(u64, 1) << bit);
                }
                event_idx += 1;
            }

            var slab = Slab.init(self.allocator, x_left, x_right);
            errdefer slab.deinit();

            // Collect active shapes in z-order
            var active_total: usize = 0;
            for (active_bits) |word| {
                active_total += @popCount(word);
            }
            try slab.active_shapes.ensureTotalCapacity(self.allocator, active_total);

            var word_idx: usize = 0;
            while (word_idx < active_bits.len) : (word_idx += 1) {
                var word = active_bits[word_idx];
                while (word != 0) {
                    const bit = @ctz(word);
                    const z_idx = word_idx * 64 + bit;
                    if (z_idx < n) {
                        try slab.active_shapes.append(self.allocator, z_order_to_shape[z_idx]);
                    }
                    word &= word - 1;
                }
            }

            try self.slabs.append(self.allocator, slab);
        }
        if (self.profile_enabled) {
            const now = nanoTimestamp();
            self.last_build_profile.slabs_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            self.last_build_profile.slabs_built = @intCast(self.slabs.items.len);
            t_last = now;
        }

        // Coalesce adjacent slabs with identical active shapes
        try self.coalesceSlabs();

        // Build each slab's segment tree. A shape's y-interval changes across
        // the slab, so the tree is built over each shape's full y-extent; the
        // containment test in query() makes the answer exact.
        for (self.slabs.items) |*slab| {
            try slab.y_intervals.ensureTotalCapacity(self.allocator, slab.active_shapes.items.len);
            for (slab.active_shapes.items) |shape_idx| {
                const bbox = shapes[shape_idx].boundingBox();
                slab.y_intervals.appendAssumeCapacity(YInterval.init(bbox.y_min, bbox.y_max));
            }
            try slab.buildTree();
        }
        if (self.profile_enabled) {
            const now = nanoTimestamp();
            self.last_build_profile.coalesce_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            self.last_build_profile.total_ms = @as(f64, @floatFromInt(now - t_start)) / 1_000_000.0;
        }
    }

    /// Merge adjacent slabs that have identical active shape sets
    fn coalesceSlabs(self: *StabbingTree) !void {
        if (self.slabs.items.len < 2) return;

        var write_idx: usize = 0;
        var read_idx: usize = 1;

        while (read_idx < self.slabs.items.len) {
            const current = &self.slabs.items[write_idx];
            const next = &self.slabs.items[read_idx];

            // Check if active shapes are identical
            if (slabsHaveSameActiveShapes(current, next)) {
                // Merge: extend current slab's x_right, free next slab
                current.x_right = next.x_right;
                current.x_mid = (current.x_left + current.x_right) / 2.0;
                next.deinit();
            } else {
                // Move to next position
                write_idx += 1;
                if (write_idx != read_idx) {
                    self.slabs.items[write_idx] = self.slabs.items[read_idx];
                }
            }
            read_idx += 1;
        }

        // Shrink the list
        self.slabs.shrinkRetainingCapacity(write_idx + 1);
    }

    /// Check if two slabs have identical active shape sets
    fn slabsHaveSameActiveShapes(a: *const Slab, b: *const Slab) bool {
        if (a.active_shapes.items.len != b.active_shapes.items.len) return false;
        for (a.active_shapes.items, b.active_shapes.items) |ai, bi| {
            if (ai != bi) return false;
        }
        return true;
    }

    /// Find intersection x-coordinates between shape boundaries
    fn findIntersections(self: *StabbingTree, shapes: []Stabbable, x_list: *std.ArrayList(f32)) !void {
        // Simplified approach: for overlapping bounding boxes, sample at regular intervals
        const n = shapes.len;
        if (n < 2) return;

        var i: usize = 0;
        while (i < n) : (i += 1) {
            const bbox_i = shapes[i].boundingBox();

            var j: usize = i + 1;
            while (j < n) : (j += 1) {
                const bbox_j = shapes[j].boundingBox();

                // Check bounding box overlap
                if (!bbox_i.overlaps(bbox_j)) continue;

                // Overlapping boxes - add boundary points and sample x values
                const x_min = @max(bbox_i.x_min, bbox_j.x_min);
                const x_max = @min(bbox_i.x_max, bbox_j.x_max);

                if (x_max <= x_min) continue;

                // Sample for potential intersection points
                const num_samples: u32 = 10;
                var k: u32 = 0;
                while (k <= num_samples) : (k += 1) {
                    const t = @as(f32, @floatFromInt(k)) / @as(f32, @floatFromInt(num_samples));
                    const x = x_min + t * (x_max - x_min);

                    const interval_i = shapes[i].yIntervalAtX(x);
                    const interval_j = shapes[j].yIntervalAtX(x);

                    if (interval_i != null and interval_j != null) {
                        const yi = interval_i.?;
                        const yj = interval_j.?;

                        // Check for y-interval overlap change (potential crossing)
                        const overlaps = !(yi.y_max < yj.y_min or yj.y_max < yi.y_min);
                        if (overlaps) {
                            try x_list.append(self.allocator, x);
                        }
                    }
                }
            }
        }
    }

    /// Query for shapes containing point (x, y)
    pub fn query(self: *const StabbingTree, x: f32, y: f32, results: *std.ArrayList(u32)) !void {
        if (self.slabs.items.len == 0) return;

        // Find slab by binary search
        var lo: usize = 0;
        var hi: usize = self.slabs.items.len;

        while (lo < hi) {
            const mid = lo + (hi - lo) / 2;
            if (self.slabs.items[mid].x_left <= x) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }

        const slab_idx = if (lo > 0) lo - 1 else 0;
        if (slab_idx >= self.slabs.items.len) return;

        const slab = &self.slabs.items[slab_idx];

        // Get candidates from segment tree
        var candidates: std.ArrayList(u32) = .empty;
        defer candidates.deinit(self.allocator);

        try slab.queryCandidates(y, self.allocator, &candidates);

        // Final containment test
        var seen = std.AutoHashMap(u32, void).init(self.allocator);
        defer seen.deinit();

        for (candidates.items) |shape_idx| {
            if (seen.contains(shape_idx)) continue;
            try seen.put(shape_idx, {});

            if (shape_idx < self.shapes.len) {
                if (self.shapes[shape_idx].containsPoint(x, y)) {
                    try results.append(self.allocator, shape_idx);
                }
            }
        }
    }

    /// Get number of slabs
    pub fn slabCount(self: *const StabbingTree) usize {
        return self.slabs.items.len;
    }
};

/// GPU-friendly packed data for the stabbing tree
pub const PackedTreeData = struct {
    allocator: Allocator,

    /// Per-segment data: [x_a, x_b, x_c, x_d, y_a, y_b, y_c, y_d, x_min, x_max, _, _]
    segment_data: std.ArrayList(f32),

    /// Per-shape: [segment_start, segment_count]
    shape_segments: std.ArrayList(u32),

    /// Per-shape RGBA color
    colors: std.ArrayList(f32),

    /// Per-shape bounding box: [x_min, y_min, x_max, y_max]
    shape_bboxes: std.ArrayList(f32),

    /// Per-shape translation offset + z-order: [x, y, z_order, bin_min, bin_max, _]
    shape_offsets: std.ArrayList(f32),

    /// Left x of each slab
    slab_xs: std.ArrayList(f32),

    /// Per-slab info. If num_bins == 0: [active_start, active_count, tree_root, _]
    /// If num_bins > 0: [global_start, global_count] + per-bin [start, count] packed consecutively for each slab.
    slab_info: std.ArrayList(i32),

    /// All active indices (or bin indices when num_bins > 0) concatenated
    active_indices: std.ArrayList(u32),
    /// Max active list size packed this frame (for shader loop cap)
    max_active: u32,
    /// Draw order for instanced rendering (back-to-front)
    draw_order: std.ArrayList(u32),

    /// Optional reusable bin scratch buffers
    bin_counts_buf: std.ArrayList(u32),
    bin_offsets_buf: std.ArrayList(u32),
    node_write_buf: std.ArrayList(u32),
    node_lists_buf: std.ArrayList(u32),
    bin_path_nodes_buf: std.ArrayList(u32),
    bin_path_len_buf: std.ArrayList(u8),
    slab_node_counts_buf: std.ArrayList(u32),
    slab_bin_counts_buf: std.ArrayList(u32),
    slab_bin_start_buf: std.ArrayList(u32),
    slab_total_active_buf: std.ArrayList(u32),
    slab_base_buf: std.ArrayList(u32),
    shape_bin_min_buf: std.ArrayList(u32),
    shape_bin_max_buf: std.ArrayList(u32),

    pub fn init(allocator: Allocator) PackedTreeData {
        return .{
            .allocator = allocator,
            .segment_data = .empty,
            .shape_segments = .empty,
            .colors = .empty,
            .shape_bboxes = .empty,
            .shape_offsets = .empty,
            .slab_xs = .empty,
            .slab_info = .empty,
            .active_indices = .empty,
            .max_active = 0,
            .draw_order = .empty,
            .bin_counts_buf = .empty,
            .bin_offsets_buf = .empty,
            .node_write_buf = .empty,
            .node_lists_buf = .empty,
            .bin_path_nodes_buf = .empty,
            .bin_path_len_buf = .empty,
            .slab_node_counts_buf = .empty,
            .slab_bin_counts_buf = .empty,
            .slab_bin_start_buf = .empty,
            .slab_total_active_buf = .empty,
            .slab_base_buf = .empty,
            .shape_bin_min_buf = .empty,
            .shape_bin_max_buf = .empty,
        };
    }

    pub fn deinit(self: *PackedTreeData) void {
        self.segment_data.deinit(self.allocator);
        self.shape_segments.deinit(self.allocator);
        self.colors.deinit(self.allocator);
        self.shape_bboxes.deinit(self.allocator);
        self.shape_offsets.deinit(self.allocator);
        self.slab_xs.deinit(self.allocator);
        self.slab_info.deinit(self.allocator);
        self.active_indices.deinit(self.allocator);
        self.draw_order.deinit(self.allocator);
        self.bin_counts_buf.deinit(self.allocator);
        self.bin_offsets_buf.deinit(self.allocator);
        self.node_write_buf.deinit(self.allocator);
        self.node_lists_buf.deinit(self.allocator);
        self.bin_path_nodes_buf.deinit(self.allocator);
        self.bin_path_len_buf.deinit(self.allocator);
        self.slab_node_counts_buf.deinit(self.allocator);
        self.slab_bin_counts_buf.deinit(self.allocator);
        self.slab_bin_start_buf.deinit(self.allocator);
        self.slab_total_active_buf.deinit(self.allocator);
        self.slab_base_buf.deinit(self.allocator);
        self.shape_bin_min_buf.deinit(self.allocator);
        self.shape_bin_max_buf.deinit(self.allocator);
    }

    /// Pack tree and shapes into GPU-ready buffers
    pub fn pack(self: *PackedTreeData, tree: *const StabbingTree, shapes: []Stabbable, y_min: f32, y_max: f32) !void {
        var t_start: i128 = 0;
        var t_last: i128 = 0;
        if (tree.profile_enabled) {
            t_start = nanoTimestamp();
            t_last = t_start;
            @constCast(tree).last_pack_profile = .{};
        }

        // Clear existing data (optional capacity retention to avoid reallocs)
        if (tree.pack_retain) {
            self.segment_data.clearRetainingCapacity();
            self.shape_segments.clearRetainingCapacity();
            self.colors.clearRetainingCapacity();
            self.shape_bboxes.clearRetainingCapacity();
            self.shape_offsets.clearRetainingCapacity();
            self.slab_xs.clearRetainingCapacity();
            self.slab_info.clearRetainingCapacity();
            self.active_indices.clearRetainingCapacity();
            self.draw_order.clearRetainingCapacity();
        } else {
            self.segment_data.clearAndFree(self.allocator);
            self.shape_segments.clearAndFree(self.allocator);
            self.colors.clearAndFree(self.allocator);
            self.shape_bboxes.clearAndFree(self.allocator);
            self.shape_offsets.clearAndFree(self.allocator);
            self.slab_xs.clearAndFree(self.allocator);
            self.slab_info.clearAndFree(self.allocator);
            self.active_indices.clearAndFree(self.allocator);
            self.draw_order.clearAndFree(self.allocator);
        }
        self.max_active = 0;

        const n_shapes: usize = shapes.len;
        const num_bins: u32 = if (tree.slabs.items.len == 0) 0 else tree.num_bins;
        const z_order_to_shape = try self.allocator.alloc(u32, n_shapes);
        defer self.allocator.free(z_order_to_shape);
        for (z_order_to_shape, 0..) |*slot, idx| slot.* = @intCast(idx);

        const SortCtx = struct {
            shapes: []Stabbable,
        };
        const sort_ctx = SortCtx{ .shapes = shapes };
        const lessThan = struct {
            fn lt(ctx: SortCtx, a: u32, b: u32) bool {
                const za = ctx.shapes[a].getZIndex();
                const zb = ctx.shapes[b].getZIndex();
                if (za != zb) return za > zb; // front-to-back
                return a > b;
            }
        }.lt;
        std.mem.sort(u32, z_order_to_shape, sort_ctx, lessThan);

        // Draw order for instanced rendering (back-to-front)
        var i_order: usize = n_shapes;
        while (i_order > 0) {
            i_order -= 1;
            try self.draw_order.append(self.allocator, z_order_to_shape[i_order]);
        }

        const shape_to_zorder = try self.allocator.alloc(u32, n_shapes);
        defer self.allocator.free(shape_to_zorder);
        for (z_order_to_shape, 0..) |shape_idx, order| {
            shape_to_zorder[shape_idx] = @intCast(order);
        }

        // Pack segments + shape_segments
        for (shapes) |shape| {
            const seg_count = shape.segmentCount();
            const seg_start: u32 = @intCast(self.segment_data.items.len / 12);

            try self.shape_segments.append(self.allocator, seg_start);
            try self.shape_segments.append(self.allocator, @intCast(seg_count));

            // Allocate space for segments
            const old_len = self.segment_data.items.len;
            try self.segment_data.resize(self.allocator, old_len + seg_count * 12);
            shape.packSegments(self.segment_data.items[old_len..]);
        }
        if (tree.profile_enabled) {
            const now = nanoTimestamp();
            const prof = &@constCast(tree).last_pack_profile;
            prof.segments_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            t_last = now;
        }

        // Pack colors
        for (shapes) |shape| {
            const color = shape.getColor();
            try self.colors.append(self.allocator, color[0]);
            try self.colors.append(self.allocator, color[1]);
            try self.colors.append(self.allocator, color[2]);
            try self.colors.append(self.allocator, color[3]);
        }
        if (tree.profile_enabled) {
            const now = nanoTimestamp();
            const prof = &@constCast(tree).last_pack_profile;
            prof.colors_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            t_last = now;
        }

        // Pack bounding boxes
        for (shapes) |shape| {
            const bbox = shape.boundingBox();
            try self.shape_bboxes.append(self.allocator, bbox.x_min);
            try self.shape_bboxes.append(self.allocator, bbox.y_min);
            try self.shape_bboxes.append(self.allocator, bbox.x_max);
            try self.shape_bboxes.append(self.allocator, bbox.y_max);
        }
        if (tree.profile_enabled) {
            const now = nanoTimestamp();
            const prof = &@constCast(tree).last_pack_profile;
            prof.bboxes_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            t_last = now;
        }

        var shape_bin_min: []u32 = &[_]u32{};
        var shape_bin_max: []u32 = &[_]u32{};
        var have_bins = false;
        var free_shape_bins = false;
        defer if (free_shape_bins) {
            if (shape_bin_min.len != 0) self.allocator.free(shape_bin_min);
            if (shape_bin_max.len != 0) self.allocator.free(shape_bin_max);
        };
        if (num_bins > 0) {
            const y_range = y_max - y_min;
            if (y_range <= 1e-6) {
                return error.InvalidRange;
            }
            const bin_height = y_range / @as(f32, @floatFromInt(num_bins));
            const inv_bin_height = 1.0 / bin_height;

            have_bins = true;
            if (tree.pack_reuse_shape_bins) {
                try self.shape_bin_min_buf.resize(self.allocator, n_shapes);
                try self.shape_bin_max_buf.resize(self.allocator, n_shapes);
                shape_bin_min = self.shape_bin_min_buf.items;
                shape_bin_max = self.shape_bin_max_buf.items;
            } else {
                shape_bin_min = try self.allocator.alloc(u32, n_shapes);
                shape_bin_max = try self.allocator.alloc(u32, n_shapes);
                free_shape_bins = true;
            }

            const max_bin_i: i32 = @intCast(num_bins - 1);
            // Pad by one bin to avoid missing shapes due to rounding/tolerance.
            const bin_pad: i32 = 1;
            var idx: usize = 0;
            while (idx < n_shapes) : (idx += 1) {
                const base = idx * 4;
                const bb_y_min = self.shape_bboxes.items[base + 1];
                const bb_y_max = self.shape_bboxes.items[base + 3];
                const bin_min_f = (bb_y_min - y_min) * inv_bin_height;
                const bin_max_f = (bb_y_max - y_min) * inv_bin_height;
                var bin_min: i32 = @intFromFloat(@floor(bin_min_f));
                var bin_max: i32 = @intFromFloat(@floor(bin_max_f));
                bin_min -= bin_pad;
                bin_max += bin_pad;
                if (bin_min < 0) bin_min = 0;
                if (bin_max < 0) bin_max = 0;
                if (bin_min > max_bin_i) bin_min = max_bin_i;
                if (bin_max > max_bin_i) bin_max = max_bin_i;
                if (bin_max < bin_min) {
                    const tmp = bin_min;
                    bin_min = bin_max;
                    bin_max = tmp;
                }
                shape_bin_min[idx] = @intCast(bin_min);
                shape_bin_max[idx] = @intCast(bin_max);
            }
        }

        // Pack offsets + z-order + bin ranges
        for (shapes, 0..) |shape, idx| {
            const offset = shape.getOffset();
            const z_order: u32 = shape_to_zorder[idx];
            const bin_min: u32 = if (have_bins) shape_bin_min[idx] else 0;
            const bin_max: u32 = if (have_bins) shape_bin_max[idx] else 0;
            try self.shape_offsets.append(self.allocator, offset.x);
            try self.shape_offsets.append(self.allocator, offset.y);
            try self.shape_offsets.append(self.allocator, @floatFromInt(z_order));
            try self.shape_offsets.append(self.allocator, @floatFromInt(bin_min));
            try self.shape_offsets.append(self.allocator, @floatFromInt(bin_max));
            try self.shape_offsets.append(self.allocator, 0.0);
        }
        if (tree.profile_enabled) {
            const now = nanoTimestamp();
            const prof = &@constCast(tree).last_pack_profile;
            prof.offsets_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
            t_last = now;
        }

        // Ensure at least one segment
        if (self.segment_data.items.len == 0) {
            try self.segment_data.appendNTimes(self.allocator, 0, 12);
            try self.shape_segments.append(self.allocator, 0);
            try self.shape_segments.append(self.allocator, 0);
            try self.colors.appendNTimes(self.allocator, 0, 4);
            try self.shape_bboxes.appendNTimes(self.allocator, 0, 4);
            try self.shape_offsets.appendNTimes(self.allocator, 0, 6);
        }

        // Pack slab data
        if (tree.slabs.items.len > 0) {
            if (num_bins == 0) {
                // Legacy path: one list per slab
                for (tree.slabs.items) |slab| {
                    try self.slab_xs.append(self.allocator, slab.x_left);

                    try self.slab_info.append(self.allocator, @intCast(self.active_indices.items.len)); // active_start
                    try self.slab_info.append(self.allocator, @intCast(slab.active_shapes.items.len)); // active_count
                    try self.slab_info.append(self.allocator, 0); // unused
                    try self.slab_info.append(self.allocator, 0); // padding

                    if (slab.active_shapes.items.len > self.max_active) {
                        self.max_active = @intCast(slab.active_shapes.items.len);
                    }

                    for (slab.active_shapes.items) |idx| {
                        try self.active_indices.append(self.allocator, idx);
                    }
                }
                if (tree.profile_enabled) {
                    const now = nanoTimestamp();
                    const prof = &@constCast(tree).last_pack_profile;
                    prof.slabs_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
                    t_last = now;
                }
            } else if (tree.bin_mode == .range_filter) {
                // Range-filter mode: one list per slab, filter by bin in shader.
                if (tree.profile_enabled) {
                    const prof = &@constCast(tree).last_pack_profile;
                    prof.slabs_count_ms = 0.0;
                    prof.slabs_prefix_ms = 0.0;
                    prof.slabs_global_ms = 0.0;
                    prof.slabs_bin_ms = 0.0;
                }
                for (tree.slabs.items) |slab| {
                    try self.slab_xs.append(self.allocator, slab.x_left);

                    try self.slab_info.append(self.allocator, @intCast(self.active_indices.items.len)); // active_start
                    try self.slab_info.append(self.allocator, @intCast(slab.active_shapes.items.len)); // active_count
                    try self.slab_info.append(self.allocator, 0);
                    try self.slab_info.append(self.allocator, 0);

                    if (slab.active_shapes.items.len > self.max_active) {
                        self.max_active = @intCast(slab.active_shapes.items.len);
                    }

                    for (slab.active_shapes.items) |idx| {
                        try self.active_indices.append(self.allocator, idx);
                    }
                }
                if (tree.profile_enabled) {
                    const now = nanoTimestamp();
                    const prof = &@constCast(tree).last_pack_profile;
                    prof.slabs_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
                    t_last = now;
                }
            } else if (tree.bin_mode == .premerge) {
                // Segment-tree-backed pre-merge: build y-interval segment tree per slab, then merge lists per bin.
                if (!have_bins) return error.InvalidRange;
                const slab_count: usize = tree.slabs.items.len;

                const leaf_count: u32 = blk: {
                    var p: u32 = 1;
                    while (p < num_bins) : (p <<= 1) {}
                    break :blk p;
                };
                const node_count: u32 = leaf_count * 2;
                const node_count_usize: usize = @intCast(node_count);
                const depth: u32 = blk: {
                    var d: u32 = 0;
                    var p: u32 = leaf_count;
                    while (p > 0) : (p >>= 1) {
                        d += 1;
                    }
                    break :blk d;
                };

                std.debug.assert(num_bins <= 64);
                std.debug.assert(depth <= 64);
                const path_stride: usize = @intCast(depth);
                const path_nodes_len: usize = @as(usize, @intCast(num_bins)) * path_stride;
                try self.bin_path_nodes_buf.resize(self.allocator, path_nodes_len);
                try self.bin_path_len_buf.resize(self.allocator, @intCast(num_bins));
                const bin_path_nodes = self.bin_path_nodes_buf.items;
                const bin_path_len = self.bin_path_len_buf.items;

                var pre_bin: u32 = 0;
                while (pre_bin < num_bins) : (pre_bin += 1) {
                    var n: u32 = leaf_count + pre_bin;
                    const base: usize = @as(usize, @intCast(pre_bin)) * path_stride;
                    var len: u8 = 0;
                    while (n >= 1) {
                        bin_path_nodes[base + @as(usize, len)] = n;
                        len += 1;
                        n >>= 1;
                    }
                    bin_path_len[@intCast(pre_bin)] = len;
                }

                const slab_node_counts_len: usize = slab_count * node_count_usize;
                const slab_bin_counts_len: usize = slab_count * @as(usize, @intCast(num_bins));
                try self.slab_node_counts_buf.resize(self.allocator, slab_node_counts_len);
                try self.slab_bin_counts_buf.resize(self.allocator, slab_bin_counts_len);
                try self.slab_bin_start_buf.resize(self.allocator, slab_bin_counts_len);
                try self.slab_total_active_buf.resize(self.allocator, slab_count);
                try self.slab_base_buf.resize(self.allocator, slab_count);

                const slab_node_counts = self.slab_node_counts_buf.items;
                const slab_bin_counts = self.slab_bin_counts_buf.items;
                const slab_bin_start = self.slab_bin_start_buf.items;
                const slab_total_active = self.slab_total_active_buf.items;
                const slab_base = self.slab_base_buf.items;

                const cpu_count = std.Thread.getCpuCount() catch 1;
                var threads: [64]std.Thread = undefined;
                const worker_count = @min(cpu_count, slab_count, threads.len + 1);
                const use_threads = worker_count > 1;

                var t_count_start: i128 = 0;
                var t_prefix_start: i128 = 0;
                var t_bin_start: i128 = 0;
                if (tree.profile_enabled) {
                    t_count_start = nanoTimestamp();
                    const prof = &@constCast(tree).last_pack_profile;
                    prof.slabs_count_ms = 0.0;
                    prof.slabs_prefix_ms = 0.0;
                    prof.slabs_global_ms = 0.0;
                    prof.slabs_bin_ms = 0.0;
                }

                const CountCtx = struct {
                    tree: *const StabbingTree,
                    slab_count: usize,
                    leaf_count: u32,
                    node_count_usize: usize,
                    num_bins: u32,
                    path_stride: usize,
                    bin_path_nodes: []const u32,
                    bin_path_len: []const u8,
                    shape_bin_min: []const u32,
                    shape_bin_max: []const u32,
                    slab_node_counts: []u32,
                    slab_bin_counts: []u32,
                    slab_total_active: []u32,

                    fn run(ctx: *const @This(), slab_index: usize) void {
                        var local_node_counts: [128]u32 = undefined;
                        @memset(local_node_counts[0..ctx.node_count_usize], 0);

                        const slab = ctx.tree.slabs.items[slab_index];
                        for (slab.active_shapes.items) |shape_idx| {
                            var l: u32 = ctx.leaf_count + ctx.shape_bin_min[shape_idx];
                            var r: u32 = ctx.leaf_count + ctx.shape_bin_max[shape_idx];
                            while (l <= r) {
                                if ((l & 1) == 1) {
                                    local_node_counts[@intCast(l)] += 1;
                                    l += 1;
                                }
                                if ((r & 1) == 0) {
                                    local_node_counts[@intCast(r)] += 1;
                                    if (r == 0) break;
                                    r -= 1;
                                }
                                l >>= 1;
                                r >>= 1;
                            }
                        }

                        const node_base = slab_index * ctx.node_count_usize;
                        std.mem.copyForwards(u32, ctx.slab_node_counts[node_base .. node_base + ctx.node_count_usize], local_node_counts[0..ctx.node_count_usize]);

                        const bin_base = slab_index * @as(usize, @intCast(ctx.num_bins));
                        var total_active: u32 = 0;
                        var bin_idx: u32 = 0;
                        while (bin_idx < ctx.num_bins) : (bin_idx += 1) {
                            var count: u32 = 0;
                            const base = @as(usize, @intCast(bin_idx)) * ctx.path_stride;
                            const len: usize = @intCast(ctx.bin_path_len[@intCast(bin_idx)]);
                            var i: usize = 0;
                            while (i < len) : (i += 1) {
                                count += local_node_counts[@intCast(ctx.bin_path_nodes[base + i])];
                            }
                            ctx.slab_bin_counts[bin_base + @as(usize, @intCast(bin_idx))] = count;
                            total_active += count;
                        }
                        ctx.slab_total_active[slab_index] = total_active;
                    }

                    fn worker(ctx: *const @This(), next: *std.atomic.Value(usize)) void {
                        while (true) {
                            const idx = next.fetchAdd(1, .monotonic);
                            if (idx >= ctx.slab_count) break;
                            @This().run(ctx, idx);
                        }
                    }
                };

                const count_ctx = CountCtx{
                    .tree = tree,
                    .slab_count = slab_count,
                    .leaf_count = leaf_count,
                    .node_count_usize = node_count_usize,
                    .num_bins = num_bins,
                    .path_stride = path_stride,
                    .bin_path_nodes = bin_path_nodes,
                    .bin_path_len = bin_path_len,
                    .shape_bin_min = shape_bin_min,
                    .shape_bin_max = shape_bin_max,
                    .slab_node_counts = slab_node_counts,
                    .slab_bin_counts = slab_bin_counts,
                    .slab_total_active = slab_total_active,
                };

                if (use_threads) {
                    var next = std.atomic.Value(usize).init(0);
                    var spawned: usize = 0;
                    while (spawned + 1 < worker_count) : (spawned += 1) {
                        threads[spawned] = std.Thread.spawn(.{}, CountCtx.worker, .{ &count_ctx, &next }) catch break;
                    }
                    CountCtx.worker(&count_ctx, &next);
                    for (threads[0..spawned]) |thread| thread.join();
                } else {
                    var slab_index: usize = 0;
                    while (slab_index < slab_count) : (slab_index += 1) {
                        CountCtx.run(&count_ctx, slab_index);
                    }
                }

                if (tree.profile_enabled) {
                    const now = nanoTimestamp();
                    const prof = &@constCast(tree).last_pack_profile;
                    prof.slabs_count_ms = @as(f64, @floatFromInt(now - t_count_start)) / 1_000_000.0;
                    t_prefix_start = now;
                }

                // Prefix sum for per-slab offsets (sequential)
                const old_len = self.active_indices.items.len;
                var total_active_sum: u32 = 0;
                var max_active: u32 = 0;

                var slab_index: usize = 0;
                while (slab_index < slab_count) : (slab_index += 1) {
                    const slab = tree.slabs.items[slab_index];
                    try self.slab_xs.append(self.allocator, slab.x_left);

                    const base: u32 = @as(u32, @intCast(old_len)) + total_active_sum;
                    slab_base[slab_index] = base;
                    try self.slab_info.append(self.allocator, @intCast(base));
                    try self.slab_info.append(self.allocator, 0);

                    var running: u32 = 0;
                    var bin_idx: u32 = 0;
                    const bin_base = slab_index * @as(usize, @intCast(num_bins));
                    while (bin_idx < num_bins) : (bin_idx += 1) {
                        const count = slab_bin_counts[bin_base + @as(usize, @intCast(bin_idx))];
                        const start = base + running;
                        slab_bin_start[bin_base + @as(usize, @intCast(bin_idx))] = start;
                        try self.slab_info.append(self.allocator, @intCast(start));
                        try self.slab_info.append(self.allocator, @intCast(count));
                        running += count;
                        if (count > max_active) max_active = count;
                    }

                    slab_total_active[slab_index] = running;
                    total_active_sum += running;
                }

                try self.active_indices.resize(self.allocator, old_len + total_active_sum);
                self.max_active = max_active;

                if (tree.profile_enabled) {
                    const now = nanoTimestamp();
                    const prof = &@constCast(tree).last_pack_profile;
                    prof.slabs_prefix_ms = @as(f64, @floatFromInt(now - t_prefix_start)) / 1_000_000.0;
                    t_bin_start = now;
                }

                const FillCtx = struct {
                    tree: *const StabbingTree,
                    slab_count: usize,
                    leaf_count: u32,
                    node_count: u32,
                    node_count_usize: usize,
                    num_bins: u32,
                    path_stride: usize,
                    bin_path_nodes: []const u32,
                    bin_path_len: []const u8,
                    shape_bin_min: []const u32,
                    shape_bin_max: []const u32,
                    slab_node_counts: []const u32,
                    slab_bin_counts: []const u32,
                    slab_bin_start: []const u32,
                    active_indices: []u32,
                    shape_to_zorder: []const u32,

                    fn run(ctx: *const @This(), slab_i: usize, node_lists: *std.ArrayList(u32)) void {
                        var node_offsets: [128]u32 = undefined;
                        var node_write: [128]u32 = undefined;

                        const node_base = slab_i * ctx.node_count_usize;
                        const node_counts = ctx.slab_node_counts[node_base .. node_base + ctx.node_count_usize];

                        var total: u32 = 0;
                        var node: u32 = 1;
                        while (node < ctx.node_count) : (node += 1) {
                            node_offsets[@intCast(node)] = total;
                            total += node_counts[@intCast(node)];
                        }
                        if (total == 0) return;

                        node_lists.resize(std.heap.page_allocator, total) catch {
                            @panic("node_lists resize failed");
                        };

                        const node_list_items = node_lists.items;
                        std.mem.copyForwards(u32, node_write[1..ctx.node_count_usize], node_offsets[1..ctx.node_count_usize]);

                        const slab = ctx.tree.slabs.items[slab_i];
                        for (slab.active_shapes.items) |shape_idx| {
                            var l: u32 = ctx.leaf_count + ctx.shape_bin_min[shape_idx];
                            var r: u32 = ctx.leaf_count + ctx.shape_bin_max[shape_idx];
                            while (l <= r) {
                                if ((l & 1) == 1) {
                                    const li: usize = @intCast(l);
                                    const pos = node_write[li];
                                    node_list_items[@intCast(pos)] = shape_idx;
                                    node_write[li] = pos + 1;
                                    l += 1;
                                }
                                if ((r & 1) == 0) {
                                    const ri: usize = @intCast(r);
                                    const pos = node_write[ri];
                                    node_list_items[@intCast(pos)] = shape_idx;
                                    node_write[ri] = pos + 1;
                                    if (r == 0) break;
                                    r -= 1;
                                }
                                l >>= 1;
                                r >>= 1;
                            }
                        }

                        const bin_base = slab_i * @as(usize, @intCast(ctx.num_bins));
                        var bi: u32 = 0;
                        while (bi < ctx.num_bins) : (bi += 1) {
                            const count = ctx.slab_bin_counts[bin_base + @as(usize, @intCast(bi))];
                            if (count == 0) continue;
                            const start = ctx.slab_bin_start[bin_base + @as(usize, @intCast(bi))];

                            var path_pos: [64]u32 = undefined;
                            var path_end: [64]u32 = undefined;
                            const path_base: usize = @as(usize, @intCast(bi)) * ctx.path_stride;
                            const path_len: usize = @as(usize, @intCast(ctx.bin_path_len[@as(usize, @intCast(bi))]));
                            var i: usize = 0;
                            while (i < path_len) : (i += 1) {
                                const node_idx = ctx.bin_path_nodes[path_base + i];
                                const off = node_offsets[@intCast(node_idx)];
                                const cnt = node_counts[@intCast(node_idx)];
                                path_pos[i] = off;
                                path_end[i] = off + cnt;
                            }

                            const HeapItem = struct {
                                z: u32,
                                path: u8,
                                shape: u32,
                            };
                            var heap: [64]HeapItem = undefined;
                            var heap_len: usize = 0;

                            i = 0;
                            while (i < path_len) : (i += 1) {
                                const pos = path_pos[i];
                                if (pos >= path_end[i]) continue;
                                const shape_idx = node_list_items[@intCast(pos)];
                                const item = HeapItem{
                                    .z = ctx.shape_to_zorder[shape_idx],
                                    .path = @intCast(i),
                                    .shape = shape_idx,
                                };
                                var idx = heap_len;
                                heap_len += 1;
                                heap[idx] = item;
                                while (idx > 0) {
                                    const parent = (idx - 1) / 2;
                                    const a = heap[parent];
                                    const b = heap[idx];
                                    const less = (b.z < a.z) or (b.z == a.z and b.shape < a.shape);
                                    if (!less) break;
                                    heap[parent] = b;
                                    heap[idx] = a;
                                    idx = parent;
                                }
                            }

                            var out_i: u32 = 0;
                            while (out_i < count and heap_len > 0) : (out_i += 1) {
                                const best = heap[0];
                                ctx.active_indices[@intCast(start + out_i)] = best.shape;

                                const p: usize = @intCast(best.path);
                                path_pos[p] += 1;

                                if (heap_len == 1) {
                                    heap_len = 0;
                                } else {
                                    heap_len -= 1;
                                    heap[0] = heap[heap_len];
                                    var idx: usize = 0;
                                    while (true) {
                                        const left = idx * 2 + 1;
                                        if (left >= heap_len) break;
                                        var smallest = left;
                                        const right = left + 1;
                                        if (right < heap_len) {
                                            const a = heap[smallest];
                                            const b = heap[right];
                                            const less = (b.z < a.z) or (b.z == a.z and b.shape < a.shape);
                                            if (less) smallest = right;
                                        }
                                        const a = heap[idx];
                                        const b = heap[smallest];
                                        const less = (b.z < a.z) or (b.z == a.z and b.shape < a.shape);
                                        if (!less) break;
                                        heap[idx] = b;
                                        heap[smallest] = a;
                                        idx = smallest;
                                    }
                                }

                                if (path_pos[p] < path_end[p]) {
                                    const shape_idx = node_list_items[@intCast(path_pos[p])];
                                    const item = HeapItem{
                                        .z = ctx.shape_to_zorder[shape_idx],
                                        .path = @intCast(p),
                                        .shape = shape_idx,
                                    };
                                    var idx = heap_len;
                                    heap_len += 1;
                                    heap[idx] = item;
                                    while (idx > 0) {
                                        const parent = (idx - 1) / 2;
                                        const a2 = heap[parent];
                                        const b2 = heap[idx];
                                        const less2 = (b2.z < a2.z) or (b2.z == a2.z and b2.shape < a2.shape);
                                        if (!less2) break;
                                        heap[parent] = b2;
                                        heap[idx] = a2;
                                        idx = parent;
                                    }
                                }
                            }
                        }
                    }

                    fn worker(ctx: *const @This(), next: *std.atomic.Value(usize)) void {
                        var node_lists = std.ArrayList(u32).empty;
                        defer node_lists.deinit(std.heap.page_allocator);
                        while (true) {
                            const idx = next.fetchAdd(1, .monotonic);
                            if (idx >= ctx.slab_count) break;
                            @This().run(ctx, idx, &node_lists);
                        }
                    }
                };

                const fill_ctx = FillCtx{
                    .tree = tree,
                    .slab_count = slab_count,
                    .leaf_count = leaf_count,
                    .node_count = node_count,
                    .node_count_usize = node_count_usize,
                    .num_bins = num_bins,
                    .path_stride = path_stride,
                    .bin_path_nodes = bin_path_nodes,
                    .bin_path_len = bin_path_len,
                    .shape_bin_min = shape_bin_min,
                    .shape_bin_max = shape_bin_max,
                    .slab_node_counts = slab_node_counts,
                    .slab_bin_counts = slab_bin_counts,
                    .slab_bin_start = slab_bin_start,
                    .active_indices = self.active_indices.items,
                    .shape_to_zorder = shape_to_zorder,
                };

                if (use_threads) {
                    var next = std.atomic.Value(usize).init(0);
                    var spawned: usize = 0;
                    while (spawned + 1 < worker_count) : (spawned += 1) {
                        threads[spawned] = std.Thread.spawn(.{}, FillCtx.worker, .{ &fill_ctx, &next }) catch break;
                    }
                    FillCtx.worker(&fill_ctx, &next);
                    for (threads[0..spawned]) |thread| thread.join();
                } else {
                    var slab_index2: usize = 0;
                    var node_lists = std.ArrayList(u32).empty;
                    defer node_lists.deinit(std.heap.page_allocator);
                    while (slab_index2 < slab_count) : (slab_index2 += 1) {
                        FillCtx.run(&fill_ctx, slab_index2, &node_lists);
                    }
                }

                if (tree.profile_enabled) {
                    const now = nanoTimestamp();
                    const prof = &@constCast(tree).last_pack_profile;
                    prof.slabs_bin_ms = @as(f64, @floatFromInt(now - t_bin_start)) / 1_000_000.0;
                    prof.slabs_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
                    t_last = now;
                }
            } else {
                // Bin lists path: global + per-bin indices
                if (!have_bins) return error.InvalidRange;

                var bin_counts: []u32 = undefined;
                var bin_offsets: []u32 = undefined;
                var free_bins = false;
                if (tree.pack_reuse_bins) {
                    const bins_len: usize = @intCast(num_bins);
                    try self.bin_counts_buf.resize(self.allocator, bins_len);
                    try self.bin_offsets_buf.resize(self.allocator, bins_len);
                    bin_counts = self.bin_counts_buf.items;
                    bin_offsets = self.bin_offsets_buf.items;
                } else {
                    bin_counts = try self.allocator.alloc(u32, num_bins);
                    bin_offsets = try self.allocator.alloc(u32, num_bins);
                    free_bins = true;
                }
                defer if (free_bins) {
                    self.allocator.free(bin_counts);
                    self.allocator.free(bin_offsets);
                };

                var t_count_start: i128 = 0;
                var t_prefix_start: i128 = 0;
                var t_global_start: i128 = 0;
                var t_bin_start: i128 = 0;
                if (tree.profile_enabled) {
                    t_count_start = nanoTimestamp();
                    const prof = &@constCast(tree).last_pack_profile;
                    prof.slabs_count_ms = 0.0;
                    prof.slabs_prefix_ms = 0.0;
                    prof.slabs_global_ms = 0.0;
                    prof.slabs_bin_ms = 0.0;
                }

                for (tree.slabs.items) |slab| {
                    try self.slab_xs.append(self.allocator, slab.x_left);

                    @memset(bin_counts, 0);

                    var global_count: u32 = 0;
                    const wide_threshold: u32 = tree.wide_threshold;

                    // Count shapes per bin
                    for (slab.active_shapes.items) |shape_idx| {
                        const bin_min = shape_bin_min[shape_idx];
                        const bin_max = shape_bin_max[shape_idx];
                        const span = bin_max - bin_min + 1;
                        if (span > wide_threshold) {
                            global_count += 1;
                        } else {
                            var b: u32 = bin_min;
                            while (b <= bin_max) : (b += 1) {
                                bin_counts[b] += 1;
                            }
                        }
                    }

                    if (tree.profile_enabled) {
                        const now = nanoTimestamp();
                        const prof = &@constCast(tree).last_pack_profile;
                        prof.slabs_count_ms += @as(f64, @floatFromInt(now - t_count_start)) / 1_000_000.0;
                        t_prefix_start = now;
                    }

                    // Prefix sums to compute starts
                    var total: u32 = global_count;
                    var b: usize = 0;
                    while (b < num_bins) : (b += 1) {
                        bin_offsets[b] = total;
                        total += bin_counts[b];
                    }

                    if (tree.profile_enabled) {
                        const now = nanoTimestamp();
                        const prof = &@constCast(tree).last_pack_profile;
                        prof.slabs_prefix_ms += @as(f64, @floatFromInt(now - t_prefix_start)) / 1_000_000.0;
                        t_global_start = now;
                    }

                    const base = @as(u32, @intCast(self.active_indices.items.len));
                    const old_len = self.active_indices.items.len;
                    try self.active_indices.resize(self.allocator, old_len + total);

                    // Write slab info (global + bins)
                    try self.slab_info.append(self.allocator, @intCast(base)); // global_start
                    try self.slab_info.append(self.allocator, @intCast(global_count)); // global_count

                    if (global_count > self.max_active) {
                        self.max_active = global_count;
                    }

                    b = 0;
                    while (b < num_bins) : (b += 1) {
                        const start = base + bin_offsets[b];
                        const count = bin_counts[b];
                        try self.slab_info.append(self.allocator, @intCast(start));
                        try self.slab_info.append(self.allocator, @intCast(count));
                        // convert offsets to absolute write positions
                        bin_offsets[b] = start;

                        if (count > self.max_active) {
                            self.max_active = count;
                        }
                    }

                    // Fill global list
                    var global_write = base;
                    for (slab.active_shapes.items) |shape_idx| {
                        const bin_min = shape_bin_min[shape_idx];
                        const bin_max = shape_bin_max[shape_idx];
                        const span = bin_max - bin_min + 1;
                        if (span > wide_threshold) {
                            self.active_indices.items[@intCast(global_write)] = shape_idx;
                            global_write += 1;
                        }
                    }

                    if (tree.profile_enabled) {
                        const now = nanoTimestamp();
                        const prof = &@constCast(tree).last_pack_profile;
                        prof.slabs_global_ms += @as(f64, @floatFromInt(now - t_global_start)) / 1_000_000.0;
                        t_bin_start = now;
                    }

                    // Fill bin indices
                    for (slab.active_shapes.items) |shape_idx| {
                        const bin_min = shape_bin_min[shape_idx];
                        const bin_max = shape_bin_max[shape_idx];
                        const span = bin_max - bin_min + 1;
                        if (span > wide_threshold) continue;
                        var b2: u32 = bin_min;
                        while (b2 <= bin_max) : (b2 += 1) {
                            const bi: usize = @intCast(b2);
                            const pos = bin_offsets[bi];
                            self.active_indices.items[@intCast(pos)] = shape_idx;
                            bin_offsets[bi] = pos + 1;
                        }
                    }

                    if (tree.profile_enabled) {
                        const now = nanoTimestamp();
                        const prof = &@constCast(tree).last_pack_profile;
                        prof.slabs_bin_ms += @as(f64, @floatFromInt(now - t_bin_start)) / 1_000_000.0;
                        t_count_start = now;
                    }
                }
                if (tree.profile_enabled) {
                    const now = nanoTimestamp();
                    const prof = &@constCast(tree).last_pack_profile;
                    prof.slabs_ms = @as(f64, @floatFromInt(now - t_last)) / 1_000_000.0;
                    t_last = now;
                }
            }
        }
        if (tree.profile_enabled) {
            const now = nanoTimestamp();
            const prof = &@constCast(tree).last_pack_profile;
            prof.total_ms = @as(f64, @floatFromInt(now - t_start)) / 1_000_000.0;
        }

        // Ensure at least one of each
        if (self.slab_xs.items.len == 0) {
            try self.slab_xs.append(self.allocator, 0);
            try self.slab_info.appendSlice(self.allocator, &[_]i32{ 0, 0, 0, 0 });
        }
        if (self.active_indices.items.len == 0) {
            try self.active_indices.append(self.allocator, 0);
        }
    }
};

// Tests
test "stabbing tree basic" {
    const allocator = std.testing.allocator;
    const shapes_mod = @import("../geometry/shapes.zig");

    // Create some shapes
    var shape1 = shapes_mod.ConvexBezierShape.init(allocator);
    defer shape1.deinit();

    try shape1.addCurve(shapes_mod.CubicBezier.init(
        shapes_mod.Vec2.init(0, 0),
        shapes_mod.Vec2.init(0.5, 0),
        shapes_mod.Vec2.init(1, 0.5),
        shapes_mod.Vec2.init(1, 1),
    ));
    try shape1.addCurve(shapes_mod.CubicBezier.init(
        shapes_mod.Vec2.init(1, 1),
        shapes_mod.Vec2.init(1, 1.5),
        shapes_mod.Vec2.init(0.5, 2),
        shapes_mod.Vec2.init(0, 2),
    ));
    try shape1.addCurve(shapes_mod.CubicBezier.init(
        shapes_mod.Vec2.init(0, 2),
        shapes_mod.Vec2.init(-0.5, 2),
        shapes_mod.Vec2.init(-1, 1.5),
        shapes_mod.Vec2.init(-1, 1),
    ));
    try shape1.addCurve(shapes_mod.CubicBezier.init(
        shapes_mod.Vec2.init(-1, 1),
        shapes_mod.Vec2.init(-1, 0.5),
        shapes_mod.Vec2.init(-0.5, 0),
        shapes_mod.Vec2.init(0, 0),
    ));
    try shape1.precompute();

    var stabbables = [_]Stabbable{shape1.asStabbable()};

    var tree_instance = StabbingTree.init(allocator);
    defer tree_instance.deinit();

    try tree_instance.build(&stabbables);

    // Query point inside shape
    var results: std.ArrayList(u32) = .empty;
    defer results.deinit(allocator);

    try tree_instance.query(0, 1, &results);
    try std.testing.expect(results.items.len > 0);
}
