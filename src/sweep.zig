//! The persistent sweep behind the trees in this package.
//!
//! Shapes are described by the x-monotone arcs of their boundaries. A vertical
//! line sweeps the plane left to right, and the arcs it cuts are kept in a
//! treap ordered by y. Every arc carries the set of shapes covering the face
//! directly above it: crossing an arc upwards toggles that arc's shape, so a
//! point is inside a shape when an odd number of its arcs lie below.
//!
//! The order only changes where arcs start, end, or cross, so there is one
//! version of the treap per such x. Both the treap and the sets are persistent,
//! which makes each new version cost O(log n) instead of a full copy: treap
//! nodes hold one timestamped override before they are copied, and the sets
//! are path-copied.
//!
//! A query binary searches the versions by x, walks one version by y, and
//! reads the set it lands on.

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const nil: u32 = std.math.maxInt(u32);

/// Bijective, so distinct keys get distinct treap priorities.
fn mix(key: u32) u32 {
    var h = key;
    h ^= h >> 16;
    h *%= 0x85ebca6b;
    h ^= h >> 13;
    h *%= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}

/// What a geometry tells the sweep. Arcs are identified by index.
pub const Event = struct {
    x: f64,
    a: u32,
    b: u30,
    kind: Kind,

    pub const Kind = enum(u2) {
        /// Arcs `a` and `b` begin, next to each other.
        start,
        /// Arcs `a` and `b` end, next to each other.
        end,
        /// Arcs `a` and `b` may swap places near `x`. This is only a hint to
        /// re-examine them: the sweep decides the order by evaluating the arcs,
        /// so extra hints are harmless and a hint may be off by `slack`.
        cross,
    };

    pub const max_arcs = std.math.maxInt(u30);

    fn before(_: void, l: Event, r: Event) bool {
        return l.x < r.x;
    }
};

/// An axis-aligned bounding box, for `forEachOverlap`.
pub const Box = struct {
    x: [2]f64,
    y: [2]f64,
    id: u32,

    fn before(_: void, l: Box, r: Box) bool {
        return l.x[0] < r.x[0];
    }
};

/// Calls `context.pair(a, b)` with the ids of every two boxes that overlap,
/// once per pair. Reorders `boxes`.
///
/// The boxes are swept in x, and the ones the sweep line is inside are kept in
/// rows by y so that a box is only compared with those near it.
pub fn forEachOverlap(gpa: Allocator, boxes: []Box, context: anytype) Allocator.Error!void {
    if (boxes.len < 2) return;
    std.mem.sortUnstable(Box, boxes, {}, Box.before);

    var bottom = std.math.inf(f64);
    var top = -std.math.inf(f64);
    var heights: f64 = 0;
    for (boxes) |box| {
        bottom = @min(bottom, box.y[0]);
        top = @max(top, box.y[1]);
        heights += box.y[1] - box.y[0];
    }
    // Rows about as tall as the average box.
    const wanted = (top - bottom) * @as(f64, @floatFromInt(boxes.len)) / heights;
    const row_count: usize = if (wanted >= 1) @intFromFloat(@min(wanted, @as(f64, @floatFromInt(boxes.len)))) else 1;
    const per_y = @as(f64, @floatFromInt(row_count)) / (top - bottom);

    const rows = try gpa.alloc(std.ArrayList(Box), row_count);
    defer gpa.free(rows);
    @memset(rows, .empty);
    defer for (rows) |*row| row.deinit(gpa);

    const Rows = struct {
        bottom: f64,
        per_y: f64,
        last: usize,

        fn of(self: @This(), y: f64) usize {
            const row = (y - self.bottom) * self.per_y;
            return if (row >= 0) @min(self.last, @as(usize, @intFromFloat(row))) else 0;
        }
    };
    const grid: Rows = .{ .bottom = bottom, .per_y = if (std.math.isFinite(per_y)) per_y else 0, .last = row_count - 1 };

    for (boxes) |box| {
        const first = grid.of(box.y[0]);
        const last = grid.of(box.y[1]);
        for (rows[first .. last + 1], first..) |*row, index| {
            var k: usize = 0;
            while (k < row.items.len) {
                const other = row.items[k];
                if (other.x[1] < box.x[0]) {
                    _ = row.swapRemove(k);
                    continue;
                }
                // Boxes sharing several rows are paired in the first of them.
                if (other.y[0] <= box.y[1] and box.y[0] <= other.y[1] and
                    @max(first, grid.of(other.y[0])) == index)
                {
                    try context.pair(other.id, box.id);
                }
                k += 1;
            }
            try row.append(gpa, box);
        }
    }
}

/// Node of a persistent treap over the arcs cut by the sweep line, in y order.
/// A node absorbs one change in place; the next one copies it, which changes
/// its parent in turn.
const SeqNode = struct {
    first: State,
    /// Version from which `second` replaces `first`, or `nil`.
    since: u32 = nil,
    second: State = undefined,

    const State = struct {
        left: u32,
        right: u32,
        arc: u32,
        /// Shapes covering the face directly above this arc.
        set: u32,
    };

    fn at(n: *const SeqNode, version: u32) State {
        return if (n.since <= version) n.second else n.first;
    }
};

/// A stabbing tree over the shapes of a `Geometry`, carrying a `T` per shape.
///
/// `combine(a, b)` must be associative. It is always called with the shapes
/// behind `a` at lower indices than those behind `b`, so it need not be
/// commutative.
///
/// `Geometry` provides:
///
///     const Input: type
///     fn init(gpa, input: Input, events: *std.ArrayList(Event)) Allocator.Error!Geometry
///     fn deinit(*Geometry, gpa) void
///     fn shapeCount(*const Geometry) usize
///     fn arcCount(*const Geometry) usize
///     fn shapeOf(*const Geometry, arc: u32) u32
///     fn arcY(*const Geometry, arc: u32, x: f64) f64
///     fn slack(*const Geometry, arc: u32) f64
///     fn bytes(*const Geometry) usize
///
/// `init` appends a start and an end event for every arc, in pairs, and a
/// cross event wherever two arcs may cross. `arcY` is the arc's height at x,
/// with x clamped to the arc's extent. `slack` is how far a cross event
/// involving the arc may be from the true crossing.
///
/// A geometry whose `arcY` is expensive can also provide
///
///     fn arcBounds(*const Geometry, arc: u32, x: f64) [2]f64
///
/// with cheap lower and upper bounds on `arcY(arc, x)`, which settle most
/// comparisons without evaluating the arc.
pub fn Tree(comptime Geometry: type, comptime T: type, comptime combine: fn (T, T) T) type {
    return struct {
        const Self = @This();

        /// Node of a persistent treap of shape ids.
        const SetNode = struct {
            left: u32,
            right: u32,
            shape: u32,
            /// The values of this subtree's shapes, combined in index order.
            folded: T,
        };

        const Values = if (T == void) void else []const T;

        geometry: Geometry,
        /// Left edge of each slab, ascending.
        xs: []const f64,
        /// Root in `seq` of the arc order inside each slab. The slab's index is
        /// the version to read the nodes at.
        roots: []const u32,
        seq: []const SeqNode,
        sets: []const SetNode,

        /// Time and space are O((n + c) log n) for n arcs and c crossings
        /// between them, after the geometry has found the crossings.
        pub const init = if (T == void) initPlain else initFolded;

        fn initPlain(gpa: Allocator, input: Geometry.Input) Allocator.Error!Self {
            return build(gpa, input, {});
        }

        /// `values[i]` belongs to shape i. The values are only read here.
        fn initFolded(gpa: Allocator, input: Geometry.Input, values: []const T) Allocator.Error!Self {
            return build(gpa, input, values);
        }

        fn build(gpa: Allocator, input: Geometry.Input, values: Values) Allocator.Error!Self {
            var events: std.ArrayList(Event) = .empty;
            defer events.deinit(gpa);
            var geometry = try Geometry.init(gpa, input, &events);
            errdefer geometry.deinit(gpa);
            std.debug.assert(geometry.arcCount() <= Event.max_arcs);
            if (T != void) std.debug.assert(values.len == geometry.shapeCount());
            std.mem.sortUnstable(Event, events.items, {}, Event.before);

            var sweep = try Sweep.init(gpa, &geometry, values);
            defer sweep.deinit();
            var xs: std.ArrayList(f64) = .empty;
            defer xs.deinit(gpa);
            var roots: std.ArrayList(u32) = .empty;
            defer roots.deinit(gpa);

            var i: usize = 0;
            while (i < events.items.len) {
                const x = events.items[i].x;
                var j = i + 1;
                while (j < events.items.len and events.items[j].x == x) j += 1;
                const x_eval = if (j < events.items.len) x + (events.items[j].x - x) / 2 else x;

                if (try sweep.step(events.items[i..j], x, x_eval)) |root| {
                    try xs.append(gpa, x);
                    try roots.append(gpa, root);
                }
                i = j;
            }

            const xs_owned = try xs.toOwnedSlice(gpa);
            errdefer gpa.free(xs_owned);
            const roots_owned = try roots.toOwnedSlice(gpa);
            errdefer gpa.free(roots_owned);
            const seq_owned = try sweep.seq.toOwnedSlice(gpa);
            errdefer gpa.free(seq_owned);
            const sets_owned = try sweep.sets.toOwnedSlice(gpa);

            return .{
                .geometry = geometry,
                .xs = xs_owned,
                .roots = roots_owned,
                .seq = seq_owned,
                .sets = sets_owned,
            };
        }

        pub fn deinit(self: *Self, gpa: Allocator) void {
            self.geometry.deinit(gpa);
            gpa.free(self.xs);
            gpa.free(self.roots);
            gpa.free(self.seq);
            gpa.free(self.sets);
            self.* = undefined;
        }

        /// Ids of the shapes containing (x, y), ascending. Caller frees.
        ///
        /// A point within rounding error of a shape's boundary may be reported
        /// either way for that shape.
        pub fn query(self: *const Self, gpa: Allocator, x: f64, y: f64) Allocator.Error![]u32 {
            var out: std.ArrayList(u32) = .empty;
            errdefer out.deinit(gpa);
            try self.queryAppend(gpa, &out, x, y);
            return out.toOwnedSlice(gpa);
        }

        /// Like `query`, but appends to `out` so a buffer can be reused.
        pub fn queryAppend(
            self: *const Self,
            gpa: Allocator,
            out: *std.ArrayList(u32),
            x: f64,
            y: f64,
        ) Allocator.Error!void {
            try self.appendSet(gpa, out, self.locate(x, y));
        }

        /// The values of the shapes containing (x, y), combined in index order, or
        /// null if there are none. O(log n) however many shapes that is.
        pub fn fold(self: *const Self, x: f64, y: f64) ?T {
            const set = self.locate(x, y);
            return if (set == nil) null else self.sets[set].folded;
        }

        /// Whether the arc is at or below y at x.
        fn arcBelow(geometry: *const Geometry, arc: u32, x: f64, y: f64) bool {
            if (@hasDecl(Geometry, "arcBounds")) {
                const bounds = geometry.arcBounds(arc, x);
                if (bounds[1] <= y) return true;
                if (bounds[0] > y) return false;
            }
            return geometry.arcY(arc, x) <= y;
        }

        /// The set of shapes covering the face that (x, y) lies in.
        fn locate(self: *const Self, x: f64, y: f64) u32 {
            var lo: usize = 0;
            var hi: usize = self.xs.len;
            while (lo < hi) {
                const mid = lo + (hi - lo) / 2;
                if (self.xs[mid] <= x) lo = mid + 1 else hi = mid;
            }
            if (lo == 0) return nil;

            const version: u32 = @intCast(lo - 1);
            var node = self.roots[version];
            var set = nil;
            while (node != nil) {
                const n = self.seq[node].at(version);
                if (arcBelow(&self.geometry, n.arc, x, y)) {
                    set = n.set;
                    node = n.right;
                } else {
                    node = n.left;
                }
            }
            return set;
        }

        fn appendSet(self: *const Self, gpa: Allocator, out: *std.ArrayList(u32), set: u32) Allocator.Error!void {
            var t = set;
            while (t != nil) {
                const n = self.sets[t];
                try self.appendSet(gpa, out, n.left);
                try out.append(gpa, n.shape);
                t = n.right;
            }
        }

        /// Heap bytes held by the structure.
        pub fn bytes(self: *const Self) usize {
            return self.geometry.bytes() + std.mem.sliceAsBytes(self.xs).len +
                std.mem.sliceAsBytes(self.roots).len + std.mem.sliceAsBytes(self.seq).len +
                std.mem.sliceAsBytes(self.sets).len;
        }

        /// The sweep line. `slots` is a mutable treap with the same shape as the
        /// newest persistent version, plus parent and neighbour links so an arc can be
        /// found and moved without searching. An arc's slot is created under the
        /// arc's own index when it starts; crossings exchange the arcs held by two
        /// slots rather than moving the slots.
        const Sweep = struct {
            gpa: Allocator,
            geometry: *const Geometry,
            values: Values,
            seq: std.ArrayList(SeqNode) = .empty,
            sets: std.ArrayList(SetNode) = .empty,
            slots: []Slot,
            arcs: []Arc,
            root: u32 = nil,
            head: u32 = nil,
            batch: u32 = 0,
            /// Index of the slab being built.
            version: u32 = 0,
            last_root: u32 = nil,
            changed: bool = false,
            /// Arcs are compared at this x, strictly inside the slab being built.
            x_eval: f64 = 0,
            /// Slots whose persistent node is out of date.
            touched: std.ArrayList(u32) = .empty,
            /// Arcs that may be out of order in this slab.
            dirty: std.ArrayList(u32) = .empty,
            /// Arcs to re-examine in the next slab as well.
            carry: std.ArrayList(u32) = .empty,

            const Slot = struct {
                left: u32 = nil,
                right: u32 = nil,
                parent: u32 = nil,
                prev: u32 = nil,
                next: u32 = nil,
                arc: u32 = nil,
                set: u32 = nil,
                /// This slot's node in `seq` as of the newest version.
                node: u32 = nil,
                mark: u32 = 0,
            };

            const Arc = struct {
                slot: u32 = nil,
                stamp: u32 = 0,
                expiry: f64 = -std.math.inf(f64),
                height_batch: u32 = 0,
                height: f64 = 0,
            };

            fn init(gpa: Allocator, geometry: *const Geometry, values: Values) Allocator.Error!Sweep {
                const slots = try gpa.alloc(Slot, geometry.arcCount());
                errdefer gpa.free(slots);
                const arcs = try gpa.alloc(Arc, geometry.arcCount());
                @memset(slots, .{});
                @memset(arcs, .{});
                return .{ .gpa = gpa, .geometry = geometry, .values = values, .slots = slots, .arcs = arcs };
            }

            fn deinit(self: *Sweep) void {
                self.seq.deinit(self.gpa);
                self.sets.deinit(self.gpa);
                self.touched.deinit(self.gpa);
                self.dirty.deinit(self.gpa);
                self.carry.deinit(self.gpa);
                self.gpa.free(self.slots);
                self.gpa.free(self.arcs);
            }

            /// Applies all events at `x`. Returns the new version's root, or null if
            /// nothing changed.
            fn step(self: *Sweep, events: []const Event, x: f64, x_eval: f64) Allocator.Error!?u32 {
                self.batch += 1;
                self.x_eval = x_eval;

                self.dirty.clearRetainingCapacity();
                for (self.carry.items) |arc| try self.markDirty(arc, x);

                for (events) |e| if (e.kind == .end) try self.endPair(e.a, e.b);
                for (events) |e| if (e.kind == .start) try self.startPair(e.a, e.b, x);
                for (events) |e| if (e.kind == .cross) {
                    try self.markDirty(e.a, x);
                    try self.markDirty(e.b, x);
                };
                try self.settle();

                self.carry.clearRetainingCapacity();
                for (self.dirty.items) |arc| {
                    const state = self.arcs[arc];
                    if (state.slot != nil and state.expiry > x_eval) try self.carry.append(self.gpa, arc);
                }

                for (self.touched.items) |slot| {
                    if (!self.alive(slot)) continue;
                    var t = slot;
                    while (t != nil and self.slots[t].mark != self.batch) : (t = self.slots[t].parent) {
                        self.slots[t].mark = self.batch;
                    }
                }
                self.touched.clearRetainingCapacity();

                const root = try self.commit(self.root);
                if (!self.changed and root == self.last_root) return null;
                self.changed = false;
                self.last_root = root;
                self.version += 1;
                return root;
            }

            /// Brings the persistent nodes of the marked slots, which are closed under
            /// taking parents, up to date. Returns the node now standing for `slot`.
            fn commit(self: *Sweep, slot: u32) Allocator.Error!u32 {
                if (slot == nil) return nil;
                const s = self.slots[slot];
                if (s.mark != self.batch) return s.node;
                self.slots[slot].mark = 0;
                const want: SeqNode.State = .{
                    .left = try self.commit(s.left),
                    .right = try self.commit(s.right),
                    .arc = s.arc,
                    .set = s.set,
                };

                if (s.node != nil) {
                    const n = &self.seq.items[s.node];
                    if (std.meta.eql(n.at(self.version), want)) return s.node;
                    if (n.since == nil) {
                        n.since = self.version;
                        n.second = want;
                        self.changed = true;
                        return s.node;
                    }
                }
                const node: u32 = @intCast(self.seq.items.len);
                try self.seq.append(self.gpa, .{ .first = want });
                self.slots[slot].node = node;
                self.changed = true;
                return node;
            }

            fn alive(self: *const Sweep, slot: u32) bool {
                const arc = self.slots[slot].arc;
                return arc != nil and self.arcs[arc].slot == slot;
            }

            fn less(self: *Sweep, a: u32, b: u32) bool {
                if (@hasDecl(Geometry, "arcBounds")) {
                    const ba = self.geometry.arcBounds(a, self.x_eval);
                    const bb = self.geometry.arcBounds(b, self.x_eval);
                    if (ba[1] < bb[0]) return true;
                    if (bb[1] < ba[0]) return false;
                }
                const ya = self.height(a);
                const yb = self.height(b);
                return ya < yb or (ya == yb and a < b);
            }

            /// The arc's height at `x_eval`, computed once per batch.
            fn height(self: *Sweep, arc: u32) f64 {
                const state = &self.arcs[arc];
                if (state.height_batch != self.batch) {
                    state.height_batch = self.batch;
                    state.height = self.geometry.arcY(arc, self.x_eval);
                }
                return state.height;
            }

            /// Crossing coordinates carry rounding error, so an arc touched at `x`
            /// stays dirty until the sweep is clear of where the event may really be.
            fn markDirty(self: *Sweep, arc: u32, x: f64) Allocator.Error!void {
                const state = &self.arcs[arc];
                if (state.slot == nil) return;
                state.expiry = @max(state.expiry, x + self.geometry.slack(arc));
                if (state.stamp == self.batch) return;
                state.stamp = self.batch;
                try self.dirty.append(self.gpa, arc);
            }

            /// Restores y order at `x_eval` by moving dirty arcs past their
            /// neighbours. Arcs that are not dirty are already in order relative to
            /// each other, so this sorts the whole sequence.
            fn settle(self: *Sweep) Allocator.Error!void {
                var changed = true;
                while (changed) {
                    changed = false;
                    for (self.dirty.items) |arc| {
                        var slot = self.arcs[arc].slot;
                        if (slot == nil) continue;
                        while (true) {
                            const prev = self.slots[slot].prev;
                            if (prev == nil or !self.less(arc, self.slots[prev].arc)) break;
                            try self.swap(prev);
                            slot = prev;
                            changed = true;
                        }
                        while (true) {
                            const next = self.slots[slot].next;
                            if (next == nil or !self.less(self.slots[next].arc, arc)) break;
                            try self.swap(slot);
                            slot = next;
                            changed = true;
                        }
                    }
                }
            }

            /// Exchanges the arc in `lower` with the one in the slot above it. Only the
            /// face between the two arcs changes.
            fn swap(self: *Sweep, lower: u32) Allocator.Error!void {
                const upper = self.slots[lower].next;
                const down = self.slots[upper].arc;
                const up = self.slots[lower].arc;
                const prev = self.slots[lower].prev;
                const below = if (prev == nil) nil else self.slots[prev].set;

                self.slots[lower].arc = down;
                self.slots[upper].arc = up;
                self.arcs[down].slot = lower;
                self.arcs[up].slot = upper;
                self.slots[lower].set = try self.setToggle(below, self.geometry.shapeOf(down));

                try self.touched.append(self.gpa, lower);
                try self.touched.append(self.gpa, upper);
            }

            /// Two arcs appear at the same x. Which of them is lower is left
            /// to `settle`.
            fn startPair(self: *Sweep, a: u32, b: u32, x: f64) Allocator.Error!void {
                var below: u32 = nil;
                var t = self.root;
                while (t != nil) {
                    if (self.less(self.slots[t].arc, a)) {
                        below = t;
                        t = self.slots[t].right;
                    } else {
                        t = self.slots[t].left;
                    }
                }
                const outside = if (below == nil) nil else self.slots[below].set;
                const between = try self.setToggle(outside, self.geometry.shapeOf(a));

                self.slots[a] = .{ .arc = a, .set = between };
                self.slots[b] = .{ .arc = b, .set = try self.setToggle(between, self.geometry.shapeOf(b)) };
                self.arcs[a].slot = a;
                self.arcs[b].slot = b;
                try self.insertAfter(below, a);
                try self.insertAfter(a, b);
                try self.markDirty(a, x);
                try self.markDirty(b, x);
            }

            /// Two arcs disappear at the same x.
            fn endPair(self: *Sweep, a: u32, b: u32) Allocator.Error!void {
                var lo = a;
                var hi = b;
                var t = self.arcs[b].slot;
                while (t != nil and t != self.arcs[a].slot) t = self.slots[t].prev;
                if (t == nil) {
                    lo = b;
                    hi = a;
                }
                // Anything still between the two arcs is passed here.
                while (self.slots[self.arcs[hi].slot].prev != self.arcs[lo].slot) {
                    try self.swap(self.slots[self.arcs[hi].slot].prev);
                }
                try self.remove(self.arcs[hi].slot);
                try self.remove(self.arcs[lo].slot);
            }

            fn insertAfter(self: *Sweep, before: u32, slot: u32) Allocator.Error!void {
                const slots = self.slots;
                const after = if (before == nil) self.head else slots[before].next;
                slots[slot].prev = before;
                slots[slot].next = after;
                if (before == nil) self.head = slot else slots[before].next = slot;
                if (after != nil) slots[after].prev = slot;

                // Either `before` has no right child or `after` has no left child.
                if (before != nil and slots[before].right == nil) {
                    slots[before].right = slot;
                    slots[slot].parent = before;
                    try self.touched.append(self.gpa, before);
                } else if (after != nil) {
                    slots[after].left = slot;
                    slots[slot].parent = after;
                    try self.touched.append(self.gpa, after);
                } else {
                    self.root = slot;
                }
                try self.touched.append(self.gpa, slot);

                while (slots[slot].parent != nil and mix(slot) > mix(slots[slot].parent)) {
                    try self.rotateUp(slot);
                }
            }

            fn remove(self: *Sweep, slot: u32) Allocator.Error!void {
                const slots = self.slots;
                while (slots[slot].left != nil and slots[slot].right != nil) {
                    const l = slots[slot].left;
                    const r = slots[slot].right;
                    try self.rotateUp(if (mix(l) > mix(r)) l else r);
                }

                const s = slots[slot];
                const child = if (s.left != nil) s.left else s.right;
                if (child != nil) slots[child].parent = s.parent;
                if (s.parent == nil) {
                    self.root = child;
                } else {
                    if (slots[s.parent].left == slot) slots[s.parent].left = child else slots[s.parent].right = child;
                    try self.touched.append(self.gpa, s.parent);
                }

                if (s.prev == nil) self.head = s.next else slots[s.prev].next = s.next;
                if (s.next != nil) slots[s.next].prev = s.prev;
                self.arcs[s.arc].slot = nil;
            }

            fn rotateUp(self: *Sweep, slot: u32) Allocator.Error!void {
                const slots = self.slots;
                const parent = slots[slot].parent;
                const grand = slots[parent].parent;

                if (slots[parent].left == slot) {
                    const inner = slots[slot].right;
                    slots[parent].left = inner;
                    if (inner != nil) slots[inner].parent = parent;
                    slots[slot].right = parent;
                } else {
                    const inner = slots[slot].left;
                    slots[parent].right = inner;
                    if (inner != nil) slots[inner].parent = parent;
                    slots[slot].left = parent;
                }
                slots[parent].parent = slot;
                slots[slot].parent = grand;
                if (grand == nil) {
                    self.root = slot;
                } else if (slots[grand].left == parent) {
                    slots[grand].left = slot;
                } else {
                    slots[grand].right = slot;
                }

                try self.touched.append(self.gpa, parent);
                try self.touched.append(self.gpa, slot);
            }

            fn setNode(self: *Sweep, left: u32, right: u32, shape: u32) Allocator.Error!u32 {
                var folded: T = if (T == void) {} else self.values[shape];
                if (T != void) {
                    if (left != nil) folded = combine(self.sets.items[left].folded, folded);
                    if (right != nil) folded = combine(folded, self.sets.items[right].folded);
                }
                const node: u32 = @intCast(self.sets.items.len);
                try self.sets.append(self.gpa, .{ .left = left, .right = right, .shape = shape, .folded = folded });
                return node;
            }

            /// Adds `shape` if it is absent and removes it if it is present.
            fn setToggle(self: *Sweep, set: u32, shape: u32) Allocator.Error!u32 {
                var t = set;
                while (t != nil) {
                    const n = self.sets.items[t];
                    if (n.shape == shape) return self.setRemove(set, shape);
                    t = if (shape < n.shape) n.left else n.right;
                }
                return self.setInsert(set, shape);
            }

            fn setInsert(self: *Sweep, set: u32, shape: u32) Allocator.Error!u32 {
                if (set == nil) return self.setNode(nil, nil, shape);
                const n = self.sets.items[set];
                if (mix(shape) > mix(n.shape)) {
                    var left: u32 = nil;
                    var right: u32 = nil;
                    try self.setSplit(set, shape, &left, &right);
                    return self.setNode(left, right, shape);
                }
                if (shape < n.shape) return self.setNode(try self.setInsert(n.left, shape), n.right, n.shape);
                return self.setNode(n.left, try self.setInsert(n.right, shape), n.shape);
            }

            fn setSplit(self: *Sweep, set: u32, shape: u32, left: *u32, right: *u32) Allocator.Error!void {
                if (set == nil) {
                    left.* = nil;
                    right.* = nil;
                    return;
                }
                const n = self.sets.items[set];
                var inner: u32 = nil;
                if (n.shape < shape) {
                    try self.setSplit(n.right, shape, &inner, right);
                    left.* = try self.setNode(n.left, inner, n.shape);
                } else {
                    try self.setSplit(n.left, shape, left, &inner);
                    right.* = try self.setNode(inner, n.right, n.shape);
                }
            }

            fn setRemove(self: *Sweep, set: u32, shape: u32) Allocator.Error!u32 {
                if (set == nil) return nil;
                const n = self.sets.items[set];
                if (n.shape == shape) return self.setMerge(n.left, n.right);
                if (shape < n.shape) return self.setNode(try self.setRemove(n.left, shape), n.right, n.shape);
                return self.setNode(n.left, try self.setRemove(n.right, shape), n.shape);
            }

            fn setMerge(self: *Sweep, low: u32, high: u32) Allocator.Error!u32 {
                if (low == nil) return high;
                if (high == nil) return low;
                const l = self.sets.items[low];
                const h = self.sets.items[high];
                if (mix(l.shape) > mix(h.shape)) return self.setNode(l.left, try self.setMerge(l.right, high), l.shape);
                return self.setNode(try self.setMerge(low, h.left), h.right, h.shape);
            }
        };
    };
}
