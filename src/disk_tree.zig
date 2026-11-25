//! Disk stabbing: preprocess a set of disks so that "which disks contain this
//! point?" is answered in O(log n + k) for k reported disks.
//!
//! A vertical line sweeps the plane left to right. The upper and lower arcs of
//! the circles it cuts are kept in a treap ordered by y, and every arc carries
//! the set of disks covering the face directly above it. That order only
//! changes where a circle starts, ends, or crosses another one, so there is one
//! version of the treap per such x. Both the treap and the disk sets are
//! persistent, which makes each new version cost O(log n) instead of a full
//! copy: treap nodes hold one timestamped override before they are copied, and
//! the disk sets are path-copied.
//!
//! A query binary searches the versions by x, walks one version by y, and
//! enumerates the disk set it lands on.

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const Circle = struct {
    cx: f64,
    cy: f64,
    r: f64,

    pub fn contains(c: Circle, x: f64, y: f64) bool {
        const dx = x - c.cx;
        const dy = y - c.cy;
        return dx * dx + dy * dy <= c.r * c.r;
    }

    /// Circles that fail this are ignored: they are never reported.
    fn valid(c: Circle) bool {
        return std.math.isFinite(c.cx) and std.math.isFinite(c.cy) and
            std.math.isFinite(c.r * c.r) and c.cx - c.r < c.cx + c.r;
    }
};

const nil: u32 = std.math.maxInt(u32);

/// An arc is `2 * disk` for the lower half of a circle, `2 * disk + 1` for the
/// upper half.
fn arcY(circles: []const Circle, arc: u32, x: f64) f64 {
    const c = circles[arc >> 1];
    const dx = x - c.cx;
    const h = @sqrt(@max(0.0, c.r * c.r - dx * dx));
    return if (arc & 1 == 0) c.cy - h else c.cy + h;
}

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
        /// Disks covering the face directly above this arc.
        set: u32,
    };

    fn at(n: *const SeqNode, version: u32) State {
        return if (n.since <= version) n.second else n.first;
    }
};

const Event = struct {
    x: f64,
    /// A disk id for start and end events, an arc for crossings.
    a: u32,
    /// `start`, `end`, or the other arc of a crossing.
    b: u32,

    const start = nil;
    const end = nil - 1;

    fn before(_: void, l: Event, r: Event) bool {
        return l.x < r.x;
    }
};

/// Emits a start and an end event per circle, and an event per point where two
/// circles cross. Crossings are only hints telling the sweep which arcs to
/// re-examine, so near-tangent pairs are included rather than risk missing one.
fn collectEvents(gpa: Allocator, circles: []const Circle, events: *std.ArrayList(Event)) Allocator.Error!void {
    const Box = struct {
        lo: f64,
        hi: f64,
        /// Centre and radius along the other axis.
        c: f64,
        r: f64,
        id: u32,

        fn before(_: void, l: @This(), r: @This()) bool {
            return l.lo < r.lo;
        }
    };

    var min_x = std.math.inf(f64);
    var max_x = -std.math.inf(f64);
    var min_y = std.math.inf(f64);
    var max_y = -std.math.inf(f64);
    for (circles) |c| {
        if (!c.valid()) continue;
        min_x = @min(min_x, c.cx - c.r);
        max_x = @max(max_x, c.cx + c.r);
        min_y = @min(min_y, c.cy - c.r);
        max_y = @max(max_y, c.cy + c.r);
    }
    // Scanning along the longer axis keeps fewer boxes open at once.
    const along_x = max_x - min_x >= max_y - min_y;

    var boxes: std.ArrayList(Box) = .empty;
    defer boxes.deinit(gpa);
    try boxes.ensureTotalCapacity(gpa, circles.len);
    for (circles, 0..) |c, i| {
        if (!c.valid()) continue;
        const id: u32 = @intCast(i);
        try events.append(gpa, .{ .x = c.cx - c.r, .a = id, .b = Event.start });
        try events.append(gpa, .{ .x = c.cx + c.r, .a = id, .b = Event.end });
        const along = if (along_x) c.cx else c.cy;
        const across = if (along_x) c.cy else c.cx;
        boxes.appendAssumeCapacity(.{ .lo = along - c.r, .hi = along + c.r, .c = across, .r = c.r, .id = id });
    }
    std.mem.sortUnstable(Box, boxes.items, {}, Box.before);

    var open: std.ArrayList(Box) = .empty;
    defer open.deinit(gpa);
    for (boxes.items) |box| {
        var k: usize = 0;
        while (k < open.items.len) {
            const other = open.items[k];
            if (other.hi < box.lo) {
                _ = open.swapRemove(k);
                continue;
            }
            if (@abs(other.c - box.c) <= other.r + box.r) {
                try appendCrossings(gpa, circles, other.id, box.id, events);
            }
            k += 1;
        }
        try open.append(gpa, box);
    }
}

fn appendCrossings(gpa: Allocator, circles: []const Circle, ia: u32, ib: u32, events: *std.ArrayList(Event)) Allocator.Error!void {
    const a = circles[ia];
    const b = circles[ib];
    const dx = b.cx - a.cx;
    const dy = b.cy - a.cy;
    const d2 = dx * dx + dy * dy;
    if (d2 == 0) return;

    const slack = 1.0 + 1e-9;
    const sum = a.r + b.r;
    const diff = a.r - b.r;
    if (d2 > sum * sum * slack or d2 * slack < diff * diff) return;

    // The crossing points are the foot of the radical line on the centre
    // line, offset by h along the perpendicular.
    const t = (a.r * a.r - b.r * b.r + d2) / (2 * d2);
    const h = @sqrt(@max(0.0, a.r * a.r / d2 - t * t));
    const mx = a.cx + t * dx;
    const my = a.cy + t * dy;

    inline for (.{ -h, h }) |s| {
        const px = mx + s * dy;
        const py = my - s * dx;
        try events.append(gpa, .{
            .x = px,
            .a = 2 * ia + @intFromBool(py > a.cy),
            .b = 2 * ib + @intFromBool(py > b.cy),
        });
    }
}

/// A disk tree without per-disk values.
pub const DiskTree = Folded(void, struct {
    fn combine(_: void, _: void) void {}
}.combine);

/// A disk tree that also carries a value of type `T` per disk and can combine
/// the values of all disks containing a point without visiting them.
///
/// `combine(a, b)` must be associative. It is always called with the disks
/// behind `a` at lower indices than those behind `b`, so it need not be
/// commutative: with disks listed back to front, alpha compositing works.
pub fn Folded(comptime T: type, comptime combine: fn (T, T) T) type {
    return struct {
        const Self = @This();

        /// Node of a persistent treap of disk ids.
        const SetNode = struct {
            left: u32,
            right: u32,
            disk: u32,
            /// The values of this subtree's disks, combined in index order.
            folded: T,
        };

        const Values = if (T == void) void else []const T;

        circles: []const Circle,
        /// Left edge of each slab, ascending.
        xs: []const f64,
        /// Root in `seq` of the arc order inside each slab. The slab's index is
        /// the version to read the nodes at.
        roots: []const u32,
        seq: []const SeqNode,
        sets: []const SetNode,

        /// Copies `circles`; ids reported by queries are indices into it.
        ///
        /// Time and space are O((n + c) log n) for c circle-circle crossings, after
        /// a pair search that is O(n log n) plus the number of circles overlapping
        /// each other along the longer axis.
        pub const init = if (T == void) initPlain else initFolded;

        fn initPlain(gpa: Allocator, circles: []const Circle) Allocator.Error!Self {
            return build(gpa, circles, {});
        }

        /// `values[i]` belongs to `circles[i]`. The values are only read here.
        fn initFolded(gpa: Allocator, circles: []const Circle, values: []const T) Allocator.Error!Self {
            std.debug.assert(values.len == circles.len);
            return build(gpa, circles, values);
        }

        fn build(gpa: Allocator, circles: []const Circle, values: Values) Allocator.Error!Self {
            std.debug.assert(circles.len < nil / 2);

            const owned = try gpa.dupe(Circle, circles);
            errdefer gpa.free(owned);

            var events: std.ArrayList(Event) = .empty;
            defer events.deinit(gpa);
            try collectEvents(gpa, owned, &events);
            std.mem.sortUnstable(Event, events.items, {}, Event.before);

            var sweep = try Sweep.init(gpa, owned, values);
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
                .circles = owned,
                .xs = xs_owned,
                .roots = roots_owned,
                .seq = seq_owned,
                .sets = sets_owned,
            };
        }

        pub fn deinit(self: *Self, gpa: Allocator) void {
            gpa.free(self.circles);
            gpa.free(self.xs);
            gpa.free(self.roots);
            gpa.free(self.seq);
            gpa.free(self.sets);
            self.* = undefined;
        }

        /// Ids of the disks containing (x, y), ascending. Caller frees.
        ///
        /// A point within rounding error of a circle's boundary may be reported
        /// either way for that circle.
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

        /// The values of the disks containing (x, y), combined in index order, or
        /// null if there are none. O(log n) however many disks that is.
        pub fn fold(self: *const Self, x: f64, y: f64) ?T {
            const set = self.locate(x, y);
            return if (set == nil) null else self.sets[set].folded;
        }

        /// The set of disks covering the face that (x, y) lies in.
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
                if (arcY(self.circles, n.arc, x) <= y) {
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
                try out.append(gpa, n.disk);
                t = n.right;
            }
        }

        /// Heap bytes held by the structure.
        pub fn bytes(self: *const Self) usize {
            return std.mem.sliceAsBytes(self.circles).len + std.mem.sliceAsBytes(self.xs).len +
                std.mem.sliceAsBytes(self.roots).len + std.mem.sliceAsBytes(self.seq).len +
                std.mem.sliceAsBytes(self.sets).len;
        }

        /// The sweep line. `slots` is a mutable treap with the same shape as the
        /// newest persistent version, plus parent and neighbour links so an arc can be
        /// found and moved without searching. Slots `2 * d` and `2 * d + 1` are created
        /// when disk d starts; crossings exchange the arcs held by two slots rather
        /// than moving the slots.
        const Sweep = struct {
            gpa: Allocator,
            circles: []const Circle,
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
            };

            fn init(gpa: Allocator, circles: []const Circle, values: Values) Allocator.Error!Sweep {
                const slots = try gpa.alloc(Slot, 2 * circles.len);
                errdefer gpa.free(slots);
                const arcs = try gpa.alloc(Arc, 2 * circles.len);
                @memset(slots, .{});
                @memset(arcs, .{});
                return .{ .gpa = gpa, .circles = circles, .values = values, .slots = slots, .arcs = arcs };
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

                for (events) |e| if (e.b == Event.end) try self.endDisk(e.a);
                for (events) |e| if (e.b == Event.start) try self.startDisk(e.a, x);
                for (events) |e| if (e.b != Event.start and e.b != Event.end) {
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

            fn less(self: *const Sweep, a: u32, b: u32) bool {
                const ya = arcY(self.circles, a, self.x_eval);
                const yb = arcY(self.circles, b, self.x_eval);
                return ya < yb or (ya == yb and a < b);
            }

            /// Crossing coordinates carry rounding error, so an arc touched at `x`
            /// stays dirty until the sweep is clear of where the event may really be.
            fn markDirty(self: *Sweep, arc: u32, x: f64) Allocator.Error!void {
                const state = &self.arcs[arc];
                if (state.slot == nil) return;
                const c = self.circles[arc >> 1];
                const scale = @max(c.r, @abs(c.cx), @abs(c.cy));
                state.expiry = @max(state.expiry, x + 1e-7 * @sqrt(c.r * scale));
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
                self.slots[lower].set = if (down & 1 == 0)
                    try self.setInsert(below, down >> 1)
                else
                    try self.setRemove(below, down >> 1);

                try self.touched.append(self.gpa, lower);
                try self.touched.append(self.gpa, upper);
            }

            fn startDisk(self: *Sweep, disk: u32, x: f64) Allocator.Error!void {
                const lo = 2 * disk;
                const hi = lo + 1;

                var below: u32 = nil;
                var t = self.root;
                while (t != nil) {
                    if (self.less(self.slots[t].arc, lo)) {
                        below = t;
                        t = self.slots[t].right;
                    } else {
                        t = self.slots[t].left;
                    }
                }
                const outside = if (below == nil) nil else self.slots[below].set;

                self.slots[lo] = .{ .arc = lo, .set = try self.setInsert(outside, disk) };
                self.slots[hi] = .{ .arc = hi, .set = outside };
                self.arcs[lo].slot = lo;
                self.arcs[hi].slot = hi;
                try self.insertAfter(below, lo);
                try self.insertAfter(lo, hi);
                try self.markDirty(lo, x);
                try self.markDirty(hi, x);
            }

            fn endDisk(self: *Sweep, disk: u32) Allocator.Error!void {
                const lo = self.arcs[2 * disk].slot;
                // Anything still between the two arcs leaves the disk here.
                while (self.slots[self.arcs[2 * disk + 1].slot].prev != lo) {
                    try self.swap(self.slots[self.arcs[2 * disk + 1].slot].prev);
                }
                const hi = self.arcs[2 * disk + 1].slot;
                try self.remove(hi);
                try self.remove(lo);
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

            fn setNode(self: *Sweep, left: u32, right: u32, disk: u32) Allocator.Error!u32 {
                var folded: T = if (T == void) {} else self.values[disk];
                if (T != void) {
                    if (left != nil) folded = combine(self.sets.items[left].folded, folded);
                    if (right != nil) folded = combine(folded, self.sets.items[right].folded);
                }
                const node: u32 = @intCast(self.sets.items.len);
                try self.sets.append(self.gpa, .{ .left = left, .right = right, .disk = disk, .folded = folded });
                return node;
            }

            fn setInsert(self: *Sweep, set: u32, disk: u32) Allocator.Error!u32 {
                if (set == nil) return self.setNode(nil, nil, disk);
                const n = self.sets.items[set];
                if (mix(disk) > mix(n.disk)) {
                    var left: u32 = nil;
                    var right: u32 = nil;
                    try self.setSplit(set, disk, &left, &right);
                    return self.setNode(left, right, disk);
                }
                if (disk < n.disk) return self.setNode(try self.setInsert(n.left, disk), n.right, n.disk);
                return self.setNode(n.left, try self.setInsert(n.right, disk), n.disk);
            }

            fn setSplit(self: *Sweep, set: u32, disk: u32, left: *u32, right: *u32) Allocator.Error!void {
                if (set == nil) {
                    left.* = nil;
                    right.* = nil;
                    return;
                }
                const n = self.sets.items[set];
                var inner: u32 = nil;
                if (n.disk < disk) {
                    try self.setSplit(n.right, disk, &inner, right);
                    left.* = try self.setNode(n.left, inner, n.disk);
                } else {
                    try self.setSplit(n.left, disk, left, &inner);
                    right.* = try self.setNode(inner, n.right, n.disk);
                }
            }

            fn setRemove(self: *Sweep, set: u32, disk: u32) Allocator.Error!u32 {
                if (set == nil) return nil;
                const n = self.sets.items[set];
                if (n.disk == disk) return self.setMerge(n.left, n.right);
                if (disk < n.disk) return self.setNode(try self.setRemove(n.left, disk), n.right, n.disk);
                return self.setNode(n.left, try self.setRemove(n.right, disk), n.disk);
            }

            fn setMerge(self: *Sweep, low: u32, high: u32) Allocator.Error!u32 {
                if (low == nil) return high;
                if (high == nil) return low;
                const l = self.sets.items[low];
                const h = self.sets.items[high];
                if (mix(l.disk) > mix(h.disk)) return self.setNode(l.left, try self.setMerge(l.right, high), l.disk);
                return self.setNode(try self.setMerge(low, h.left), h.right, h.disk);
            }
        };
    };
}

const testing = std.testing;

/// Checks `tree` against a linear scan at `queries` random points in `box`
/// (min x, min y, max x, max y), skipping points too close to a boundary to
/// have a well-defined answer.
fn expectMatchesScan(circles: []const Circle, box: [4]f64, queries: usize, random: std.Random) !void {
    const gpa = testing.allocator;
    var tree = try DiskTree.init(gpa, circles);
    defer tree.deinit(gpa);

    var got: std.ArrayList(u32) = .empty;
    defer got.deinit(gpa);
    var want: std.ArrayList(u32) = .empty;
    defer want.deinit(gpa);

    var done: usize = 0;
    points: while (done < queries) {
        const x = box[0] + random.float(f64) * (box[2] - box[0]);
        const y = box[1] + random.float(f64) * (box[3] - box[1]);

        want.clearRetainingCapacity();
        for (circles, 0..) |c, i| {
            if (!c.valid()) continue;
            const d = std.math.hypot(x - c.cx, y - c.cy);
            if (@abs(d - c.r) < 1e-9 * @max(1.0, c.r)) continue :points;
            if (d < c.r) try want.append(gpa, @intCast(i));
        }
        done += 1;

        got.clearRetainingCapacity();
        try tree.queryAppend(gpa, &got, x, y);
        try testing.expectEqualSlices(u32, want.items, got.items);
    }
}

fn expectQuery(tree: *const DiskTree, x: f64, y: f64, want: []const u32) !void {
    const got = try tree.query(testing.allocator, x, y);
    defer testing.allocator.free(got);
    try testing.expectEqualSlices(u32, want, got);
}

test "no circles" {
    var tree = try DiskTree.init(testing.allocator, &.{});
    defer tree.deinit(testing.allocator);
    try expectQuery(&tree, 0, 0, &.{});
}

test "overlapping, nested and disjoint" {
    const circles = [_]Circle{
        .{ .cx = 0, .cy = 0, .r = 2 },
        .{ .cx = 1, .cy = 0.5, .r = 1 },
        .{ .cx = -1, .cy = 0.5, .r = 1 },
        .{ .cx = 0, .cy = 0, .r = 5 },
        .{ .cx = 20, .cy = 0, .r = 1 },
    };
    var tree = try DiskTree.init(testing.allocator, &circles);
    defer tree.deinit(testing.allocator);

    try expectQuery(&tree, 0.5, 0.2, &.{ 0, 1, 3 });
    try expectQuery(&tree, -0.5, 0.3, &.{ 0, 2, 3 });
    try expectQuery(&tree, 0, 0.6, &.{ 0, 3 });
    try expectQuery(&tree, 0, -1.5, &.{ 0, 3 });
    try expectQuery(&tree, 3, 3, &.{3});
    try expectQuery(&tree, 20.5, 0.5, &.{4});
    try expectQuery(&tree, 10, 0, &.{});
    try expectQuery(&tree, -1000, 0, &.{});
    try expectQuery(&tree, 1000, 0, &.{});
    try expectQuery(&tree, 0, 1000, &.{});
    try expectQuery(&tree, 0, -1000, &.{});
    try expectQuery(&tree, std.math.nan(f64), 0, &.{});
    try expectQuery(&tree, 0, std.math.nan(f64), &.{});
}

test "invalid circles are ignored" {
    const circles = [_]Circle{
        .{ .cx = 0, .cy = 0, .r = 0 },
        .{ .cx = 0, .cy = 0, .r = -1 },
        .{ .cx = 0, .cy = 0, .r = std.math.nan(f64) },
        .{ .cx = 0, .cy = 0, .r = std.math.inf(f64) },
        .{ .cx = std.math.inf(f64), .cy = 0, .r = 1 },
        .{ .cx = 1e300, .cy = 0, .r = 1 },
        .{ .cx = 0, .cy = 0, .r = 1 },
    };
    var tree = try DiskTree.init(testing.allocator, &circles);
    defer tree.deinit(testing.allocator);
    try expectQuery(&tree, 0.1, 0.1, &.{6});
}

test "random circles" {
    var prng: std.Random.DefaultPrng = .init(1);
    const random = prng.random();
    var circles: [300]Circle = undefined;

    for (0..40) |round| {
        // Sparse to heavily overlapping.
        const max_r = 0.5 + @as(f64, @floatFromInt(round));
        const n = 1 + random.uintLessThan(usize, circles.len);
        for (circles[0..n]) |*c| c.* = .{
            .cx = random.float(f64) * 100,
            .cy = random.float(f64) * 100,
            .r = 0.01 + random.float(f64) * max_r,
        };
        try expectMatchesScan(circles[0..n], .{ -5, -5, 105, 105 }, 2000, random);
    }
}

test "tall and narrow" {
    var prng: std.Random.DefaultPrng = .init(2);
    const random = prng.random();
    var circles: [300]Circle = undefined;
    for (&circles) |*c| c.* = .{
        .cx = random.float(f64),
        .cy = random.float(f64) * 1000,
        .r = 1 + random.float(f64) * 20,
    };
    try expectMatchesScan(&circles, .{ -25, -25, 25, 1025 }, 5000, random);
}

test "mixed scales" {
    var prng: std.Random.DefaultPrng = .init(3);
    const random = prng.random();
    var circles: [300]Circle = undefined;
    for (&circles) |*c| c.* = .{
        .cx = 1e6 + random.float(f64) * 10,
        .cy = -1e6 + random.float(f64) * 10,
        .r = std.math.pow(f64, 10, random.float(f64) * 6 - 4),
    };
    try expectMatchesScan(&circles, .{ 1e6 - 1, -1e6 - 1, 1e6 + 11, -1e6 + 11 }, 5000, random);
}

test "integer grid" {
    // Many circles through common points, tangencies, and shared event xs.
    var prng: std.Random.DefaultPrng = .init(4);
    const random = prng.random();
    var circles: std.ArrayList(Circle) = .empty;
    defer circles.deinit(testing.allocator);

    for ([_]f64{ 0.5, 1, @sqrt(2.0), 2, 2.5 }) |r| {
        circles.clearRetainingCapacity();
        for (0..10) |i| for (0..10) |j| {
            try circles.append(testing.allocator, .{ .cx = @floatFromInt(i), .cy = @floatFromInt(j), .r = r });
        };
        try expectMatchesScan(circles.items, .{ -3, -3, 12, 12 }, 5000, random);
    }

    circles.clearRetainingCapacity();
    for (0..200) |_| try circles.append(testing.allocator, .{
        .cx = @floatFromInt(random.uintLessThan(u32, 12)),
        .cy = @floatFromInt(random.uintLessThan(u32, 12)),
        .r = @floatFromInt(1 + random.uintLessThan(u32, 4)),
    });
    try expectMatchesScan(circles.items, .{ -5, -5, 16, 16 }, 5000, random);
}

test "concentric and duplicate" {
    var prng: std.Random.DefaultPrng = .init(5);
    const random = prng.random();
    var circles: [200]Circle = undefined;
    for (&circles, 0..) |*c, i| c.* = .{ .cx = 3, .cy = -2, .r = @floatFromInt(1 + i / 2) };
    try expectMatchesScan(&circles, .{ -100, -105, 106, 101 }, 5000, random);
}

test "rings through one point" {
    // Every circle passes through the origin.
    var prng: std.Random.DefaultPrng = .init(6);
    const random = prng.random();
    var circles: [64]Circle = undefined;
    for (&circles, 0..) |*c, i| {
        const angle = std.math.tau * @as(f64, @floatFromInt(i)) / circles.len;
        c.* = .{ .cx = 5 * @cos(angle), .cy = 5 * @sin(angle), .r = 5 };
    }
    try expectMatchesScan(&circles, .{ -11, -11, 11, 11 }, 5000, random);
    try expectMatchesScan(&circles, .{ -0.01, -0.01, 0.01, 0.01 }, 5000, random);
}

test "fold matches a scan" {
    // Composition of affine maps mod a prime: associative, not commutative.
    const Affine = struct {
        a: u64,
        b: u64,

        const p = 1_000_003;

        fn then(f: @This(), g: @This()) @This() {
            return .{ .a = f.a * g.a % p, .b = (f.b * g.a + g.b) % p };
        }
    };
    const Tree = Folded(Affine, Affine.then);
    const gpa = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(7);
    const random = prng.random();

    var circles: [300]Circle = undefined;
    var values: [300]Affine = undefined;
    for ([_]f64{ 2, 10, 40 }) |max_r| {
        for (&circles, &values) |*c, *v| {
            c.* = .{ .cx = random.float(f64) * 100, .cy = random.float(f64) * 100, .r = 0.01 + random.float(f64) * max_r };
            v.* = .{ .a = 1 + random.uintLessThan(u64, Affine.p - 1), .b = random.uintLessThan(u64, Affine.p) };
        }
        var tree = try Tree.init(gpa, &circles, &values);
        defer tree.deinit(gpa);

        for (0..5000) |_| {
            const x = random.float(f64) * 100;
            const y = random.float(f64) * 100;
            const ids = try tree.query(gpa, x, y);
            defer gpa.free(ids);

            var want: ?Affine = null;
            for (ids) |id| want = if (want) |w| w.then(values[id]) else values[id];
            try testing.expectEqual(want, tree.fold(x, y));
        }
    }
}
