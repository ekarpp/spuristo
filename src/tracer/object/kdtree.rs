use super::*;
use std::time::Instant;
use std::f64::INFINITY;
use crate::Axis;

#[cfg(test)]
mod kdtree_tests;

/// Triangle mesh constructed as a kD-tree
pub type Mesh = KdTree<Triangle>;

/// A k dimensional tree used to accelerate ray intersection calculations.
/// Implements a binary tree that splits a large mesh of objects to smaller
/// subobjects.
pub struct KdTree<T> {
    objects: Vec<T>,
    boundary: AaBoundingBox,
    root: Box<KdNode>,
    material: Material,
}

/// Implementation of a SAH based kd-tree
/// References:
/// [ekzhang/rpt](https://github.com/ekzhang/rpt/blob/master/src/kdtree.rs),
/// [fogleman/pt](https://github.com/fogleman/pt/blob/master/pt/tree.go),
/// [Article by Amsallem](https://www.flomonster.fr/articles/kdtree.html)
impl<T: Bounded> KdTree<T> {
    /// Constructs a kD-tree of the given objects with the given material.
    /// Should each object have their own material instead?
    pub fn new(objects: Vec<T>, material: Material) -> Self {
        let start = Instant::now();

        let indices = (0..objects.len()).collect();
        let bounds: Vec<AaBoundingBox> = objects
            .iter()
            .map(|obj| obj.bounding_box()).collect();
        let boundary = bounds
            .iter()
            .fold(AaBoundingBox::default(), |b1, b2| b1.merge(b2));
        let root = KdNode::construct(&objects, &bounds, &boundary, indices);

        println!(
            "Constructed kd-tree of {} triangles in {:#?}",
            objects.len(),
            start.elapsed()
        );

        Self {
            root,
            objects,
            boundary,
            material,
        }
    }

    /// Returns self uniformly scaled as an instance with largest dimension
    /// of bounding box scaled to 1.0
    pub fn to_unit_size(self) -> Box<Instance<Self>> {
        let AaBoundingBox { ax_min, ax_max } = self.bounding_box();

        let bb_dim = ax_max - ax_min;
        let s = 1.0 / bb_dim.max_element();
        self.scale(s, s, s)
    }

    fn hit_subtree(
        &self,
        node: &KdNode,
        r: &Ray,
        t_min: f64,
        t_max: f64,
        aabb: &AaBoundingBox,
    ) -> Option<Hit> {
        // extract split info or check for hit at leaf node
        let (axis, point, mut node_first, mut node_second) = match node {
            KdNode::Split(axis, point, left, right) => (*axis, *point, left, right),
            KdNode::Leaf(indices) => {
                let mut tt = t_max;
                let mut h = None;
                for idx in indices {
                    h = self.objects[*idx].hit(r, t_min, tt).or(h);
                    tt = h.as_ref().map_or(tt, |hit| hit.t);
                }
                return h.map(|mut h| {
                    h.object = self;
                    h
                });
            }
        };

        let t_split = (point - r.o(axis)) / r.d(axis);

        let (mut aabb_first, mut aabb_second) = aabb.split(axis, point);

        let left_first = r.o(axis) < point || (r.o(axis) == point && r.d(axis) <= 0.0);
        // intersect first the AABB that we reach first
        if !left_first {
            std::mem::swap(&mut aabb_first, &mut aabb_second);
            std::mem::swap(&mut node_first, &mut node_second);
        }

        let (t_start, t_end) = aabb.intersect(r);
        let (t_start, t_end) = (t_start.max(t_min), t_end.min(t_max));

        // PBR Figure 4.19 (a). we hit only the first aabb.
        if t_split > t_end || t_split <= 0.0 {
            self.hit_subtree(node_first, r, t_start, t_end, &aabb_first)
        // PBR Figure 4.19 (b). we hit only the second aabb.
        } else if t_split < t_start {
            self.hit_subtree(node_second, r, t_start, t_end, &aabb_second)
        } else {
            let h1 = self.hit_subtree(node_first, r, t_start, t_end, &aabb_first);

            /* if we hit something in the first AABB before the split, there
             * is no need to process the other subtree. */
            if h1.as_ref().filter(|h| h.t < t_split).is_some() {
                h1
            } else {
                let h2 = self.hit_subtree(node_second, r, t_start, t_end, &aabb_second);
                if h1.is_some() && h2.is_some() {
                    if h1 < h2 {
                        h1
                    } else {
                        h2
                    }
                } else {
                    h1.or(h2)
                }
            }
        }
    }
}

impl<T: Bounded> Bounded for KdTree<T> {
    fn bounding_box(&self) -> AaBoundingBox {
        self.boundary
    }
}

impl<T: Bounded> Object for KdTree<T> {
    fn material(&self) -> &Material {
        &self.material
    }

    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<Hit> {
        let (t_start, t_end) = self.boundary.intersect(r);
        let (t_start, t_end) = (t_start.max(t_min), t_end.min(t_max));
        // box missed / is behind
        if t_start > t_end {
            None
        } else {
            self.hit_subtree(&self.root, r, t_start, t_end, &self.boundary)
        }
    }

    fn sample_on(&self, _rand_sq: DVec2) -> DVec3 {
        // pick at random a triangle and do it
        todo!()
    }

    fn sample_towards(&self, _xo: DVec3, _rand_sq: DVec2) -> Ray {
        // pick at random a triangle and do it
        todo!()
    }
    fn sample_towards_pdf(&self, _ri: &Ray) -> f64 {
        // pick at random a triangle and do it
        todo!()
    }
}

const COST_TRAVERSE: f64 = 15.0;
const COST_INTERSECT: f64 = 20.0;
const EMPTY_BONUS: f64 = 0.2;

/// A node in the kD-tree. Can be either a plane split or a leaf node.
pub enum KdNode {
    /// X-split, axis (x = 0, y = 1, z = 2), split point and child nodes
    Split(Axis, f64, Box<KdNode>, Box<KdNode>),
    /// Stores indices to the object vector in the kD-tree
    Leaf(Vec<usize>),
}

impl KdNode {
    /// Computes the cost for split along `axis` at `point`.
    fn cost(
        boundary: &AaBoundingBox,
        axis: Axis,
        point: f64,
        num_left: usize,
        num_right: usize,
    ) -> f64 {
        if !boundary.cuts(axis, point) {
            INFINITY
        } else {
            let (left, right) = boundary.split(axis, point);

            let cost = COST_TRAVERSE + COST_INTERSECT *
                (num_left as f64 * left.area() / boundary.area()
                 + num_right as f64 * right.area() / boundary.area());

            if num_left == 0 || num_right == 0 {
                (1.0 - EMPTY_BONUS) * cost
            } else {
                cost
            }
        }
    }

    /// Finds the best split according to SAH.
    fn find_best_split(aabbs: &Vec<&AaBoundingBox>, boundary: &AaBoundingBox) -> (Axis, f64, f64) {
        let mut best_cost = INFINITY;
        let mut best_point = INFINITY;
        let mut best_axis = Axis::X;

        for axis in [Axis::X, Axis::Y, Axis::Z] {
            let mut mins: Vec<f64> = Vec::new();
            let mut maxs: Vec<f64> = Vec::new();

            aabbs.iter().for_each(|aabb| {
                mins.push(aabb.min(axis));
                maxs.push(aabb.max(axis));
            });

            mins.sort_by(|a: &f64, b: &f64| a.partial_cmp(b).unwrap());
            maxs.sort_by(|a: &f64, b: &f64| a.partial_cmp(b).unwrap());

            let mut num_left = 0;
            let mut num_right = aabbs.len();

            // iterator for mins
            let mut min_idx = 0;
            // iterator for maxs
            let mut max_idx = 0;

            // add infinity to end as "null"
            mins.push(INFINITY);
            maxs.push(INFINITY);

            // do quasi merge
            while mins[min_idx] < INFINITY || maxs[max_idx] < INFINITY {
                let is_min = mins[min_idx] <= maxs[max_idx];
                let point = mins[min_idx].min(maxs[max_idx]);

                // update objects on right before cost..
                if !is_min { max_idx += 1; num_right -= 1; }

                let cost = Self::cost(boundary, axis, point, num_left, num_right);

                if cost < best_cost {
                    best_cost = cost;
                    best_axis = axis;
                    best_point = point;
                }

                // ..and objects on left after cost
                if is_min { min_idx += 1; num_left += 1; }
            }
        }

        (best_axis, best_point, best_cost)
    }

    /// Partitions indices to left and right parts along `axis` at `point`
    fn partition(
        aabbs: &[&AaBoundingBox],
        indices: Vec<usize>,
        axis: Axis,
        point: f64,
    ) -> (Vec<usize>, Vec<usize>) {
        // fix size
        let mut left: Vec<usize> = Vec::new();
        let mut right: Vec<usize> = Vec::new();
        aabbs.iter().zip(indices).for_each(|(aabb, idx)| {
            if aabb.min(axis) < point {
                left.push(idx);
            }
            // are we missing one, since both check strict?
            if aabb.max(axis) > point {
                right.push(idx);
            }
        });
        (left, right)
    }

    /// Constructs nodes of the kD-tree recursively with SAH.
    pub fn construct<T: Bounded>(
        objects: &[T],
        bounds: &[AaBoundingBox],
        boundary: &AaBoundingBox,
        indices: Vec<usize>,
    ) -> Box<Self> {
        // filter relevant AABBs
        let aabbs: Vec<&AaBoundingBox> = indices.iter()
            .map(|idx| &bounds[*idx]).collect();
        let (axis, point, cost) = Self::find_best_split(&aabbs, boundary);
        // cut not worth it, make a leaf
        if cost > COST_INTERSECT * indices.len() as f64 {
            Box::new(Self::Leaf(indices))
        } else {
            let (left_idx, right_idx) = Self::partition(&aabbs, indices, axis, point);
            let (left_bound, right_bound) = boundary.split(axis, point);
            Box::new(Self::Split(
                axis,
                point,
                Self::construct(objects, bounds, &left_bound, left_idx),
                Self::construct(objects, bounds, &right_bound, right_idx),
            ))
        }
    }
}
