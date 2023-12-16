extern crate alloc;

use alloc::collections::BTreeSet;
use core::cmp::Ordering;

#[derive(Clone, Debug)]
pub struct TheilSen {
    points: BTreeSet<Point2>,
    slopes: Vec<f64>,
}

impl TheilSen {
    pub fn new() -> Self {
        Self {
            points: BTreeSet::new(),
            slopes: vec![],
        }
    }

    pub fn of(points: &[(f64, f64)]) -> Self {
        let mut this = Self {
            points: BTreeSet::new(),
            slopes: vec![],
        };

        this.extend(points);

        this
    }

    /// How many points are actually in the data set.
    pub fn cardinality(&self) -> usize {
        self.points.len()
    }

    /// Add a new point, expressed as an `(x, y)` pair,` to the data set. If a
    /// point contains any non-finite values, it is skipped.
    pub fn push(&mut self, point: (f64, f64)) {
        let (x, y) = point;

        if !(x.is_finite() && y.is_finite()) {
            return;
        }

        let p = Point2 { x, y };

        for before in self.points.iter().filter(|other| other.x < p.x) {
            let slope = (p.y - before.y) / (p.x - before.x);
            if !slope.is_finite() {
                continue;
            }
            self.slopes.push(slope);
        }

        for after in self.points.iter().filter(|other| other.x > p.x) {
            let slope = (after.y - p.y) / (after.x - p.x);
            if !slope.is_finite() {
                continue;
            }
            self.slopes.push(slope);
        }

        self.points.insert(p);
    }

    pub fn extend(&mut self, points: &[(f64, f64)]) {
        for p in points {
            self.push(*p);
        }
    }

    // Calculate median of slopes of all potential lines.
    fn slope(&self) -> f64 {
        if self.slopes.is_empty() {
            return f64::NAN;
        }

        let mut copy = self.slopes.clone();
        copy.sort_by(|f1, f2| f1.total_cmp(f2));

        let mid = copy.len() / 2;

        if copy.len() % 2 == 0 {
            let diff = copy[mid] - copy[mid - 1];
            copy[mid] - (diff / 2.0)
        } else {
            copy[mid]
        }
    }

    /// Retrieve best fit estimate, as a `(slope, intercept)` pair.
    pub fn estimate(&self) -> (f64, f64) {
        let median_slope = self.slope();

        let median_intercept = {
            let mut intercepts: Vec<f64> = self
                .points
                .iter()
                .map(|p| p.y - (p.x * median_slope))
                .collect();

            intercepts.sort_by(|f1, f2| f1.total_cmp(f2));

            let len = intercepts.len();

            match len {
                0 => f64::NAN,
                len if len % 2 == 0 => {
                    let mid = len / 2;
                    let left = intercepts[mid];
                    let diff = left - intercepts[mid - 1];
                    left - (diff / 2.0)
                }
                len => {
                    let mid = len / 2;
                    intercepts[mid]
                }
            }
        };

        (median_slope, median_intercept)
    }
}

impl Default for TheilSen {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone, Debug)]
struct Point2 {
    x: f64,
    y: f64,
}

impl Eq for Point2 {}

impl Ord for Point2 {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.x.total_cmp(&other.x) {
            Ordering::Equal => self.y.total_cmp(&other.y),
            inequal => inequal,
        }
    }
}

impl PartialEq for Point2 {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl PartialOrd for Point2 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_points_is_nan() {
        let points = vec![];
        let mut estimator = TheilSen::new();

        for p in points {
            estimator.push(p);
        }

        let (m, b) = estimator.estimate();

        assert!(m.is_nan());
        assert!(b.is_nan());
    }

    #[test]
    fn one_point_is_nan() {
        let points = vec![(1.0, 1.0)];
        let mut estimator = TheilSen::new();

        for p in points {
            estimator.push(p);
        }

        let (m, b) = estimator.estimate();

        assert!(m.is_nan());
        assert!(b.is_nan());
    }

    #[test]
    fn two_points_is_a_line() {
        let points = vec![(0.0, 0.0), (1.0, 1.0)];
        let mut estimator = TheilSen::new();

        for p in points {
            estimator.push(p);
        }

        let (m, b) = estimator.estimate();

        assert_eq!(m, 1.0);
        assert_eq!(b, 0.0);
    }

    #[test]
    fn moderately_perturbed_line() {
        let slope = 0.1;
        let intercept = -1.0;
        let mut points: Vec<(f64, f64)> = (0..10)
            .map(|x| x as f64)
            .map(|x| (x, intercept + (slope * x)))
            .collect();

        // Moving the y-values for less than 29% of observations should not
        // make a difference.
        points[3].1 *= 0.1;
        points[7].1 *= 0.1;

        let estimator = TheilSen::of(&points);

        let (m, b) = estimator.estimate();

        assert_eq!(m, slope);
        assert_eq!(b, intercept);
    }

    use proptest::collection::vec;
    use proptest::*;

    proptest! {
        #[test]
        fn any_floats_no_crash(
            points in vec((num::f64::ANY, num::f64::ANY), 0..16)
        ) {
            let mut estimator = TheilSen::new();

            estimator.extend(&points);

            let (m, b) = estimator.estimate();

            assert!(
                (m.is_nan() && b.is_nan()) ||
                (!m.is_nan() && !b.is_nan())
            );
        }
    }
}
