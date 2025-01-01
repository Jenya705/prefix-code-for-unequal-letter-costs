use std::{cmp::Reverse, collections::BinaryHeap};

pub type TreeIndex = u128;

pub trait TreeSolver {
    fn new(letters: &[u32]) -> Self;

    fn solve(&mut self, lengths: &mut [u32], indices: &mut [TreeIndex]) -> Option<bool>;
}

pub struct FirstToFindTreeSolver {
    letters: Vec<u32>,
    heap: BinaryHeap<(Reverse<u32>, u128)>,
}

impl TreeSolver for FirstToFindTreeSolver {
    fn new(letters: &[u32]) -> Self {
        Self {
            letters: letters.to_vec(),
            heap: BinaryHeap::new(),
        }
    }

    fn solve(&mut self, mut lengths: &mut [u32], mut indices: &mut [TreeIndex]) -> Option<bool> {
        let mut length_changed = false;

        self.heap.clear();
        self.heap.push((Reverse(0), 1));

        while let Some((Reverse(cost), mut i)) = self.heap.pop() {
            if cost < lengths[0] {
                i *= self.letters.len() as TreeIndex;
                for &letter in self.letters.iter() {
                    self.heap.push((Reverse(cost + letter), i));
                    i += 1;
                }
            } else {
                length_changed = length_changed || lengths[0] != cost;
                indices[0] = i;
                lengths[0] = cost;
                lengths = &mut lengths[1..];
                indices = &mut indices[1..];
                if lengths.len() == 0 {
                    break;
                }
            }
        }

        Some(length_changed)
    }
}
