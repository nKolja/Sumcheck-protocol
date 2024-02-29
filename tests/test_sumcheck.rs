
#[cfg(test)]
mod sumcheck_tests {

use sumcheck::sumcheck::{MvPoly, Prover};
use sumcheck::sumcheck::{random_mv_poly, verify_proof, verify_interactive, assert_degrees};
use ark_poly::polynomial::multivariate::Term;
use ark_poly::polynomial::MVPolynomial;
use rand::{Rng, thread_rng};



#[test]
fn interactive_dense_poly() {
    let mut rng = thread_rng();
    let n: usize = rng.gen_range(1..=6);
    let d: Vec<usize> = (0..n).map(|_| rng.gen_range(0..=2)).collect();

    let f = random_mv_poly(&d);

    let mut prover = Prover::new(&f, &d);

    assert!(verify_interactive(&mut prover));
}

#[test]
fn interactive_sparse_poly() {
    let mut rng = thread_rng();
    let n: usize = rng.gen_range(1..=8);
    let d = vec![rng.gen_range(1..=15); n];

    let f = MvPoly::rand(d[0], n, &mut rng);
    println!("{:?}", f);

    let mut prover = Prover::new(&f, &d);

    assert!(verify_interactive(&mut prover));
}

#[test]
fn interactive_constant_poly() {

    // f = 2x0^3 + x1 x3 + x2 x3
    let t0 = Term::new(vec![(0, 3)]);            // 2 x0^3
    let t1 = Term::new(vec![(0, 1), (2, 1)]);    // x0 x2
    let t2 = Term::new(vec![(1, 1), (2, 1)]);    // x1 x2

    let cffs = vec![
        (2u64.into(), t0),    // 2 * x0^3
        (1u64.into(), t1),    // x0 * x2
        (1u64.into(), t2),    // x1 * x2
    ];

    let d = vec![3, 1, 1];      // individual max degrees

    let f = MvPoly::from_coefficients_vec(3, cffs.clone());

    let mut prover = Prover::new(&f, &d);

    assert!(verify_interactive(&mut prover));
}





#[test]
fn non_interactive_dense_poly() {
    let mut rng = thread_rng();
    let n: usize = rng.gen_range(1..=6);
    let d: Vec<usize> = (0..n).map(|_| rng.gen_range(0..=2)).collect();

    let f = random_mv_poly(&d);

    let prover = Prover::new(&f, &d);
    let proof = prover.generate_proof();

    assert!(verify_proof(&proof));
}

#[test]
fn non_interactive_sparse_poly() {
    let mut rng = thread_rng();
    let n: usize = rng.gen_range(1..=8);
    let d = vec![rng.gen_range(1..=15); n];

    let f = MvPoly::rand(d[0], n, &mut rng);

    let prover = Prover::new(&f, &d);
    let proof = prover.generate_proof();
    
    assert!(verify_proof(&proof));
}

#[test]
fn non_interactive_constant_poly() {

    // f = 2x0^3 + x1 x3 + x2 x3
    let t0 = Term::new(vec![(0, 3)]);            // 2 x0^3
    let t1 = Term::new(vec![(0, 1), (2, 1)]);    // x0 x2
    let t2 = Term::new(vec![(1, 1), (2, 1)]);    // x1 x2

    let cffs = vec![
        (2u64.into(), t0),    // 2 * x0^3
        (1u64.into(), t1),    // x0 * x2
        (1u64.into(), t2),    // x1 * x2
    ];

    let d = vec![3, 1, 1];      // individual max degrees

    let f = MvPoly::from_coefficients_vec(3, cffs.clone());

    let prover = Prover::new(&f, &d);
    let proof = prover.generate_proof();
    
    assert!(verify_proof(&proof));

}





#[test]
fn mv_poly_check_degrees() {

    let mut rng = thread_rng();
    let n: usize = rng.gen_range(1..=6);
    let d: Vec<usize> = (0..n).map(|_| rng.gen_range(0..=2)).collect();

    let f = random_mv_poly(&d);


    assert!(assert_degrees(&f, &d));
}



}

