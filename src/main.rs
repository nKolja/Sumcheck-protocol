// src/main.rs
use sumcheck::sumcheck::{hash_poly_to_field, random_mv_poly, verify_proof, verify_interactive, assert_degrees};
use sumcheck::sumcheck::{Poly, MvPoly, Prover}; // Ensure this is correctly referenced based on your actual module and types
use ark_poly::polynomial::multivariate::Term;
use ark_poly::polynomial::MVPolynomial;
use ark_bls12_381::Fr as F;

fn main() {


    // Test hashing function is ok

    // Example polynomial creation (adapt based on your Poly definition and constructors)
    let poly1 = Poly::from_coefficients_vec(vec![(0, F::from(1u64)), (1, F::from(1u64))]); // Represents 1 + x
    let poly2 = Poly::from_coefficients_vec(vec![(2, F::from(1u64)), (1, F::from(3u64))]); // Represents 1 + 2x
    let poly3 = Poly::from_coefficients_vec(vec![(0, F::from(1u64)), (1, F::from(1u64))]); // Same as poly1

    // Hash the polynomials
    let hash1 = hash_poly_to_field(&poly1);
    let hash2 = hash_poly_to_field(&poly2);
    let hash3 = hash_poly_to_field(&poly3);

    // Print the hashes to inspect them
    // println!("Hash1: {:?}", hash1);
    // println!("Hash2: {:?}", hash2);
    // println!("Hash3: {:?}", hash3);

    // Example checks (adapt these checks to your needs)
    assert_ne!(hash1, hash2, "Hashes of different polynomials should not match.");
    assert_eq!(hash1, hash3, "Hashes of identical polynomials should match.");






    // Test sumcheck proof is ok
    
    let c =  Term::new(vec![]);                          // constant term
    let x0 = Term::new(vec![(1, 1), (0, 1)]);           // x0 * x1
    let x1 = Term::new(vec![(1, 2), (0, 1), (2,2)]);    // x0 * x1^2 * x2^2
    let x2 = Term::new(vec![(2, 3)]);                   // x2^3

    let cffs = vec![
        (F::from(2u64), c), // Random constant term
        (F::from(1u64), x0),    // Coefficient of x0 * x1 is 1
        (F::from(3u64), x1),    // Coefficient of x0 * x1^2 * x2^2 is 1
        (F::from(2u64), x2),    // Coefficient of x2^3 is 2
    ];


    let f = MvPoly::from_coefficients_vec(3, cffs.clone());

    let mut prover_f = Prover::new(&f, &vec![1,2,3]);
    let proof_f = prover_f.generate_proof();

    assert!(verify_proof(&proof_f));
    assert!(verify_interactive(&mut prover_f));

    // let n:  usize = 8;
    // let d0: usize = 3;
    let d: Vec<usize> = vec![1, 3, 0, 2, 2, 3, 2, 1];

    let g = random_mv_poly(&d);
    // d[1] = 2; 
    let mut prover_g = Prover::new(&g, &d);
    let proof_g = prover_g.generate_proof();

    assert!(verify_proof(&proof_g));
    assert!(verify_interactive(&mut prover_g));
    assert!(assert_degrees(&g, &d));

}
