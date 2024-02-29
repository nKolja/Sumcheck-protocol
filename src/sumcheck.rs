#[allow(unused_imports)]
use ark_ff::{Field, UniformRand};
use ark_bls12_381::Fr as F;
use ark_std::{One, Zero};

// use ark_ff_optimized::fp64::Fp as F; //doesn;t work

use ark_poly::polynomial::univariate::SparsePolynomial as UniSparsePolynomial;
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
#[allow(unused_imports)]
use ark_poly::polynomial::{MVPolynomial, UVPolynomial};
use ark_poly::polynomial::Polynomial;


pub type Poly = UniSparsePolynomial<F>;
pub type MvPoly = SparsePolynomial<F, SparseTerm>;

// use std::collections::hash_map::DefaultHasher;
// use std::hash::Hash;
use blake3::Hasher;


// ADDED FOR EASIER PRINTING AND DEBUGGING

use std::fmt;

fn to_subscript(n: usize) -> String {
    let subs = ["₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"];
    n.to_string().chars().map(|c| subs[c.to_digit(10).unwrap() as usize]).collect()
}

fn to_superscript(n: usize) -> String {
    let sups = ["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"];
    n.to_string().chars().map(|c| sups[c.to_digit(10).unwrap() as usize]).collect()
}

pub struct DisplayMvPoly<'a>(&'a MvPoly);
pub struct DisplayPoly(Poly);

impl<'a> fmt::Display for DisplayMvPoly<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (coeff, term) in self.0.terms() {
            write!(f, "{} ", coeff)?;
            for (var, power) in term.iter() {
                write!(f, "x{}{} ", to_subscript(*var), to_superscript(*power))?;
            }
            write!(f, "+ \n")?;
        }
        Ok(())
    }
}

impl fmt::Display for DisplayPoly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (degree, coeff) in self.0.iter() {
            write!(f, "{} ", coeff)?;
            write!(f, "x{} ", to_superscript(*degree))?;
            write!(f, "+ \n")?;
        }
        Ok(())
    }
}






fn legendre_interpolation(x_values: &[F], y_values: &[F]) -> Poly {
    assert_eq!(x_values.len(), y_values.len());

    let dd = x_values.len();
    let mut interpolated_poly = Poly::zero();

    for i in 0..dd {
        // Polynomial that is 1 at i, and 0 elsewhere
        let mut delta_i = Poly::from_coefficients_vec(vec![(0, F::one())]);
        let mut denom = F::one();
        for j in 0..dd {
            if i != j {
                let x_min_xj = Poly::from_coefficients_vec(vec![(0, -x_values[j]), (1, F::one())]);
                denom *= x_values[i] - x_values[j];
                delta_i = delta_i.mul(&x_min_xj);
            }

        }
        delta_i = &delta_i * denom.inverse().unwrap();

        interpolated_poly = interpolated_poly + (&delta_i * y_values[i]);
    }

    interpolated_poly

}

// Integer into vector of field elements 0 or 1
// reverse representation
// 1 = 0b00...001 -> [1, 0, 0, ...]
// assumed i bounded by usize, so n < 64
fn int_to_vec(i: usize, n: usize) -> Vec<F> {
    (0..n)
    .map(|bit| 
    {
        if ((i >> bit) & 0x01) == 0x01
            {F::one()}
        else
            {F::zero()}
    }
    ).collect()
}




fn sum_over_cube(f: &MvPoly) -> F {
    let n = f.num_vars();
    let mut out = F::zero();

    for i in 0..(1 << n) {
        out += f.evaluate(&int_to_vec(i, n));
    }

    out
}



fn eval_sum_to_uni(f: &MvPoly, r: &Vec<F>, d: &usize) -> Poly {
    let mut evals = vec![F::zero(); d+1];
    let n = f.num_vars();

    for a in 0..=*d {
        for i in 0..(1 << (n - r.len() - 1)) {

            let mut temp_eval = Vec::with_capacity(n);
            let bin_rep = int_to_vec(i, n - r.len() - 1);

            temp_eval.extend_from_slice(&r);
            temp_eval.push(F::from(a as u64));
            temp_eval.extend_from_slice(&bin_rep);

            evals[a] += f.evaluate(&temp_eval);
        }
    }

    let domain: Vec<F> = (0..=*d).map(|i| F::from(i as u64)).collect();
    let polynomial = legendre_interpolation(&domain, &evals);

    polynomial
}



pub fn hash_poly_to_field(g: &Poly) -> F {
    // Convert the polynomial g into a string
    let g_string = format!("{:?}", g);

    // Compute the hash of g using Blake3
    let mut hasher = Hasher::new();
    hasher.update(g_string.as_bytes());
    let g_hash_bytes = hasher.finalize().as_bytes().to_vec();

    // Convert the hash into a field element, retrying with a new hash of the hash if necessary
    let mut g_hash_f: Option<F> = None;
    let mut current_hash_bytes = g_hash_bytes.clone();
    while g_hash_f.is_none() {
        g_hash_f = F::from_random_bytes(&current_hash_bytes);
        if g_hash_f.is_none() {
            let mut temp_hasher = Hasher::new();
            temp_hasher.update(&current_hash_bytes);
            current_hash_bytes = temp_hasher.finalize().as_bytes().to_vec();
        }
    }
    g_hash_f.unwrap()
}


pub struct Prover {
    pub f: MvPoly,      // Multivariate polynomial
    pub d: Vec<usize>,  // Maximal degree for each variable
    pub c: F,           // The output c = sum_{b in {0,1}^n} f(b)
    pub r: Vec<F>,      // Vector of random values that is populated during protocol
}

// f(pt) = c
pub struct Proof {
    pub f: MvPoly,      // Multivariate polynomial
    pub d: Vec<usize>,  // Maximal degree for each variable
    pub c: F,           // The output c = sum_{b in {0,1}^n} f(b)

    pub g: Vec<Poly>,   // Vector of univariate polynomial outputs
}


impl Prover {
    pub fn new(f: &MvPoly, d: &Vec<usize>) -> Self {
        Self{
            f: f.clone(), 
            d: d.clone(),
            c: sum_over_cube(&f),
            r: vec![],
        }
    }

    fn response_gi(&mut self, ri: Option<&F>) -> Poly {
        if let Some(val) = ri {
            self.r.push(*val);
        }
        let gi = eval_sum_to_uni(&self.f, &self.r, &self.d[self.r.len()]);
        gi
    }


    pub fn generate_proof(&self) -> Proof{
        let n = self.f.num_vars();
        let mut r = vec![];
        let mut g = vec![];

        for i in 0..n {
            let gi = eval_sum_to_uni(&self.f, &r, &self.d[i]);
            r.push(hash_poly_to_field(&gi));
            g.push(gi);

        }

        Proof {
            f: self.f.clone(),
            d: self.d.clone(),
            c: self.c.clone(),
            g: g,
        }
    }

}







pub fn verify_proof(proof: &Proof) -> bool {
        // Retrieve the information from the proof
        let f =  &proof.f;
        let c =  &proof.c;
        let d =  &proof.d;
        let g =  &proof.g;
        let n = f.num_vars();

        let mut ret = true;
        let r: Vec<F> = g.iter().map(|gi| hash_poly_to_field(&gi)).collect();

        // Check that p0(0) + p0(1) = c
        ret = ret && (g[0].evaluate(&F::zero()) + g[0].evaluate(&F::one()) == *c) && (g[0].degree() <= d[0]);

        // Check that p_i+1(0) + p_i+1(1) = p_i(r_i) for all i
        for i in 0..n-1 {
            ret = ret && (g[i+1].evaluate(&F::zero()) + g[i+1].evaluate(&F::one()) == g[i].evaluate(&r[i]) && (g[i+1].degree() <= d[i+1]));

        }

        // Check that f(pt) = p_n(r_n) and deg(g) <= d
        ret = ret && (f.evaluate(&r) == g[n-1].evaluate(&r[n-1])) && (g[n-1].degree() <= d[n-1]);

        // If all checks pass, the proof is valid
        ret
}



pub fn verify_interactive(the_prover: &mut Prover) -> bool {
    let mut rng = ark_std::rand::thread_rng();
    
    // Start. Receive f and c
    let f =  &the_prover.f.clone();
    let d =  &the_prover.d.clone();
    let c =  &the_prover.c.clone();

    let n = f.num_vars();

    let mut g = vec![];
    let mut r = vec![];

    let mut ret = true;

    // First round. Receive first polynomial.
    g.push(the_prover.response_gi(None));
    ret = ret && (g[0].evaluate(&F::zero()) + g[0].evaluate(&F::one()) == *c) && (g[0].degree() <= d[0]);
    r.push(F::rand(&mut rng));


    // Send randomness and receive next polynomial
    for i in 0..(n-1) {
        g.push(the_prover.response_gi(Some(&r[i])));
        ret = ret && (g[i+1].evaluate(&F::zero()) + g[i+1].evaluate(&F::one()) == g[i].evaluate(&r[i])) && (g[i+1].degree() <= d[i+1]);
        r.push(F::rand(&mut rng));
    }

    ret = ret && !(g[n-1].evaluate(&r[n-1]) != f.evaluate(&r)) && (g[n-1].degree() <= d[n-1]);

    ret
}


pub fn random_mv_poly(d: &Vec<usize>) -> MvPoly {
    let mut rng = ark_std::rand::thread_rng();
    let n = d.len();
    let mut ind = vec![0; n];

    let mut monomials = Vec::new();

    // Loop through product (0..=d[0])(0..=d[1])...(0..=d[n-1])
    loop {
        let term = Term::new(ind.iter().enumerate().map(|(i, v)| (i, *v)).collect());
        monomials.push((F::rand(&mut rng), term));

        let mut i = 0;
        while i < n {
            ind[i] += 1;
            if ind[i] <= d[i] {
                break;
            }
            ind[i] = 0;
            i += 1;
        }
        if i == n {
            break;
        }
    }

    MvPoly::from_coefficients_vec(n, monomials)
}



pub fn assert_degrees(f: &MvPoly, d: &Vec<usize>) -> bool {
    let n = f.num_vars();
    let mut rng = ark_std::rand::thread_rng();

    let mut ret = true;

    for i in 0..n {
        // In theory i have to check they're all different, but this happens with negligible probability
        // would have to be adjusted for small fields, but I can't be bothered
        // Also not sure if x=0..d[i]+2 is secure or not, so i did random.
        let x_values: Vec<F> = (0..d[i]+2).map(|i| F::from(i as u64)).collect();
        let mut pt = vec![(0..n).map(|_| F::rand(&mut rng)).collect::<Vec<F>>(); d[i]+2];
        for j in 0..d[i]+2 {
            pt[j][i] = x_values[j];
        }
        let y_values: Vec<_> = pt.iter().map(|p| f.evaluate(p)).collect();

        let gi = legendre_interpolation(&x_values, &y_values);
        ret = ret && (gi.degree() <= d[i]);
    }

    ret
}



// pub fn main() {

// }



