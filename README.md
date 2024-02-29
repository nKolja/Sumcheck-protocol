# Sumcheck-protocol
An implementation of the sumcehck protocol in Rust


The implementation provides a Prover struct, and a Proof struct.
The Prover can be used to create a Proof, which is a sumcheck non-interactive proof, which can be verified with the `verify_proof` function.
There is also the possiblity of an interactive proof with the prover by using the `verify_interactive` function.

The input to the prover is a multivariate polynomial from use `ark_poly::polynomial::multivariate`, however special attention was made to not use any representation-specific properties of the polynomial apart from the fact that it can be efficiently queried.
It is assumed that both the prover and the verifier have access to the polynomial only as an oracle.
The reason for this choice is because the main benefit of the sumcheck protocol stems from such use-cases.
For example such is the case when the product of individual variable degrees is very large, but the polynomial can still be efficiently queried (even though it cannot be efficiently represented).

The only trait of the multivariate polynomial that is used is `evaluate`, while other representation-specific traits are only used when generating a random polynomial for testing purposes.

The tests cover both interactive and non-interactive protocol for random dense and sparse polynomials with number of variables `n` randomly chosen and ranging between 1 and 8, and individual variable degrees ranging between 0 and 15.
There are also constant tests for the interactive and non-interactive protocol for the degree 3 polynomial from the book `Proofs, Arguments and Zero-Knowledge` by Justin Thaler, page 36, `https://people.cs.georgetown.edu/jthaler/ProofsArgsAndZK.pdf`. 

Note that Thaler never mentions the oracle access to the polynomial as a constrant which greatly changes the protocols utility. 
In particular the sum over the hybercube can be computed much faster then the generic method of evaluating the function at all points of the hybercube and summing the result.
