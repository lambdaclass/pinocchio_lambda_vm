# README
This is a self-contained Rust implementation of [Pinocchio](https://eprint.iacr.org/2013/279.pdf). It is intented for those who want to learn more about SNARKs. Pinocchio is a proving system for efficiently verifying general computations. This source code is the companion of this [blog post](https://blog.lambdaclass.com/pinocchio-virtual-machine-nearly-practical-verifiable-computation/
) were a more friendly and high level explanation of Pinocchio can be found.

It also contains an arkworks_adapter project to use circuits made with Arkworks R1CS with this implementation of Pinocchio
# Usage

To run tests:
`make test`
`make test_arkworks_adapter`

To use a docker shell:
`make docker-shell`

To use a nix shell:
`make nix-shell`

# Pinocchio
Pinocchio is a proving system. It is composed of an initial setup, a prover and a verifier. The idea is the following. For a computer program P written in some language, the prover runs the program P and computes its output and a proof of execution. This proof convinces the verifier that he has actually run the program and the output is correct. This is useful for example when the prover is an untrusted agent.

To achieve this, the system needs to work with mathematical expressions instead of normal traces of program executions. For this reason, correct execution instances are expressed in terms of polynomials in a process called arithmetization. 

`pinocchio_lambda_vm` uses elliptic curve cryptography, polynomials and other mathematical tools common to many other SNARKs.
 
The module structure is the following:
 - ```pinocchio/```: setup, prover, verifier.
 - ```math/```: elliptic curve, polynomials and other math primitives.
 - ```circuits/```: code relevant to the arithmetization process.

These concepts are not simple, but they're probably simpler than you think. For further reading we suggest:
- [Pinocchio](https://eprint.iacr.org/2013/279.pdf): original paper with formal explanation of the system.
- [Moonmath manual](https://leastauthority.com/community-matters/moonmath-manual/): friendly resource to understand the mathematical background of zk-SNARKs in general.
- [Pairing for beginners](https://static1.squarespace.com/static/5fdbb09f31d71c1227082339/t/5ff394720493bd28278889c6/1609798774687/PairingsForBeginners.pdf): a good introduction to elliptic curves in general and pairings in particular.
- [An Introduction to Mathematical Cryptography](https://link.springer.com/book/10.1007/978-0-387-77993-5): general cryptography and useful for pairings also.
- [Master thesis from David Møller Hansen](https://www.sagemath.org/files/thesis/hansen-thesis-2009.pdf): other helpful resource to understand pairings and their implementations.

Some of the algorithms in `pinocchio_lambda_vm`  can be found in these books. You will find references to this resources throughout the code.

# Understanding the code
We encourage to start by reading the [blog post](https://www.notamonadtutorial.com/pinocchio-virtual-machine-nearly-practical-verifiable-computation/
).

 Then, in the ```tests/``` module you will find integration tests that will guide you through the happy path: the setup, the prover and the verifier. Each of these components can be located in their own files:

- `pinocchio/setup.rs`: generates the `VerificationKey` and the `EvaluationKey` with the relevant data to construct and verify proofs. This data comes from the structure of the program P encoded as `Polynomials` in a `QAP` (Quadratic arithmetic program). To hide this data, random `FieldElement`s are sampled as `ToxicWaste` and then mapped to `EllipticCurveElement`s via repeated addition of a `generator()` of the curve. 

- `pinocchio/prover.rs`: takes the circuit encoded as a `QAP` and the trace of the program as `FieldElement`s. Then, it applies the `msm(...)` operation in order to generate the proof elements, in this case, `EllipticCurveElement`s hidings.

- `pinocchio/verifier.rs`: verifies the proof by checking the conditions mentioned in the paper. This involves computing `pairing(...)` operations between `EllipticCurveElement`'s.

The rest of the files implement mathematical primitives and other auxiliary tools to complete the system.
