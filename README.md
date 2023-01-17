# Pinocchio Lambda VM

Zero Knowledge Virtual Machine from scratch implementing Pinocchio. You can read more about the implementation in our post [Pinocchio Virtual Machine: Nearly Practical Verifiable Computation
](https://www.notamonadtutorial.com/pinocchio-virtual-machine-nearly-practical-verifiable-computation/).

To run all tests:

```make test```

To start a shell in docker to develop:

```make docker-shell```
or
```make docker-cuda-shell```

To start a shell using nix to develop:

```make nix-shell```
or
```make nix-cuda-shell```
