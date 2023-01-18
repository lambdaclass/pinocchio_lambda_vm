# Arkworks Adapter for Lambda's Pinocchio 

To generate a proof from a constraint system made in arkworks in Lambda's Pinocchio, call the function to import the r1cs, and the function to import the IO and witness 

```
let r1cs = pinocchio_r1cs_from_arkworks_cs(&cs);
let (io, witness) = pinocchio_io_and_witness_from_arkworks_cs(&cs);
```
