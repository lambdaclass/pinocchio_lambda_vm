# Arkworks Adapter for Lambda's Pinocchio 

To generate a proof from a constraint system made in arkworks in Lambda's Pinocchio, call the function to import the r1cs, and the function to import the IO and witness 

```
let r1cs = arkworks_cs_to_pinocchio_r1cs(&cs);
let (io, witness) = arkworks_io_and_witness_to_pinocchio_io_and_witness(&cs);
```
