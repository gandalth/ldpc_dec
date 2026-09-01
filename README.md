# Overview

ldpc_dec is a belief-propagation decoder written in Rust.
This is a WIP decoder which shall operate on both classical LDPC codes
and quantum LDPC (qLDPC) codes.

# Installation

Install Rust on your system and invoke
```
cargo run -- classic
```
to build and execute ldpc_dec for a first test run.

In order to select the simulation type (see below), invoke one of the following:

```
cargo run -- classic
cargo run -- jd
cargo run -- css
```

# Implementation status

The decoder can work in three different modes. First, as a classical
belief-propgation decoder for classical bits. In this case, it accepts
soft information from the channel and performs iterations. It reports
whether it found a valid codeword before exceeding the maximum number
of iterations and presents either this codeword or, in case of
failure, the quantized state of variable node information. A second
mode is called joint decoding. In this mode, the decoder accepts soft
information from the channel model for both the variable nodes and the
measured syndrome. The joint decoding engine models syndrome
measurements as auxiliary variable nodes in the decoding graph. This
poses correct parity constraints while allowing for syndrome
measurements to be corrected as well. Note that this for this mode,
the channel model can either impose bit flips (X errors) or phase
flips (Z errors), but not both simultaneously. The third mode is
called css mode and allows for using css-type codes in the
decoder. The orthogonality constraint of the decoders parity-check
matrices is checked and the decoder will run two instances: one for X
errors, and a second one for Z errors.

To illustrate the decoding process, it currently uses a simple AWGN
channel model for the classical mode and a Quantum-BSC-channel which
introduces bit flips and syndrome measurement flips. the Quantum-BSC
channel model is used for both joint decoding mode and css mode. Note
that for css mode, the channel instantations for X errors and Z errors
are currently uncorrelated.

Decoder includes:

- Standard belief-propagation decoding for classical LDPC codes.
- Joint decoding of bit-flip errors and syndrome measurement errors.
- CSS decoding for CSS-type codes.

# Future Work

This is a recently started project with additional channel models and
decoding variants to follow.

