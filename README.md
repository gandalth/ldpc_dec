# Overview

ldpc_dec is a belief-propagation decoder written in Rust.
This is a WIP decoder which shall operate on both classical LDPC codes
and quantum LDPC (qLDPC) codes.

# Installation

Install Rust on your system and invoke
```
cargo run
```
to build and execute ldpc_dec.

# Implementation status

The decoder accepts soft information from the channel model for both
the variable nodes and the measured syndrome. This will allow us to
model different reliability levels for different bit positions in the
future, e.g. based on topology. Further, the soft information for the
syndrome measurement will allow for repeated syndrome measurements and
measurement updates while decoding. The joint decoding engine models
syndrome measurements as auxiliary variable nodes in the decoding
graph. This poses correct parity constraints while allowing for
syndrome measurements to be corrected as well.

To illustrate the decoding process, it currently uses a simple AWGN
channel model for the classical mode and a Quantum-BSC-channel which
introduces bit flips and syndrome measurement flips.

Decoder includes:

- Standard belief-propagation decoding for classical LDPC codes.
- Joint decoding of bit-flip errors and syndrome measurement errors.

# Future Work

This is a recently started project with additional channel models and
decoding variants to follow.

