# Historical experiment 05 outputs

These CSVs and figures use the original measurement schema and remain unchanged.
They must not be presented as corrected post-contact union-LKC results.

Euler and Crofton surface were summed over nearest-centre partitions;
mesh curvature and surface were summed over individual ball preimages.
After the balls overlap, partition boundaries introduce internal interfaces,
and ball sums double-count overlap. VF and binary spanning remain distinct
union statistics. Pre-contact sums are useful controls but retain numerical
discretisation error.

The [corrected runner and schema](../README.md) measure the union directly,
report both 6- and 26-neighbour topology, and keep per-ball and partition
diagnostics in separate columns. A new full run is required to replace the
scientific interpretation of the post-contact curves; its default destination
is `results_union_v2/`.
