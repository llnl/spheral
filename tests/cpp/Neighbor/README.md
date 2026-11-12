# NeighborTestBase modes

`neighbor_test_base_test` runs all 20 GoogleTest cases in either mode. With no
arguments, it uses smaller populations intended to finish in under 10 seconds;
confirm the elapsed time on the machine and build configuration being used.
Pass `--full` to restore the original populations and sampling behavior.

| Fixture | Default | `--full` |
| --- | --- | --- |
| Random 1-D NodeLists | 40 / 100 / 20 nodes | 1,000 / 2,500 / 500 nodes |
| Random 2-D NodeLists | 80 / 200 / 40 nodes | 10,000 / 25,000 / 5,000 nodes |
| Random 3-D NodeLists | 64 / 160 / 96 nodes | 1,000 / 2,500 / 1,500 nodes |
| Cylindrical 2-D | 12 radial rings | 100 radial rings |

Both modes retain the spatial domains, smoothing tensor construction, and checks
of raw neighbor calls, connectivity neighbors, complete node-pair lists, overlap
neighbors, and intersections. Default node samples include both x boundaries of
each NodeList, and intersection samples include each available within-list and
cross-list category before adding random pairs. Random fixtures still check 10
nodes per list and 10 intersection pairs; the default cylindrical fixture checks
20 rather than 50. Overlap tests check two nodes per list in default mode,
including 3-D. Full mode retains the original zero sampled 3-D overlap checks.

From the existing CMake build directory, build the target:

```sh
cmake --build . --target neighbor_test_base_test
```

CTest runs the default mode:

```sh
ctest -R '^neighbor_test_base_test$' --output-on-failure
```

Run the original workload directly from that same build directory:

```sh
time ./tests/neighbor_test_base_test --full
```

GoogleTest options can be used with either mode:

```sh
./tests/neighbor_test_base_test --gtest_filter='NeighborRandom3d.*'
./tests/neighbor_test_base_test --full --gtest_filter='NeighborRandom3d.*'
```

Each run prints the selected mode. CTest does not forward arbitrary command-line
arguments to its registered test executable, so select `--full` by running the
executable directly. Node positions and sampled checks remain randomized.
