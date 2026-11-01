# Bounded finite-six-spinon catalogues and scans

[Experimental build](experimental.md) · [Selected states and limitations](experimental-xxx-six-spinon.md)

The headers `bethe/experimental/xxx_six_spinon_catalogue.hpp` and
`xxx_six_spinon_scan.hpp` add integer label enumeration and owned Async scans
around the selected-state solver. They retain its restrictions: periodic XXX,
`J=1`, a singlet reference, `N>=12` divisible by four, and the two-two-string
and three-string families. Enumeration does not make an unresolved state
solvable, or establish complete coverage of the physical six-spinon sector.

## Momentum and identity

```cpp
namespace six = bethe::experimental::xxx::six_spinon;
six::MomentumCatalogue catalogue(12, 3); // q=2*pi*3/12=pi/2
auto request = catalogue.candidate(0);
auto id = six::candidate_id(request);
auto requests = catalogue.requests(100); // Refuse before solving if more than 100.
```

The catalogue accepts the **Fourier momentum transfer** `q=P0-Pn`, as an integer
index in `[0,N)`. The supported ground states have `P0=0`. The root state's
`momentum_index` remains its absolute `Pn`, so these indices have opposite signs
modulo `N`. They are not interchangeable.

Candidate ordering is family `22` then `3`, increasing string labels, then
lexicographic holes. All six holes and the one or two string labels are retained.
IDs such as `N12-22-h-6,-4,-2,0,2,6-b-1,3` store **twice** each quantum number;
they do not depend on precision, energy ordering, momentum selection or completion
order. Every candidate has multiplicity one: no reflection or parity folding is
applied. A scan-local line `state_id` is its one-based record position, not this
physical identity.

Integer suffix counts permit direct random access at one momentum without
enumerating other momenta. For example, `MomentumCatalogue(64,16)` counts
166,849 nominal candidates (103,817 in `22`, 63,032 in `3`) without solving any
roots or materializing every request. Its default `requests()` limit is 10,000,
so that materialization is refused. Call `candidate(i)` to select a bounded
subset or disjoint ordinal sets; those retain the same physical IDs.

The constructor's `max_table_entries` limits stored `size_t` cells, with a
default of 4,000,000. It checks the required `7*N*(N+5)` cells before allocation.
The recurrence takes comparable work; integer overflow raises rather than
saturating a count. This storage budget is independent of the request-count
budget. No large-N numerical capability follows from being able to count labels.

## Owned asynchronous evaluation

```cpp
namespace finite = bethe::experimental::xxx;
auto reference = bethe::heisenberg::sector_ground_state<double>(12, uni20::half_int{0});
finite::PreparedSinglet<double> ground(12, reference);
uni20::async::TbbScheduler scheduler(2);
auto pending = six::scan_async(std::move(ground), std::move(requests),
                               {&scheduler, 1, 2});
auto const& result = pending.get_wait(scheduler);
```

Include `uni20/async/tbb_scheduler.hpp` and link `TBB::tbb` for this scheduler.
The returned
`Async<ScanResult<Real>>` can also feed another task through `read()`; synchronous
`scan(...)` is a convenience wrapper. Both own their prepared-ground snapshot
and canonical requests. The scheduler must outlive completion. Dropping the
result handle does not cancel work. Concurrent runners pull new requests as
they finish, using the shared execution facility without intermediate batches.

Invalid requests, mixed system sizes, duplicate labels and an exceeded candidate
budget are rejected before submission. The explicit task-local precision of the
prepared singlet controls solving, scratch and form factors, including MPFR.
The scan does not change ambient precision or silently retry in another type.

Each record retains its physical ID, the solver request/result and diagnostic.
`factor` is absent if no form factor was attempted; `line` is absent for any
rejected state or form factor. Accepted lines use the positive excitation gap,
Fourier transfer and normalized **raising** weight. For a longitudinal singlet
response, divide those weights and their moments by two. Failed states are not
zero-weight lines. The N-bin `moments` contain only accepted lines; they do not
estimate missing intensity. `converged()` means all submitted requests produced
lines, including the vacuous empty scan, not full physical-sector completeness.

Cancellation stops new admission. Already active numerical kernels finish before
the operation reports `ExecutionCancelled`; execution exceptions similarly wait
for admitted work before propagation. An interrupted operation does not publish
a partially filled `ScanResult`. Ordinary numerical failures instead remain in
the completed result as explicit records. This is an in-memory scan, not a
persistent restart controller.

## Fresh example and checks

In a clean public source build with `BETHE_BUILD_EXPERIMENTAL=ON` and
`BETHE_BUILD_EXAMPLES=ON`, build and run:

```sh
cmake --build . --target bethe_example_six_spinon_scan
./experimental/bethe_example_six_spinon_scan
```

It enumerates and evaluates the 11 nominal `N=12`, `q=pi/2` candidates, prints
every record and the accepted moments, and returns 2 if any candidate is rejected
(0 if all requests succeed). With the current fp64 solver, nine candidates
produce lines and two retain unresolved-seed failures. No archive or cluster
account is needed.

Tests compare the integer catalogue against independent exhaustive labels and
occupied-quantum-number momenta for `N=12,16,20`; random-access counts/IDs at
`N=64` also match the independent integer reference. Numerical checks reuse
independent high-precision fixtures in both families and the two mapped
non-self-conjugate `N=12` cases. Serial/concurrent results, explicit MPFR working
precision, a known unresolved `N=20` state, duplicate/budget refusal, task
composition and admitted-work drainage are tested separately. The selected-state
guide describes the original-equation and form-factor acceptance criteria.
