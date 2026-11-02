# Resuming spinon density controls

[Spinon density recovery](xxx-spinon-density.md) · [Building examples](building.md)

`XXXSpinonDensityRecovery` can resume a completed prefix of its declared
low/guard-precision controls. It uses the same density evaluator, arithmetic
screens and cross-control agreement test as an ordinary call. No new sampling
point is drawn. This is the numerical progress interface; persistent files,
producer identity and process ownership belong to the caller.

## Stop after a complete control

The [complete example](../examples/spinon_recovery_progress.cpp) makes a fresh
two-spinon calculation, cancels after its first control, and computes only the
remaining control on resume. It uses a native binary64/long-double pair when
the latter has a wider significand. Enable `BETHE_BUILD_EXAMPLES=ON`, then run
from your build directory:

```sh
cmake --build . --target bethe-spinon-recovery-example
./examples/bethe-spinon-recovery-example
```

The new overload is:

```cpp
auto partial = recovery.resume(beta_decimal_strings, {}, observer, stop_token);
auto finished = recovery.resume(beta_decimal_strings, partial, observer);
auto checked = recovery.restore(finished);
```

The observer receives a **const view** of the current result after each complete
low/guard control, before the stage's final cross-control assessment. Copy it
if it must outlive the callback. A callback exception propagates to the caller;
it is not converted to a numerical retry. The normal two-argument call remains
available and uses the same loop with an empty observer.

A complete control means both precisions have been evaluated and screened, or
a numerical failure has been recorded. Cancellation and wall limits are checked
between special-function kernels; they do not interrupt a running kernel.
An unfinished control can therefore require recomputation. It is not counted
as a completed control, and a failed control is not a zero density.

## Restore evidence, not an acceptance flag

Each result retains its full numerical policy and original decimal input
spellings. Changing the policy, budget or decimal inputs is an error on resume,
even when differently spelled inputs represent the same number. A generic
`recover(source, previous, observer, stop)` overload is also available, but
its factory must reproduce the same physical point at every precision and
across restarts; the library cannot inspect an arbitrary closure's identity.

`restore` validates the ordered control prefix, declared precision of every
stored scalar and amplitude, contour budgets, polarization dimensions,
failure codes, numerical diagnostics and monotonic consumed wall time.
It retains the low-precision value, log scale and arithmetic indicator at the
corresponding guard precision, allowing it to reconstruct the precision and
arithmetic screens. Requested tolerances are recomputed from the bound policy.

The overall status and estimate are rebuilt from this evidence. Every earlier
passing control still participates in the agreement check, including one
followed by a failed control. A later pair cannot erase a contradiction.
Histories continuing beyond an already terminal stage are rejected.

Already accepted, inconsistent or fully exhausted histories perform no new
numerical kernels. A partial cancelled history can continue with a new stop
token. An over-budget history remains a time-limit result with no estimate.
These checks establish consistency with the saved numerical evidence; they
do not independently recalculate form factors or provide rigorous error bounds.

## Account for time across process boundaries

`wall_seconds` is cumulative consumed time. Each completed control records
its cumulative completion time, which must not decrease or exceed the saved
total. A resumed call adds elapsed time to that total. Observer execution time
is included in subsequent budget checks. Saving a progress view does not capture
the remainder of a running control or time spent after that publication.

An external restart layer must therefore bind the producer, point identity,
policy and original deadline, and charge the elapsed interval after the last
publication before admitting more work. It may increase the saved total to
account for that interval; it must never reset it. Accepted results can be
validated and reused without new computation. If an operating-system process
dies inside a control, that unfinished control may be computed again if budget
remains. There is no exactly-once execution guarantee.

This interface does not itself lock files, serialize numbers, publish an atomic
snapshot or ensure power-loss durability. Those operations remain a separate
experimental adapter. It also does not change Async admission or cancellation:
an Async consumer must still drain admitted sibling tasks before reporting a
publication failure or cancellation.

Tests interrupt between controls and compare resumed results with fresh ones;
count source calls to detect repeated work; retain contradictions and numerical
failures; reject changed inputs, policies, diagnostics and budgets; and preserve
native binary128 and mixed-stage MPFR evidence. No private archive is required.
