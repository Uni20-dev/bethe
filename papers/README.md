# Research papers

[Overview](../README.md) | [Bibliography and provenance](../CITATIONS.md) |
[Model catalogue](../docs/models.md)

This is a small, curated archive of papers used to develop the solvers. The
complete bibliography remains in [data/citations.json](../data/citations.json);
entries below use its stable reference IDs. Archiving a paper does not imply
that all its methods or models are implemented.

For a connected pedagogical account alongside the individual papers, see
Jean-Sébastien Caux's [The Bethe Ansatz](https://integrability.org/).
Our [bibliography](../CITATIONS.md) also links the particular sections used
for equations, spinons, and state classification.

## Archived PDFs

### rylands-2022

Colin Rylands, Bruno Bertini, and Pasquale Calabrese,
*Integrable quenches in the Hubbard model*.

- [Local PDF](rylands-2022-arxiv-v1.pdf), unchanged arXiv:2206.07985v1.
- [Source and license evidence](https://arxiv.org/abs/2206.07985v1): the
  abstract page's **view license** link identifies
  [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
- [Exact download](https://arxiv.org/pdf/2206.07985v1); retrieved 2026-09-22.
- [Bibliography entry](../CITATIONS.md#rylands-2022), including the published
  article's DOI. This file is the preprint, not the publisher's edition.
- Relevant here: Hubbard Bethe equations and conventions. The repository
  does not implement the paper's quench dynamics.

The PDF is 1,112,865 bytes (28 pages). Its SHA-256 digest is recorded in
[SHA256SUMS](SHA256SUMS). From the repository root, verify archived bytes with:

```sh
sha256sum --check papers/SHA256SUMS
```

## Linked editions, not archived

These editions have arXiv distribution notices rather than an explicit license
allowing us to redistribute them. Keep the links; add a PDF only after finding
a suitable licensed edition or obtaining permission. This list records the
editions checked on 2026-09-22, not a claim that no other edition is available.

| Reference | Edition checked | Relevant topic |
| --- | --- | --- |
| [lieb-wu-2003](../CITATIONS.md#lieb-wu-2003) | [cond-mat/0207529v2](https://arxiv.org/abs/cond-mat/0207529v2) | Hubbard equations and ground state |
| [mei-2017](../CITATIONS.md#mei-2017) | [1609.08045v3](https://arxiv.org/abs/1609.08045v3) | Open XXX/XXZ equations |
| [vlijm-2016](../CITATIONS.md#vlijm-2016) | [1606.09516v2](https://arxiv.org/abs/1606.09516v2) | Periodic XXZ equations and excited states |
| [dugave-2015](../CITATIONS.md#dugave-2015) | [1412.8217v2](https://arxiv.org/abs/1412.8217v2) | Massive periodic XXZ equations |
| [grijalva-2019](../CITATIONS.md#grijalva-2019) | [1901.10932v4](https://arxiv.org/abs/1901.10932v4) | Massive open XXZ boundary roots |

The supplied root-level `0207529v2.pdf` is left untouched and untracked.
Additional personal reading copies can go in `papers/local/`, which Git ignores.
Neither location is part of the public archive; do not force-add personal copies
without checking redistribution permission.

## Adding a paper

1. Add or reuse its reference ID in `data/citations.json`, following the
   [citation guide](../docs/citations.md). Keep bibliographic corrections there.
2. Check permission for the **exact edition** being archived. Public download
   access alone is not enough: arXiv's
   [nonexclusive distribution license](https://arxiv.org/licenses/nonexclusive-distrib/1.0/license.html)
   grants distribution rights to arXiv, and its
   [older assumed license](https://arxiv.org/licenses/assumed-1991-2003/license.html)
   has the same distinction. Do not transfer a publisher edition's license to
   a different preprint without evidence.
3. Use a versioned filename such as `<reference-id>-arxiv-v2.pdf` or
   `<reference-id>-published.pdf`. Preserve the original bytes and notices;
   put reading notes in a separate Markdown document.
4. Add an entry here with author attribution, source, exact edition, license
   link and evidence, retrieval date, and whether the file was changed.
   Add its digest to `SHA256SUMS` using a repository-root-relative path.
5. Verify the file opens, the title/version matches, and the checksum passes.
   Review staged files explicitly so local reading copies are not included.

Papers retain their individual licenses and copyright notices; they are **not
relicensed under the code's GPL-3.0-or-later license**. No author endorsement is
implied. The PDFs are reading material, not build or runtime dependencies.

For this small collection, ordinary Git keeps PDFs available in a normal clone.
Avoid repeatedly replacing large binaries; revisit a separate archive or Git LFS
if the collection becomes substantial.
