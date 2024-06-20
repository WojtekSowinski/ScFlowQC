<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-scflowqc_logo_dark.png">
    <img alt="nf-core/scflowqc" src="docs/images/nf-core-scflowqc_logo_light.png">
  </picture>
</h1>

## Introduction

(This is part of the submission for ny final year university project - "scFlow: Optimizing Software Engineering Practices for Modular and Extensible Bioinformatics Pipelines")

This repository contains a small pipeline replicating the functionality of the single-sample quality control stage of [the ScFlow pipeline](https://github.com/combiz/nf-core-scflow).
The code published here informed sections 5.1 and 6.2 and Chapter 8 of the project report and was used to test my implementation and refactorings of [ScFlowQC](https://github.com/WojtekSowinski/scflow-modules/tree/main/scflowqc).

## Usage

Please refer to [the scFlow pipeline repository](https://github.com/combiz/nf-core-scflow) for instruction on using the pipeline

To run the pipeline using the original ScFlow R package, checkout commit 52be7b1b2b56e138f68f6f4f57a5a464468288b2.

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
