# Individual Project for Computational Mathematics II - Code Files

> MATH2731 Individual Project, December 2025

## The Code

There code here is used to produce figures for the report, alongside animations for the video component

## Figures Produced

| Figure Number | Final Code File | Description |
| --- | --- | --- |
| 1 | [BasicUniformTorus](BasicUniformTorus.m) | Round torus, strong convergence |
| 2 (top) | [BasicUniformTorus](BasicUniformTorus.m) | Round torus, weak convergence due to low n |
| 2 (bottom) | [BasicUniformTorus](BasicUniformTorus.m) | Round torus, weak convergence due to low sigma |
| 3 | [EllipticUniformTorus](EllipticUniformTorus.m) | Torus with circular cross-section but elliptic birds-eye view |
| 4 | [ChainedElliptic](ChainedElliptic.m) | 5 Chained Tori, each with circular cross-section but elliptic birds-eye view, good convergence |
| 5 | [ChainedElliptic](ChainedElliptic.m) | 5 Chained Tori, each with circular cross-section but elliptic birds-eye view, bad convergence |
| 6 | [ChainedElliptic_NoMHForConsiderations](ChainedElliptic_NoMHForConsiderations.m) | 5 Chained Tori, each with circular cross-section but elliptic birds-eye view, bad convergence, chain considerations sampled separately randomly |

## Animations Produced

> In consecutive order of being shown

| Animation Clip | Final Code File | Description |
| --- | --- | --- |
| [torus_animation](torus_animation.mp4) | [BasicUniformTorus](BasicUniformTorus.m) | Forming a uniform torus |
| [torus_animation_elliptic](torus_animation_elliptic.mp4) | [EllipticUniformTorus](EllipticUniformTorus.m) | Forming a torus  with circular cross-section but elliptic birds-eye view |
| [torus_animation_elliptic_chained](torus_animation_elliptic_chained.mp4) |[ChainedElliptic](ChainedElliptic.m) | Forming 5 Chained Tori, each with circular cross-section but elliptic birds-eye view |

> **Warning.** The animations each have 450 frames, so will take long to create. If you want to not do an animation,the last lines of each of the 3 code files calls `plotter` - if the last argument is `true` change it to `false` to not have an animation be created.

> **Warning.** All code uses the *Mapping Toolbox* - this needs to be installed before running.

> **Remark.** *Repository made public on Day DD Month YYYY, at around XXpm (UTC+1)*
