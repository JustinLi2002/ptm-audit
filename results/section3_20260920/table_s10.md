# Supplementary Table S10

Delta AUROC of the added channel under natural negatives, as a function of assumed annotation incompleteness. Threshold-sampled training, protein-disjoint evaluation, mean over three partitions.

### Language model channel

| PTM task | ceiling | burden at α=0 | ΔAUROC α=1 | λ=1.1 | λ=1.25 | λ=1.5 | ΔAUROC α=0 | retained |
|---|---|---|---|---|---|---|---|---|
| Methylation K/R | 1.16× | 2.6% | -0.2488 | -0.2199 | -0.2052† | -0.2052† | -0.2052 | 82% |
| Methylation R | 1.23× | 3.8% | -0.2238 | -0.1970 | -0.1674† | -0.1674† | -0.1674 | 75% |
| Acetylation K | 1.73× | 12.2% | -0.1776 | -0.1563 | -0.1289 | -0.0926 | -0.0666 | 37% |
| N-glycosylation | 1.30× | 5.0% | -0.0666 | -0.0525 | -0.0357† | -0.0331† | -0.0331 | 50% |
| Sumoylation K | 2.04× | 17.1% | -0.0663 | -0.0581 | -0.0473 | -0.0327 | -0.0090 | 14% |
| Phosphorylation S/T | 1.70× | 11.9% | -0.0266 | -0.0210 | -0.0122 | +0.0029 | +0.0147 | — |
| Ubiquitination K | 1.71× | 16.1% | -0.0114 | -0.0073 | -0.0013 | +0.0076 | +0.0141 | — |
| Phosphorylation Y | 1.14× | 2.3% | +0.0149 | +0.0193 | +0.0204† | +0.0204† | +0.0204 | — |

Mean over 7 units (methylation K and R collapsed): α=1 -0.0814, λ=1.1 -0.0692 (85%), λ=1.25 -0.0559 (69%), λ=1.5 -0.0448 (55%), α=0 -0.0351 (43%).

† target above this task's ceiling; the value shown is at α = 0, where the achieved inflation is the ceiling rather than the target.

### Interaction embedding channel

| PTM task | ceiling | burden at α=0 | ΔAUROC α=1 | λ=1.1 | λ=1.25 | λ=1.5 | ΔAUROC α=0 | retained |
|---|---|---|---|---|---|---|---|---|
| Methylation K/R | 1.16× | 2.6% | -0.2585 | -0.2314 | -0.2179† | -0.2179† | -0.2179 | 84% |
| Methylation R | 1.23× | 3.8% | -0.2817 | -0.2515 | -0.2184† | -0.2184† | -0.2184 | 78% |
| Acetylation K | 1.73× | 12.2% | -0.2616 | -0.2334 | -0.1985 | -0.1545 | -0.1242 | 47% |
| N-glycosylation | 1.30× | 5.0% | -0.0339 | -0.0271 | -0.0186† | -0.0176† | -0.0176 | 52% |
| Sumoylation K | 2.04× | 17.1% | -0.1196 | -0.1069 | -0.0909 | -0.0693 | -0.0365 | 30% |
| Phosphorylation S/T | 1.70× | 11.9% | -0.0152 | -0.0121 | -0.0078 | -0.0012 | +0.0034 | — |
| Ubiquitination K | 1.71× | 16.1% | -0.0294 | -0.0255 | -0.0204 | -0.0133 | -0.0086 | 29% |
| Phosphorylation Y | 1.14× | 2.3% | -0.0001 | +0.0028 | +0.0038† | +0.0038† | +0.0038 | — |

Mean over 7 units (methylation K and R collapsed): α=1 -0.1043, λ=1.1 -0.0919 (88%), λ=1.25 -0.0786 (75%), λ=1.5 -0.0672 (64%), α=0 -0.0568 (54%).

† target above this task's ceiling; the value shown is at α = 0, where the achieved inflation is the ceiling rather than the target.
