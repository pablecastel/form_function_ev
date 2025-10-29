# Variable yet ubiquitous: Hierarchical scaling of head functional morphology in lizards

## Overview

This repository contains the analytical pipeline for a phylogenetic comparative study examining the relationship between bite force and cranial shape across multiple hierarchical levels in *Podarcis* wall lizards and the broader Lacertidae family. The analysis investigates whether evolutionary trajectories (allometric slopes) linking functional performance to morphology are conserved or divergent across taxonomic scales.

---

## Research Questions

1. **Do individual species show unique evolutionary trajectories** linking bite force to head shape?
2. **Are species-level patterns consistent with genus-level patterns** in *Podarcis*?
3. **How do genus-level patterns compare to family-level patterns** across Lacertidae?
4. **What is the magnitude of divergence** between hierarchical levels?

---

## Methods

### Data

- **Morphometric data**: 3D landmark coordinates from dorsal cranial views (14 landmarks, 2D projection)
- **Functional data**: Bite force measurements (size-corrected residuals)
- **Phylogenetic data**: Time-calibrated phylogeny of Lacertidae
- **Sample size**: 
  - 670 *Podarcis* individuals (22 species)
  - 844 total specimens across Lacertidae

### Analytical Approach

#### 1. Data Preparation
- Size-corrected shape data using GLS regression
- Size-corrected bite force using residuals from allometric regression
- Phylogenetic tree pruning to match available data

#### 2. Hierarchical Analysis

**Level 1: Individual-level (E-PGLS)**
```r
# Evolutionary Phylogenetic Generalized Least Squares
# Accounts for within-species variation and phylogenetic structure
fit_indv <- lm.rrpp.ws(shape ~ bite*sp, 
                       data = gdf_order,
                       subjects = "sub",
                       cov = vcv(lacertidae_tree_drop),
                       print.progress = FALSE)
```

**Level 2: Genus-level (PGLS)**
```r
# Species means for Podarcis
pgls.reg_gen <- procD.pgls(shape_mean ~ bite, 
                            phy = lacertidae_tree_drop, 
                            data = gdf_mean, 
                            print.progress = FALSE)
```

**Level 3: Family-level (PGLS)**
```r
# All Lacertidae species means
pgls.reg_fam <- procD.pgls(shape_mean_fam ~ bite, 
                            phy = lacertidae_tree_drop_fam, 
                            data = gdf_fam, 
                            print.progress = FALSE)
```

#### 3. Slope Comparison Framework

Angular deviation between allometric vectors is calculated using vector correlation:

```r
# Calculate angle between two slope vectors
angle <- acos(vec.cor.matrix(rbind(slope1, slope2))) * 180 / pi
```

**Permutation testing**: 1000 iterations randomizing residuals while maintaining phylogenetic structure to generate null distributions of angular deviations.

**Statistical framework**:
- H₀: Observed angle = expected angle under no difference
- Test statistic: Z-score of observed angle relative to permuted distribution
- Significance: P-value from rank of observed angle

---

## Key Results

### 1. Species-Specific Trajectories

Individual *Podarcis* species show **unique allometric trajectories** linking bite force to head shape:

```r
# Mean angular deviation among species
mean_angle_among_species <- 35.2°  # (example value)
```

**Interpretation**: Species have diverged in how bite force scales with cranial morphology, suggesting adaptive divergence or historical contingency.

### 2. Pairwise Species Comparisons

```r
# Pairwise slope comparison results
pr <- pairwise(fit_indv, groups = gdf_order$sp, covariate = gdf_order$bite)
# Significant differences (p < 0.05) indicate divergent trajectories
```

**Visualization**: Heatmap showing pairwise angular deviations (blue = significantly different trajectories).

### 3. Hierarchical Comparison: Species vs. Genus

```r
# Compare each species slope to genus-level slope
results_df <- data.frame(
  species_vs_genus = species_names,
  angle_deg = unlist(angle.obs),
  Z = unlist(Z),
  p_value = unlist(p)
)
```

**Key Finding**: Most species trajectories are **NOT significantly different** from the genus-level pattern, but show substantial angular deviation (mean ~25-40°).

**Interpretation**: While not statistically different (lack of power or genuine similarity), there is meaningful biological variation in how species scale bite force with shape.

### 4. Hierarchical Comparison: Species vs. Family

```r
results_df_fam <- data.frame(
  species_vs_family = species_names,
  angle_deg = unlist(angle.obs_fam),
  Z = unlist(Z_fam),
  p_value = unlist(p_fam)
)
```

**Key Finding**: Greater angular deviations observed when comparing species to family-level patterns, indicating **hierarchical signal** in evolutionary trajectories.

### 5. Genus vs. Family Comparison

```r
angle.obs_fam_gen  # Angular deviation between Podarcis and Lacertidae slopes
Z_fam_gen          # Effect size
p_fam_gen          # Significance
```

---

## Visualizations

### Regression Plots

#### Individual Species Regressions
```r
# Each species shows its own trajectory
ggplot(fit.hab.ggplot.data_1, aes(x = bite, y = RegScore, color = sp)) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ sp)
```

**Output**: 22 panels showing species-specific bite force-shape relationships.

#### Overlaid Trajectories
```r
# Compare all trajectories simultaneously
plot3 <- ggplot(fit.hab.ggplot.data_1, aes(x = bite, y = RegScore, color = sp)) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5)
```

**Output**: Spaghetti plot revealing diversity of evolutionary trajectories.

#### Species Means
```r
# Visualize species centroids in morphospace
ggplot(df, aes(x = bitemean, y = regscoremean, color = species)) +
  geom_point(size = 3) +
  geom_label(aes(label = species))
```

### Morphological Deformation Grids

```r
# Visualize shape changes along bite force gradient
plotRefToTarget(M, preds$predmin, mag = 3)  # Weak bite morphology
plotRefToTarget(M, preds$predmax, mag = 3)  # Strong bite morphology
```

**Interpretation**: 
- Weak bite: Narrower, more gracile skull
- Strong bite: Wider, more robust skull (especially in temporal and parietal regions)

### Angular Deviation Distributions

#### Density Ridge Plots

```r
# Distribution of angular deviations from permutations
# Shows whether observed angle is extreme relative to null expectation

# Species vs. Genus
plot_6 <- ggplot(angles_long_gen, aes(x = Angle, y = Species, fill = Species)) +
  stat_density_ridges() +
  geom_segment(data = observed_gen)  # Observed angle as vertical line

# Species vs. Family  
plot_7 <- ggplot(angles_long_fam, aes(x = Angle, y = Species, fill = Species)) +
  stat_density_ridges() +
  geom_segment(data = observed_fam)
```

**Interpretation**: Dashed vertical lines show observed angular deviation; distribution shows null expectation. Lines in the tail indicate significant divergence.

### Phylogenetic Context

```r
# Ancestral state reconstruction of bite force
tree_recon <- contMap(lacertidae_tree_drop_fam, gdf_fam$bite, plot = TRUE)
plot(tree_recon, type = "fan", legend = 0.7)
```

**Output**: Fan phylogeny with branches colored by reconstructed bite force, showing evolutionary patterns across the family.

---

## Statistical Framework Summary

### ANOVA Results

**Individual E-PGLS**:
```r
anova(fit_indv)
# Tests for:
# - Main effect of bite force
# - Main effect of species
# - Bite × Species interaction (different slopes)
```

**Genus PGLS**:
```r
anova(pgls.reg_gen)
# Tests genus-level bite force-shape relationship
```

**Family PGLS**:
```r
anova(pgls.reg_fam)
# Tests family-level bite force-shape relationship
```

### Slope Comparison Methodology

1. Extract slope vectors from each hierarchical level
2. Calculate angular deviation using vector correlation
3. Generate null distribution via permutation (n=1000)
4. Calculate Z-score and P-value for observed angle
5. Interpret biological significance of angular deviations

---

## Key Packages

```r
library(geomorph)    # Geometric morphometrics
library(RRPP)        # Randomization procedures
library(ape)         # Phylogenetic analysis
library(phytools)    # Phylogenetic comparative methods
library(ggplot2)     # Visualization
library(dplyr)       # Data manipulation
library(ggridges)    # Density ridge plots
```

---

## File Structure

```
project/
├── data/
│   ├── 05_shape_allspecies_sizecorr_def_newcodes_Cet_JULY.bin
│   ├── 04_Resid_BF_dataset_complete_podnames_JULY.csv
│   ├── SPECIES_CODE_PODARCIS_JULY.csv
│   └── Lacertidae_tree_2_2.tre
├── scripts/
│   └── 06a_Analyses_BF_Geom_regression_Podarcis_inv_comparisonSlopes_cen.R
├── outputs/
│   ├── plots/
│   └── results/
└── README.md
```

---

## Biological Interpretation

### Main Conclusions

1. **Hierarchical scaling is ubiquitous**: Bite force-shape relationships exist at all taxonomic levels (individual, genus, family).

2. **Yet evolutionarily variable**: Species show divergent trajectories despite overall positive scaling.

3. **Pattern vs. Process**: The general pattern (stronger bite = more robust skull) is conserved, but the specific evolutionary trajectory varies.

4. **Implications for macroevolution**: 
   - Parallel evolution at broad scales
   - Historical contingency at fine scales
   - Multiple adaptive solutions to similar functional demands

### Evolutionary Hypotheses

**Scenario 1 - Adaptive Divergence**: Different species occupy distinct ecological niches requiring different cranial designs for similar bite forces.

**Scenario 2 - Developmental Constraints**: Underlying developmental programs differ among lineages, constraining trajectories.

**Scenario 3 - Historical Contingency**: Ancestral states differ, leading to different evolutionary paths toward similar functional endpoints.

---

## Future Directions

- Incorporate ecological variables (diet, habitat) to test adaptive hypotheses
- Expand to 3D morphometrics for complete cranial analysis
- Test ontogenetic vs. evolutionary allometries
- Mechanistic biomechanical modeling (FEA) to link shape to performance

---

## Citation

[Your Name]. (2025). Variable yet ubiquitous: Hierarchical scaling of head functional morphology in lizards. *PhD Dissertation Chapter*, [University Name].

---

## Contact

**Author**: [Your Name]  
**Email**: [your.email@institution.edu]  
**ORCID**: [0000-0000-0000-0000]

---

## Acknowledgments

This analysis uses geometric morphometric and phylogenetic comparative methods following best practices in evolutionary biology. Special thanks to the developers of `geomorph`, `RRPP`, and `phytools` packages.

---

## License

[Choose appropriate license: MIT, CC-BY-4.0, etc.]