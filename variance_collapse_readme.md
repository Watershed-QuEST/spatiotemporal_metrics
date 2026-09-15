# Running variance collapse

### 1. Load the helpers

Run from your repository’s main folder:

```r
source("helpers/variance_breakpoints_helper.R")
source("helpers/variance_collapse_fx_helper.R")
```

PELT detection requires `changepoint`. Install it once if needed:

```r
install.packages("changepoint")
```

### 2. Prepare your data table

Use one row per site and sampling event, with numeric concentration and
**cumulative contributing watershed area** columns. Resolve duplicates and
quality flags first.

### 3. Detect breakpoints

Replace the column names below with yours:

```r
segmentation <- detect_variance_breakpoints_by_group(
  data = samples,
  concentration = c(NO3 = "nitrate", SRP = "phosphate"),
  area = "contributing_area",
  event = "Event",
  watershed = "watershed_id",
  min_sites = 6,
  min_segment_length = 4
)
```

- Use one concentration column or several named columns.
- Omit `watershed` if your data contain one watershed.
- Omit `event` if your data contain one sampling event.

The helper removes non-finite concentration/area pairs, orders observations by
area, and runs PELT separately for each group. Keep its output unchanged for the
next step.

### 4. Calculate variance collapse

```r
collapse <- variance_collapse_by_group(segmentation)
```

### 5. View and save results

```r
collapse$by_event

write.csv(
  collapse$by_event,
  "variance_collapse_results.csv",
  row.names = FALSE
)
```

Each row represents one watershed, event, and constituent:

| Result | Meaning |
|---|---|
| `threshold` | Lower bounding area of the first variance decrease |
| `threshold_next_area` | Upper bounding area |
| `threshold_midpoint` | Midpoint between those areas |
| `status` | Whether collapse was detected, no decrease was found, or detection was skipped |

Thresholds are `NA` when no collapse is identified.

To inspect one group’s breakpoints and variance comparisons:

```r
segmentation$groups             # Identifies groups in result order
segmentation$results[[1]]$breakpoints
collapse$results[[1]]$changepoints
```

A skipped group has a `NULL` collapse result. Save both complete outputs if you
want to retain detection settings and detailed results:

```r
saveRDS(
  list(segmentation = segmentation, collapse = collapse),
  "variance_collapse_details.rds"
)
```
