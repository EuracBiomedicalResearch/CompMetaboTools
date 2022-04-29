# CompMetaboTools 0.4

## CompMetaboTools 0.4.2

- Pass `...` also to the `plot` call in `plot_pca`.

## CompMetaboTools 0.4.1

- Remove `plotOverlay` function (replaced by `plotChromatogramsOverlay` from
  `xcms`.

## CompMetaboTools 0.4.0

- Remove all feature grouping code (added to `xcms`).

# CompMetaboTools 0.3

## CompMetaboTools 0.3.3

- Add parameter `type` to `chromPeakArea`.

## CompMetaboTools 0.3.2

- Add `chromPeakArea` function.
- Add `groupFeatures,XCMSnExp,AbundanceSimilarityParam`.

## CompMetaboTools 0.3.1

- Add `groupFeatures,XCMSnExp,SimilarRtimeParam`.

## CompMetaboTools 0.3.0

- Move `groupClosest` and `groupSimilarityMatrix` to the `MsFeatures` package.


# CompMetaboTools 0.2

## CompMetaboTools 0.2.5

- Add parameter `mad` to `rsd` and `rowRsd` allowing to base the RSD calculation
  on the median absolute deviation instead of the standard deviation.

## CompMetaboTools 0.2.4

- Add `moreAreValidThan` to test for proportion of non-missing values in rows
  of a matrix.

## CompMetaboTools 0.2.3

- Add `joyPlot` to plot stacked (partially overlapping) chromatograms.

## CompMetaboTools 0.2.2

- Remove the `extractSpectraData` function as it was added to `MSnbase`.


# CompMetaboTools 0.1

## CompMetaboTools 0.1.0

- Change from `Chromatograms` to `MChromatograms` and from `Spectra` to
  `MSpectra` following the changes introduced in `MSnbase` 2.15.3.


## CompMetaboTools 0.0.8

- Fix issue in `.group_correlation_matrix` that would not group elements even if
  their correlation is larger than the threshold in some special cases.
- Add `featureGroupSpectra` function to extract a spectrum for each feature
  group which can be either a *pseudo spectrum* or a full scan MS1 spectrum.


## CompMetaboTools 0.0.7

- Add `matchRtMz` function to match features based on retention time and m/z 
  values.
- Fix `groupFeatures` for older xcms result objects that don't record MS level.
- Add `plotFeatureGroups` to draw feature groups into the m/z-rt space.

## CompMetaboTools 0.0.6

- Add the `plotOverlay` function to create a overlay plot of EICs.
- Use a less *greedy* grouping algorithm: create only groups of features with
  a correlation between each other `>= threshold`.
- `groupFeatures` adds a *process history* step to the `XCMSnExp` with the
  parameter object.
- Add `AbundanceCorrelationParam` method to group features based on feature
  abundances.
- Add `groupFeatures` method for `XCMSnExp` objects allowing to group features.
- Add `SimilarRtimeParam` and `EicCorrelationParam` feature grouping methods.

## CompMetaboTools 0.0.5

- Add `groupClosest` function to group values into small groups based on values
  being smaller than a threshold.
- Add `groupToSinglePolarityPairs` function.
- Add `c` method to (row-wise) combine `Chromatograms` objects.
- Add `groupByCorrelation` to allow grouping of rows within a numeric matrix
  by pairwise correlation with each other.
- Add `groupEicCorrelation` to allow grouping of (extracted ion) chromatograms 
  (EICs) based on pairwise correlation.

## CompMetaboTools 0.0.4

- Add `extractSpectraData` to support converting data from the `MSnbase` to the
  new `Spectra` package.

## CompMetaboTools 0.0.2

- Add `withinBatchFit`, `withinBatchAdjust` and `dropModels` function to perform
  per-batch separate model-based adjustments.

## CompMetaboTools 0.0.2

- Add `flag_*` functions to flag/identify potentially problematic model fits
  e.g. fitted to experimental data to explain an injection-order-dependent
  signal drift.

## CompMetaboTools 0.0.1

- Add `extract_time_stamp`.
- Add `plot_pca`.
- Add `rsd` and `rowRsd`.
- Add `sync_files_local` and `remove_local_files`.
