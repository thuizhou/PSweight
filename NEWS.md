
\<!-NEWS.md is generated from NEWS.Rmd. Please edit that file –\>

# PSweight 1.2.0 (Release Date: 2024-03-27)

This release of PSweight version 1.1.8 includes several bug fixes and
enhancements aimed at improving the functionality and user experience.

1.Enhanced Parameter Support: The “PSmethod” and “OUTmethod” functions
have been updated to accurately incorporate parameters specified via
“ps.control” and “out.control”.

2.Propensity Score Model Output: Users can now access the fitted
propensity score model directly from the “SumStat” function output
through the “ps.fitObjects” component, facilitating deeper analysis and
diagnostics.

3.Improved Handling of Single Covariate Models: The package now includes
checks within “summary.SumStat()” and “plot.SumStat()” functions to
ensure they operate correctly even when the propensity score model is
based on a single covariate.

4.Flexible Treatment Variable Naming: The “SumStat_cl()” function has
been modified to replace the hardcoded ‘trt’ variable with a dynamic
reference, “data\[\[zname\]\]”, allowing users to specify their
treatment variable names.

These updates aim to enhance the PSweight package’s functionality,
usability, and analytical precision. Users are encouraged to explore the
new features and provide feedback for continuous improvement.
