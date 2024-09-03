# FKF v0.2.5

- Fix to class export issues on R-devel
- Change url in README

# FKF v0.2.3

- Fix bug in Kalman Smoother
- Tidy up Kalman Smoother documentation to correct notation

# FKF v0.2.2

- Fix memory leak in fks identified by CRAN valgrind tests
- Fix building on MAC arm64 (M1) platform

# FKF v0.2.1

2021-09-28 Paul Smith

- Incorporated Kalman Smoother Code by Danial McDonald
- Outputs from the `fkf` function altered for use by the smoother. This
  includes a change in the output ordering. Code
  that references by name (e.g. `out$at`} should work as expected. Code
  that references by location (e.g. `out[[1]]`) might not.
- Implimentation of USE_FC_LEN_T on calls to BLAS and LAPACK

# FKF v0.1.7

2020-06-03 Paul Smith

- Fixed ASAN error caused by error in Vignette example
- Depreciated check.input option to fkf so dimensions always checked - the reason for ongoing ASAN issues since v0.1.2.

# FKF v0.1.6

2020-06-01 Paul Smith

- Tidying NAMESPACE to fix S4 class export warning with FKF v0.1.5 on
	all linux R-devel
- Converted Vignettes to use knitr

# FKF v0.1.5

2018-07-16 Paul Smith

- Correction of incorrect submission in v0.1.4 of the fix to the check.inputs = FALSE example
- Fixed simmilar ASAN and valgrind error in unit tests
- Stopped RUnit tests writing files which caused CRAN build errors due to directory permissions

# FKF v0.1.4

2018-07-01 Paul Smith

- section Maintainer updated
- Fixed check.inputs = FALSE example which was passing incorrect format variables resulting in ASAN and valgrind errors
- Moved function documentaion to use roxygen2 for generation 

# FKF v0.1.2

2014-01-03 Aleksandar Spasojevic

- sspir::kfilter reference removed because package sspir has been archived
- section Maintaner updated since ZHAW - Institute of Data Analysis and Process Design is maintaing FKF now

# FKF v0.1.2

2012-03-20 David Luethi

- RUnit added as dependancy


# FKF v0.1.1

2010-10-17 David Luethi

- Error in plot.fkf: Standardized residuals were calculated	wrongly.
- plot-method of fkf-objects with type == 'resid.qq' puts the number of the variable in the title (main = ...)
- fkf.Rd, plot.fkf.Rd: Documentation extended, examples added.
- R unit test added.

2010-12-19 Philipp Erb

- Fixed minor inconsistency in input check messages of fkf()

# FKF v0.1.0

2009-09-28 Philipp Erb

- NA-values in the vector of response variables yt are now supported. In the filter recursions the observation vector yt with missing values will be reduced to the vector yt- without missing values and the measurement equation will be adjusted accordingly.

# FKF v0.0.3

2009-07-30 David Luethi

- c-function 'cfkf': Warning for non-finite 'logDet' or 'mahalanobis' removed (too verbose as information is contained in	the output 'status' and also documented). Potential warnings occurring in the Cholesky factorization are more informative now.

# FKF v0.0.2

2009-06-26 David Luethi

- Minor corrections in the documentation.
- Function 'fkf' scans for missing arguments now.

# FKF v0.0.1

2009-02-10 David Luethi

- Release of FKF
