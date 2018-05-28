#!/bin/sh

LOGFILE=test.log

R --vanilla --slave > ${LOGFILE} 2>&1 <<'EOF'
library(doMC)
library(RUnit)

verbose <- as.logical(Sys.getenv('FOREACH_VERBOSE', 'FALSE'))

library(doMC)
registerDoMC()

options(warn=1)
options(showWarnCalls=TRUE)

cat('Starting test at', date(), '\n')
cat(sprintf('doMC version: %s\n', getDoParVersion()))
cat(sprintf('Running with %d worker(s)\n', getDoParWorkers()))

tests <- c('options.R')

errcase <- list()
for (f in tests) {
  cat('\nRunning test file:', f, '\n')
  t <- runTestFile(f)
  e <- getErrors(t)
  if (e$nErr + e$nFail > 0) {
    errcase <- c(errcase, t)
    print(t)
  }
}

if (length(errcase) == 0) {
  cat('*** Ran all tests successfully ***\n')
} else {
  cat('!!! Encountered', length(errcase), 'problems !!!\n')
  for (t in errcase) {
    print(t)
  }
}

cat('Finished test at', date(), '\n')
EOF
