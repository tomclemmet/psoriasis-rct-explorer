This is an R shiny app for a meta-analysis of trial data, reporting the data and meta-analysis results

The app/ directory holds the app script and the sqlite database. The R/ directory includes R scripts I will run manually to produce the sqlite database, check for mistakes in the database, and run the meta-analysis

Always assume scripts are run from the project root directory, no need for complex code resolving file paths

When editing a file in the 'R/' directory, don't be overly defensive. These scripts will only ever be run by me, the priority is readability

Don't include Claude as a co-author on commits

Don't waste tokens on checks I could do myself