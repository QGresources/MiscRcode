# MiscRcode
A place to store miscellaneous, useful R code. Provide a brief description of the R code below.

* ManhattanPlot.R: A function to make Manhattan Plots. This was modified from Stephen Turner's function. Y-axis is the -log10 p value.

* ManhattanBeta.R: A function to make Manhattan Plots that show the absolute value of SNP effects. This was modified from the code above. Y-axis is the |Beta|.

* SplitCSV.sh: This code splits a csv with many columns into multiple smaller files with fewer columns. 
      To run: `./SplitCSV.sh <.csv file> <number of columns in output>`
