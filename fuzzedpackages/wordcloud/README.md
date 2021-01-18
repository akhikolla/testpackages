# wordcloud
Functionality to create pretty word clouds, visualize differences and similarity between documents, and avoid over-plotting in scatter plots with text.

# Usage

See: http://blog.fellstat.com/?cat=11


# Building and installing
Get the released version from CRAN:

```R
install.packages("wordcloud")
```

To build from this repository

```
R CMD build wordcloud
R CMD INSTALL wordcloud_*.*.tar.gz
```

Or simply

```R
# install.packages("devtools")
devtools::install_github("ifellows/wordcloud")
```