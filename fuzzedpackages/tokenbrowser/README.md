[![Travis-CI Build Status](https://travis-ci.org/kasperwelbers/tokenbrowser.svg?branch=master)](https://travis-ci.org/kasperwelbers/tokenbrowser)

Tokenbrowser
============
Create text browsers from tokenized texts to visualize text analysis annotations


Quick guide
============
This quick and dirty howto shows how to visualizate tokens as full texts using the browser template. Given tokens in a data.frame format, the create_browser() function creates an html file in which the tokens are pasted back to texts. Furthermore, it supports ways to highlight or color scale tokens.

To start, we open tokenbrowser and view the SOTU (state of the union) demo data.

```{r}
library(tokenbrowser)

head(sotu_data$tokens)
head(sotu_data$meta)
```

The SOTU data consists of two data.frames: the tokens and the document meta data. Both data.frames have a "doc_id"" column, to match the tokens to the document meta. The only required columns are the doc_id and token column in the tokens data.frame. The meta data can optionally be used as additional information in the html browser.

First, we create a simple browser using the create_browser function. By default, the function stores the html as a temporary file. The function returns the url to this location, which can then conveniently be used to open this url in the browser using the browseURL function.

```{r}
url = create_browser(sotu_data$tokens, sotu_data$meta)
browseURL(url)
```

To color or otherwise annotate the words in the browser, we can first add html tags to the tokens. This can be done manually using the tag_tokens() function and its support functions (these will be explained in another manual). For the sake of convenience, there are several standard wrappers for coloring words. Here we demonstrate the highlighted_browser and colorscaled_browser functions.

The highlighted_browser function simply takes an additional argument, value. If value is a logical vector, it specifies which tokens to highlight. Alternatively, if it is a numeric vector with values between 0 and 1, it specifies how strongly words are highlighted (tokens with a value of 0 or NA will not be tagged). For this demo we highlight words based on the number of characters.

```{r}
highlight = nchar(as.character(sotu_data$tokens$token))
highlight = highlight / max(highlight)
highlight[highlight < 0.3] = NA
url = highlighted_browser(sotu_data$tokens, value = highlight, sotu_data$meta)
browseURL(url)
```

Next, the colorscaled_browser function can similarly color words, but instead of highlighting it uses a scale ranging from -1 to 1. This is for instance usefull to visualize sentiment words or the results of a wordscaling analysis. For now, we simply use the number of characters again. 

```{r}
scale = nchar(as.character(sotu_data$tokens$token))
scale = rescale_var(sqrt(scale), -1, 1)
url = colorscaled_browser(sotu_data$tokens, value = scale, meta=sotu_data$meta)
browseURL(url)
```

Finally, the categorical_browser can color words differently, for instance to indicate different topics. For this example, we highlight certain part-of-speech tags.

```{r}
category = match(sotu_data$tokens$pos, c('N','M','V'))
url = categorical_browser(sotu_data$tokens, category=category, meta=sotu_data$meta)
browseURL(url)
```

Customization
============

It is also possible to create custom browsers, or to just create the HTML tagged texts. 

To just get a list of HTML texts with the different coloring features shown above, see the functions: ?highlight_tokens, ?colorscale_tokens and ?category_highlight_tokens

To get more dirty, please see the documentation of the following functions:
To add html tags to tokens, the ?tag_tokens function can be used in combination with ?attr_style to set style attributes. To paste tokens together, the ?wrap_documents function can be used, or the more elaborate ?create_browser function.


