context("test all stemmers")

split_words <- function(text) {
  tolower(unlist(stringi::stri_split_boundaries(str = text, opts_brkiter = stringi::stri_opts_brkiter("word", skip_word_none = TRUE))))
}

# cat(paste(italian_stemmer(words), collapse = "\", \""))

test_that("French", {
  expect_equal(object = french_stemmer(words = c("tester", "testament", "chevaux", "aromatique", "personnel", "folle", "acheteuse")),
               expected = c("test", "testament", "cheval", "aromat", "personel", "fou", "achet"))
  expect_equal(object = french_stemmer(words = c("complète", "caisière", "homicides")),
               expected = c("complet", "caisier", "homicid"))
  text <- "La deuxième situation décrite a une portée plus large: lorsque le juge est compétent au titre de la convention, il devra statuer par défaut"
  words <- split_words(text = text)
  expect_equal(object = french_stemmer(words),
               expected = c("la", "deux", "situ", "decrit", "a", "une", "porte", "plu", "larg", "lorsqu", "le", "juge", "est", "competent", "au", "titr", "de", "la", "convention", "il", "devra", "statu", "par", "defaut"))
})


test_that("German", {
  expect_equal(object = german_stemmer(words = c("kinder")),
               expected = c("kind"))
  text <- "Ist das Gericht nach dem Übereinkommen zuständig, muss es ein Versäumnisverfahren durchführen, wenn und soweit das innerstaatliche"
  words <- split_words(text = text)
  expect_equal(object = german_stemmer(words),
               expected = c("ist", "das", "gericht", "nach", "dem", "ubereinkomm", "zustandig", "muss", "es", "ein", "versaumnisverfahr", "durchfuhr", "wenn", "und", "soweit", "das", "innerstaatlich"))
})

test_that("Spanish", {
  expect_equal(object = spanish_stemmer(words = c("perros")),
               expected = c("perr"))
  text <- "De no hacerse la declaración en ese plazo, sólo el juez es competente para registrarlo"
  words <- split_words(text = text)
  expect_equal(object = spanish_stemmer(words),
               expected = c("de", "no", "hacers", "la", "declaracion", "en", "ese", "plaz", "sólo", "el", "juez", "es", "competent", "para", "registrarl"))
})

test_that("Italian", {
  expect_equal(object = italian_stemmer(c("arrivederci")),
               expected = c("arrivederc"))
  text <- "Questa disposizione dovrebbe prescindere dal domicilio del convenuto in uno Stato vincolato dalla Convenzione e applicarsi in tutti i casi in cui sussista la competenza del giudice adito ai sensi della Convenzione"
  words <- split_words(text = text)
  expect_equal(object = italian_stemmer(words),
               expected = c("quest", "disposizion", "dovrebb", "prescinder", "dal", "domicil", "del", "convenut", "in", "uno", "stato", "vincolat", "dalla", "convenzion", "e", "applicars", "in", "tutti", "i", "casi", "in", "cui", "sussist", "la", "competenz", "del", "giudic", "adito", "ai", "sensi", "della", "convenzion"))
})

test_that("Portuguese", {
  expect_equal(object = portuguese_stemmer(c("adeus")),
               expected = c("adeu"))
  text <- "A segunda situação tem um âmbito mais lato. Quando o tribunal é competente segundo a Convenção, terá de prosseguir a instância"
  words <- split_words(text = text)
  expect_equal(object = portuguese_stemmer(words),
               expected = c("a", "segund", "situaca", "tem", "um", "ambit", "mai", "lato", "quand", "o", "tribunal", "é", "competent", "segund", "a", "convenca", "tera", "de", "prosseguir", "a", "instanci"))
})

test_that("Finnish", {
  expect_equal(object = finnish_stemmer(c("taivas")),
               expected = c("taiva"))
})

test_that("Swedish", {
  expect_equal(object = swedish_stemmer(c("stiga")),
               expected = c("stig"))
})

test_that("test empty word - numbers", {
  expect_equal(object = french_stemmer(words = c("")),
               expected = c(""))
  expect_equal(object = french_stemmer(words = c("123")),
               expected = c("123"))
  expect_equal(object = german_stemmer(words = c("")),
               expected = c(""))
  expect_equal(object = spanish_stemmer(words = c("")),
               expected = c(""))
  expect_equal(object = italian_stemmer(words = c("")),
               expected = c(""))
  expect_equal(object = portuguese_stemmer(words = c("")),
               expected = c(""))
  expect_equal(object = finnish_stemmer(words = c("")),
               expected = c(""))
  expect_equal(object = swedish_stemmer(words = c("")),
               expected = c(""))
})
