# test backend version
expect_true(is.numeric_version(ced_version()))
expect_true(ced_version() >= "2.2")

# test empty input
expect_error(ced_enc_detect(NA_integer_))
expect_error(ced_enc_detect(character(3), character(2), NULL))
expect_silent(ced_enc_detect(character(3), character(1), NULL))
expect_error(ced_enc_detect(character(3), NULL, character(2)))
expect_silent(ced_enc_detect(character(3), NULL, character(1)))
expect_identical(ced_enc_detect(NULL), character(0))
expect_identical(ced_enc_detect(raw()), character(0))
expect_identical(ced_enc_detect(character()), character(0))
expect_identical(ced_enc_detect(character(1)), NA_character_)
expect_identical(ced_enc_detect(NA_character_), NA_character_)

# test hints types
expect_error(ced_enc_detect(NA_character_, NA_integer_, NULL))
expect_error(ced_enc_detect(NA_character_, NULL, NA_integer_))
expect_error(ced_enc_detect(NA_character_, numeric(2), NULL))
expect_error(ced_enc_detect(NA_character_, NULL, numeric(2)))

# test ASCII encoding
expect_identical(ced_enc_detect("Hello"), "US-ASCII")
expect_identical(ced_enc_detect(c("Hello", "World")), c("US-ASCII", "US-ASCII"))

# test preserve names
expect_identical(ced_enc_detect(c(a = "test")), c(a = "US-ASCII"))

# test hints
expect_identical(ced_enc_detect("Hello", "ASCII", "EN"), "US-ASCII")
expect_identical(ced_enc_detect(c("Hello", "\u041f\u0440\u0438\u0432\u0435\u0442"), c("ASCII", "UTF-8"), c("EN", "RU")), c("US-ASCII", "UTF-8"))

test_file <- system.file("test.txt", package = "ced")
test_txt <-  read.dcf(test_file, all = TRUE)

# test UTF-8 in various languages
expect_identical(ced_enc_detect(test_txt[["English"]]), "US-ASCII")
expect_identical(ced_enc_detect(test_txt[["Italian"]]), "US-ASCII")
expect_identical(ced_enc_detect(test_txt[["Spanish"]]), "UTF-8")
expect_identical(ced_enc_detect(test_txt[["French"]]), "UTF-8")
expect_identical(ced_enc_detect(test_txt[["Czech"]]), "UTF-8")
expect_identical(ced_enc_detect(test_txt[["Russian"]]), "UTF-8")
expect_identical(ced_enc_detect(test_txt[["Ukrainian"]]), "UTF-8")
expect_identical(ced_enc_detect(test_txt[["Chinese"]]), "UTF-8")
expect_identical(ced_enc_detect(test_txt[["Japanese"]]), "UTF-8")
expect_identical(ced_enc_detect(test_txt[["Korean"]]), "UTF-8")
expect_identical(ced_enc_detect(test_txt[["Arabic(3)"]]), "UTF-8")

# test non UTF-8 enodings
expect_identical(ced_enc_detect(iconv(test_txt[["Russian"]], "UTF-8", "WINDOWS-1251")), "windows-1251")
expect_identical(ced_enc_detect(iconv(test_txt[["Russian"]], "UTF-8", "IBM866")), "IBM866")

# test raw input
expect_identical(ced_enc_detect(charToRaw("Hello")), "US-ASCII")
expect_identical(ced_enc_detect(charToRaw("\u041f\u0440\u0438\u0432\u0435\u0442")), "UTF-8")
expect_identical(ced_enc_detect(charToRaw("Hello"), "ASCII", "EN"), "US-ASCII")
expect_identical(ced_enc_detect(charToRaw("\u041f\u0440\u0438\u0432\u0435\u0442"), "UTF-8", "RU"), "UTF-8")
