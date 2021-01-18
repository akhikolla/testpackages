# regex pattern to validate
uuid_ptrn <- "[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}"
uuid_nil <- "00000000-0000-0000-0000-000000000000"


# test uuid_generate_random -----------------------------------------------

expect_error(uuid_generate_random(NULL))
expect_error(uuid_generate_random(integer(0)))
# expect_error(uuid_generate_random(NA_integer_))
expect_error(uuid_generate_random(c(0, 0)))
expect_equal(uuid_generate_random(0), character(0))
expect_true(is.character(uuid_generate_random(1)))
expect_true(grepl(uuid_ptrn, uuid_generate_random(1)))
expect_equal(length(uuid_generate_random(1)), 1)
expect_equal(length(uuid_generate_random(5)), 5)
expect_equal(length(unique(uuid_generate_random(1000))), 1000)


# uuid_generate_nil -------------------------------------------------------

expect_error(uuid_generate_nil(NULL))
expect_error(uuid_generate_nil(integer(0)))
# expect_error(uuid_generate_nil(NA_integer_))
expect_error(uuid_generate_nil(c(0, 0)))
expect_equal(uuid_generate_nil(0), character(0))
expect_true(grepl(uuid_ptrn, uuid_generate_nil(1)))
expect_true(is.character(uuid_generate_nil(1)))
expect_equal(uuid_generate_nil(1), uuid_nil)
expect_equal(length(unique(uuid_generate_nil(100))), 1)


# test uuid_generate_name -------------------------------------------------

expect_error(uuid_generate_name(NULL))
expect_equal(uuid_generate_name(character(0)), character(0))
expect_true(grepl(uuid_ptrn, uuid_generate_name(NA_character_)))
expect_true(grepl(uuid_ptrn, uuid_generate_name("")))
expect_true(grepl(uuid_ptrn, uuid_generate_name("a")))
expect_equal(length(uuid_generate_name(letters)), length(letters))
expect_equal(length(unique(uuid_generate_name(letters))), length(letters))


# test uuid_validate ------------------------------------------------------

expect_error(uuid_validate(NULL))
expect_equal(uuid_validate(character(0)), logical(0))
expect_true(uuid_validate(uuid_generate_random(1)))
expect_true(uuid_validate(uuid_generate_nil(1)))
expect_true(uuid_validate(uuid_generate_name("a")))
expect_false(uuid_validate(NA_character_))
expect_false(uuid_validate(""))
expect_false(uuid_validate("a"))
expect_equal(length(uuid_validate(c("a", "f"))), 2)
expect_equal(length(uuid_validate(uuid_generate_random(10))), 10)
