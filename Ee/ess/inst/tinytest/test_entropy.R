# Test that the two entropy implementations give same results
expect_equal(ess:::joint_entropy(cars), ess:::joint_entropy2(cars))
