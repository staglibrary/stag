/**
 * Tests for the stag.h header. Includes configuration information about the
 * library.
 */
#include <gtest/gtest.h>
#include <stag.h>

TEST(StagTest, Version) {
  EXPECT_EQ(stag::VERSION_MAJOR, 0);
  EXPECT_EQ(stag::VERSION_MINOR, 1);
  EXPECT_EQ(stag::VERSION_PATCH, 6);
}
