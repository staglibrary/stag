/**
 * Tests for the stag.h header. Includes configuration information about the
 * library.
 *
 * This file is provided as part of the STAG library and released under the MIT
 * license.
 */
#include <gtest/gtest.h>
#include <stag.h>

TEST(StagTest, Version) {
  EXPECT_EQ(stag::VERSION_MAJOR, 1);
  EXPECT_EQ(stag::VERSION_MINOR, 1);
  EXPECT_EQ(stag::VERSION_PATCH, 0);
}
