Contributing to STAG C++
=========================

When developing on the STAG project, please follow the guidelines below.

1. Every feature or bug should have a corresponding Github issue.
2. Create a branch for each issue you work on.
3. Develop the fix or feature locally.
4. Submit a pull request and state in the request which issue it corresponds to.

Please make sure you've addressed all of the following points before submitting a pull
request.

1. The code is clean and well commented.
2. You should add unit tests to test the newly developed code.
    - Include tests for invalid inputs. Any new method should include checks for
      invalid arguments and raise appropriate exceptions.
3. You should add an entry to the CHANGELOG to document the new feature or fix.
4. Make sure that any new methods are correctly documented in the doxygen-generated
documentation.
