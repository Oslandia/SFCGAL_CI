# SFCGAL Release Procedure

This document outlines the steps and checks necessary to perform a release of the SFCGAL project on GitLab.

## Steps to Perform a Release

### 1. Create a Tag and Release on GitLab
1. **Ensure All Tests Pass**:
   - Verify that all Continuous Integration (CI) tests pass on all services:
     - **GitLab CI** (Linux: Debian, Fedora, OpenSUSE. Windows docker)
     - **GitHub Actions** (Windows: MSYS2, vcpkg)
     - **Cirrus CI** (FreeBSD, macOs)

2. **Verify Dependencies**:

#### 1. Verify pysfcgal
    - After creating the release, ensure that `pysfcgal` (Python bindings for SFCGAL) functions correctly with the new SFCGAL version.
    - Run all relevant tests for `pysfcgal` and verify that they pass.

#### 2. Verify PostGIS Compatibility
    - Confirm that PostGIS is compatible with the new release of SFCGAL.
    - Run integration tests and ensure no issues arise.

3. **Update NEWS**
    - If necessary, update NEWS file and then edit it to make it more user readable. Especially, the `BREAKING CHANGE` needs to be reviewed carefully and often to be rewritten, including migration guide for instance.

4. **Create a Tag**:
   - Navigate to the repository on GitLab.
   - Go to the "Repository" > "Tags" section.
   - Click "New Tag".
   - Name the tag according to the version you are releasing (e.g., `v1.5.1`).
   - Provide a description for the release.

5. **Create a Release**:
   - Navigate to the "Repository" > "Releases" section.
   - Click "New Release".
   - Select the tag you created.
   - Fill in the release notes with a summary of changes and enhancements.
   - Publish the release.

### 4. Release pysfcgal
- Once `pysfcgal` has been verified to work with the new SFCGAL release, proceed to release `pysfcgal`.
- Follow the necessary steps to tag and release `pysfcgal`, similar to the process for SFCGAL.

### 5. Announce the Release
- Announce the new release on relevant platforms (mailing lists, forums, social media).
- Update any relevant documentation to reflect the new release.

## Checklist Before Releasing
- [ ] All CI tests pass on GitLab, GitHub, and Cirrus.
- [ ] `pysfcgal` works with the new SFCGAL version.
- [ ] PostGIS is compatible with the new SFCGAL version.
