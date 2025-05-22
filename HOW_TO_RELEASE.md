# How to release

- edit the CHANGELOG.md. The best way is to start with commitizen for that:

```bash
cz changelog --incremental --unreleased-version v2.0.0
```

and then edit it to make it more user readable. Especially, the `BREAKING
CHANGE` needs to be reviewed carefully and often to be rewritten, including
migration guide for instance.
- edit the version in `CMakeLists.txt`
- edit the version in `sonar-project.properties`
- create a merge request with these changes
- once it is merged, create a tagged release on gitlab.


What to check after the release:

- [TODO] the docker hub page
- [TODO] the new tag in the online documentation
