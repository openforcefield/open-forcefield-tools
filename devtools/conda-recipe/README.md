This is a recipe for building the current development package into a `conda` binary.

The installation on `travis-ci` is done by building the `conda` package, installing it, running the tests, and then if successful pushing the package to `anacond cloud` (and the docs to AWS S3, if desired).
The `anaconda` auth token is an encrypted environment variable generated using:
```bash
anaconda auth -n openforcefield-travis -o omnia --max-age 22896000 -c --scopes api:write
```
and then saved in the environment variable `BINSTAR_TOKEN`.

You can set up travis to store an encrypted token via
```bash
gem install travis travis encrypt BINSTAR_TOKEN=xx
```
where `xx` is the token output by `anaconda`.
The final command should print a line (containing `secure:`) for inclusion in your ``.travis.yml` file.
