[![Maven Central](https://maven-badges.herokuapp.com/maven-central/org.broadinstitute/gatk-fermilite-jni/badge.svg)](https://maven-badges.herokuapp.com/maven-central/org.broadinstitute/gatk-fermilite-jni)
[![Build Status](https://travis-ci.org/broadinstitute/gatk-fermilite-jni.svg?branch=master)](https://travis-ci.org/broadinstitute/gatk-fermilite-jni)

# gatk-fermilite-jni
JNI code for fermi lite.

This project builds dynamic libraries with a JNI API.
It allows Java code to call Heng Li's fermi lite assembler.

To build you'll need gmake, git, gcc, and Java 8.

To build and install a snapshot locally:

```
./gradlew install
```

This will work for testing but will only include a native library for your system.

To upload a snapshot from a Broad Institute OSX machine with both OSX and Linux binaries:
```
commit your changes and push your branch to github
scripts/build_both_dylib_and_so.sh
./gradlew uploadArchives printVersion
```

To upload to maven central
```
commit your changes and push your branch to github
# might have to do these to be able to sign the tagged version
# git config --global commit.gpgsign true
# gpg --list-keys # to get your key ID
# git config --global user.signingkey <your key ID>
# export GPG_TTY=$(tty)
# git config --global gpg.program gpg2
git tag -a -s 1.0.0-rc6 -m 'tag comment' # or some similarly formatted version
scripts/build_both_dylib_and_so.sh # build on Linux and move lib to src/main/c
# source artifactory credentials
./gradlew uploadArchive -Drelease=true
# go to https://oss.sonatype.org/#stagingRepositories and close and release your artifact
```

To use this JNI binding on another architecture for which we don't provide a binary:
  Go into ```src/main/c```.
  Modify the Makefile to produce a library name appropriate to your system.
  Type ```make``` (you'll need gmake, git, and gcc).
  Move the library you built somewhere permanent on your machine.
  Use ```-DLIBFML_PATH=<that permanent location>``` when you run GATK (or other Java program).
