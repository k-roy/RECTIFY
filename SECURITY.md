# Security Policy

## Supported Versions

RECTIFY is an academic bioinformatics tool. We provide security fixes for
the most recent minor release and the previous minor release.

| Version | Supported |
|--------:|:---------:|
| 0.9.x   | yes       |
| < 0.9   | no (pre-public development) |

## Reporting a Vulnerability

If you discover a security issue (e.g. a crafted BAM that triggers code
execution, an arbitrary-file-read path traversal, or any other vulnerability
that goes beyond a routine bug), please **do not open a public GitHub
issue**.

Instead, email **kevinrjroy@gmail.com** with:

- A short description of the issue and its impact
- The minimum input or command that triggers it
- The RECTIFY version (`rectify --version`)

You should receive an acknowledgement within five business days. We will
work with you privately to confirm the issue, prepare a fix, and coordinate
disclosure. If the issue is non-exploitable but reproducible, we may ask to
discuss it in a public issue once we understand the scope.

## Scope

The maintainers consider the following in scope for security reports:

- Code execution from a crafted input file (BAM, FASTQ, GFF, FASTA)
- Path traversal or arbitrary-write via untrusted manifest TSVs / config YAMLs
- Pickle-based deserialization of untrusted files
- Credential or token leakage through logs / output files
- Denial of service triggered by inputs that pass header validation

The following are typically out of scope:

- Resource exhaustion from clearly oversized but well-formed inputs (a
  100 GB BAM will use a lot of memory — that is the user's
  responsibility to size around)
- Issues that require an attacker who can already modify the user's local
  filesystem or environment
- Issues in third-party dependencies (please report those upstream and
  optionally CC us so we can pin the version)
