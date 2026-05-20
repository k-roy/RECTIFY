#!/usr/bin/env python3
"""
Portable path envelope for cross-cluster and cross-job sidecar correctness.

The problem
-----------
Cluster paths are unstable across runs:

- ``$L_SCRATCH`` on Sherlock is per-job, ephemeral (``/lscratch/<jobid>``).
- ``$TMPDIR`` on Hoffman2 is per-job (``/scratch/job/<jobid>``).
- ``$SCRATCH`` varies by user; ``$OAK`` varies by group; ``$HOME`` varies by
  cluster.
- ``$SAMPLE_OUTPUT`` (where rectify writes its results) may be moved by
  users between runs.

A sidecar that records ``/lscratch/13446575/foo.bam`` as the path of an
output will NEVER find the file on the next run, because ``$L_SCRATCH`` has
rotated to a new job id.

The fix
-------
Store paths as :class:`PortablePath` envelopes that capture BOTH a tokenized
template (``$L_SCRATCH/foo.bam``) and the literal absolute path at write
time (``cached_absolute``). On resolve, expand the template against the
CURRENT environment. If that doesn't exist on disk, fall back to
``cached_absolute`` and log a warning. If neither exists, raise.

Tokenization at write time tries, in order:
  1. Sample-relative (path is under ``sample_output_dir``).
  2. Env-relative, with longest-prefix-first env var lookup.
  3. Absolute (no template; just the literal absolute path).

Invariant for persistent outputs
--------------------------------
A PortablePath used as a NON-TRANSIENT output of a stage MUST NOT be
env_relative against ``$L_SCRATCH`` or ``$TMPDIR``. Those env vars are
job-ephemeral; persistent outputs must live in ``$SAMPLE_OUTPUT``
(sample-relative) or rooted under a stable env (``$OAK``, ``$SCRATCH``,
``$HOME``, ``$WORK``). :meth:`PortablePath.from_path` refuses to construct
such a value when ``transient=False`` (the default).

The corresponding test scenarios are in ``tests/test_portable_path.py``
(see Commit A.5 briefing §3.8).
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional, Union

logger = logging.getLogger(__name__)


# Env vars consulted at tokenization time. Order matters: a path under both
# ``$L_SCRATCH`` and ``$HOME`` (e.g. when ``$L_SCRATCH = $HOME/lscratch``)
# tokenizes against the LONGEST-prefix match, i.e. ``$L_SCRATCH``. So put
# more-specific / longer-prefix vars first.
#
# ``$SAMPLE_OUTPUT`` is NOT in this list — sample-relative is its own kind
# and is checked before env-relative.
_ENV_VAR_PRECEDENCE = (
    "L_SCRATCH",  # Sherlock per-job ephemeral (/lscratch/<jobid>)
    "TMPDIR",     # H2 per-job ephemeral (/scratch/job/<jobid>) or generic
    "SCRATCH",    # Sherlock per-user shared (/scratch/users/<user>)
    "OAK",        # Sherlock per-group persistent (/oak/stanford/groups/...)
    "WORK",       # H2/some clusters work dir
    "HOME",       # last resort; everyone has $HOME
)


# Env vars that are PER-JOB EPHEMERAL — paths rooted in these vanish when
# the job exits. Persistent (non-transient) outputs MUST NEVER tokenize
# against these.
_EPHEMERAL_ENV_VARS = frozenset({"L_SCRATCH", "TMPDIR"})


PortablePathKind = Literal["sample_relative", "env_relative", "absolute"]


class PathResolutionError(FileNotFoundError):
    """Raised when a :class:`PortablePath` cannot be resolved on the current host.

    Carries both attempted candidates (env-expanded template + cached_absolute)
    in the message so the caller can diagnose drift.
    """


@dataclass(frozen=True)
class PortablePath:
    """A path stored portably across runs.

    Three forms:

    - ``sample_relative``: relative to ``sample_output_dir``. Moves with
      the sample dir. The preferred form for stage outputs that live
      under the sample's working directory.
    - ``env_relative``: relative to an environment variable's value
      (e.g. ``$OAK/refs/genome.fa``). Resolves against the current
      environment on read. Falls back to ``cached_absolute`` if the env
      var is unset or the resolved path doesn't exist (e.g. the user
      moved the file).
    - ``absolute``: no template; just an absolute path. Survives only if
      the path stays put. Resume across machine moves will fail.

    Use :meth:`from_path` to construct; :meth:`resolve_now` to materialize.
    """

    # Stored verbatim in the sidecar JSON:
    template: str
    cached_absolute: str
    kind: PortablePathKind
    env_var: Optional[str] = None

    def __post_init__(self) -> None:
        # Schema sanity: env_var iff kind == env_relative
        if self.kind == "env_relative":
            if not self.env_var:
                raise ValueError(
                    f"PortablePath kind='env_relative' requires env_var; "
                    f"got template={self.template!r}, env_var={self.env_var!r}"
                )
            if not self.template.startswith("$"):
                raise ValueError(
                    f"PortablePath kind='env_relative' must have $-prefixed "
                    f"template; got template={self.template!r}"
                )
        else:
            if self.env_var is not None:
                raise ValueError(
                    f"PortablePath kind={self.kind!r} cannot have env_var; "
                    f"got env_var={self.env_var!r}"
                )

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_path(
        cls,
        p: Union[Path, str],
        *,
        sample_output_dir: Optional[Union[Path, str]] = None,
        transient: bool = False,
    ) -> "PortablePath":
        """Construct a PortablePath from an absolute (or to-be-absolutized) path.

        Tokenization order at construction time:

        1. Sample-relative — if ``p`` is under ``sample_output_dir``
           (after both are realpath'd, so symlinks don't cause misses).
        2. Env-relative — try each env var in :data:`_ENV_VAR_PRECEDENCE`
           in order; first match wins. Longest-prefix-first is achieved
           by ordering: ``L_SCRATCH`` > ``TMPDIR`` > ``SCRATCH`` > ``OAK``
           > ``WORK`` > ``HOME``.
        3. Absolute fallback — store the path as kind='absolute'.

        Args:
            p: the path to tokenize. May be relative; will be made absolute
               via :meth:`Path.resolve`.
            sample_output_dir: if given, paths under it tokenize as
               ``sample_relative``. Typically the stage's output directory.
            transient: if True, env_relative against ``$L_SCRATCH`` or
               ``$TMPDIR`` is allowed (used for per-region temp files).
               If False (default), env_relative against an ephemeral env
               var raises :class:`ValueError`.

        Returns:
            The most-specific tokenization possible.

        Raises:
            ValueError: if the result would be a non-transient env_relative
                against an ephemeral env var.
        """
        p = Path(p)
        if not p.is_absolute():
            p = p.resolve()
        # Follow symlinks for prefix matching. Keep the user-supplied form
        # in cached_absolute for diagnostic continuity.
        p_real = Path(os.path.realpath(p))
        cached_abs = str(p)

        # 1. Sample-relative.
        if sample_output_dir is not None:
            sample_resolved = Path(os.path.realpath(Path(sample_output_dir)))
            try:
                rel = p_real.relative_to(sample_resolved)
                return cls(
                    template=str(rel),
                    cached_absolute=cached_abs,
                    kind="sample_relative",
                    env_var=None,
                )
            except ValueError:
                pass  # not under sample dir; fall through

        # 2. Env-relative — try each candidate env var in precedence order.
        for var in _ENV_VAR_PRECEDENCE:
            val = os.environ.get(var)
            if not val:
                continue
            val_resolved = Path(os.path.realpath(Path(val)))
            try:
                rel = p_real.relative_to(val_resolved)
            except ValueError:
                continue
            # Match. Enforce invariant for persistent outputs.
            if not transient and var in _EPHEMERAL_ENV_VARS:
                raise ValueError(
                    f"PortablePath: refuse to construct non-transient "
                    f"env_relative path against ${var} (ephemeral). "
                    f"Persistent outputs must live in $SAMPLE_OUTPUT or a "
                    f"stable env-rooted dir ($OAK / $SCRATCH / $HOME / $WORK). "
                    f"If this is a per-region temp output (sub_outputs), pass "
                    f"transient=True. Offending path: {p!s}"
                )
            return cls(
                template=f"${var}/{rel}",
                cached_absolute=cached_abs,
                kind="env_relative",
                env_var=var,
            )

        # 3. Absolute fallback.
        return cls(
            template=str(p),
            cached_absolute=cached_abs,
            kind="absolute",
            env_var=None,
        )

    # ------------------------------------------------------------------
    # Resolution
    # ------------------------------------------------------------------

    def resolve_now(
        self,
        sample_output_dir: Optional[Union[Path, str]] = None,
    ) -> Path:
        """Resolve this PortablePath to a concrete Path against the current env.

        Resolution order:

        1. Expand by kind:

           - ``sample_relative``: ``(sample_output_dir / template).resolve()``
           - ``env_relative``: expand ``template`` via :func:`os.path.expandvars`,
             then resolve. If the env var is unset on the current host,
             expansion will leave the literal ``$VAR`` in the string — that
             is treated as "primary missing" and we fall through to step 2.
           - ``absolute``: ``Path(template).resolve()``.

        2. If the primary candidate doesn't exist on disk, fall back to
           ``cached_absolute`` and log a warning. If it exists, return it.

        3. If neither candidate exists, raise :class:`PathResolutionError`
           listing both candidates in the message.

        Args:
            sample_output_dir: required for ``kind='sample_relative'``.

        Returns:
            A :class:`Path` that EXISTS on disk.

        Raises:
            PathResolutionError: if neither candidate exists.
            ValueError: if ``kind='sample_relative'`` but ``sample_output_dir``
                is None.
        """
        primary: Optional[Path] = None

        if self.kind == "sample_relative":
            if sample_output_dir is None:
                raise ValueError(
                    f"PortablePath kind='sample_relative' requires "
                    f"sample_output_dir to resolve; got None. self={self!r}"
                )
            primary = (Path(sample_output_dir) / self.template).resolve()

        elif self.kind == "env_relative":
            assert self.env_var is not None  # __post_init__ guarantees this
            val = os.environ.get(self.env_var)
            if val:
                expanded = os.path.expandvars(self.template)
                # expandvars leaves unset $X as the literal "$X". If our env
                # var IS set but expansion still contains "$", something
                # weirder is going on (unrelated $-vars in the template).
                # Treat any remaining "$" as a failure to expand → fall back.
                if "$" not in expanded:
                    primary = Path(expanded).resolve()
            # else: env var unset → primary stays None → fall back to cached.

        elif self.kind == "absolute":
            primary = Path(self.template).resolve()

        else:
            raise ValueError(f"PortablePath: unknown kind={self.kind!r}")

        if primary is not None and primary.exists():
            return primary

        # Primary missing — try cached_absolute.
        cached = Path(self.cached_absolute)
        if cached.exists():
            if primary is not None:
                logger.warning(
                    "PortablePath drift: template %r resolved to %s "
                    "(does not exist); falling back to cached_absolute=%s",
                    self.template, primary, cached,
                )
            else:
                logger.warning(
                    "PortablePath: env var $%s unset on this host; "
                    "falling back to cached_absolute=%s",
                    self.env_var, cached,
                )
            return cached

        # Neither exists.
        raise PathResolutionError(
            f"PortablePath cannot resolve: "
            f"template={self.template!r} -> {primary!s} (missing); "
            f"cached_absolute={cached!s} (missing). "
            f"self={self!r}"
        )

    # ------------------------------------------------------------------
    # JSON serialization
    # ------------------------------------------------------------------

    def to_json(self) -> dict:
        """Serialize to a JSON-compatible dict for embedding in a sidecar.

        Omits ``env_var`` when it's None to keep the JSON tidy.
        """
        d: dict = {
            "template": self.template,
            "cached_absolute": self.cached_absolute,
            "kind": self.kind,
        }
        if self.env_var is not None:
            d["env_var"] = self.env_var
        return d

    @classmethod
    def from_json(cls, d: dict) -> "PortablePath":
        """Deserialize from a dict loaded from a sidecar.

        Raises:
            ValueError: if the dict is missing required keys or has an
                invalid structure (e.g. bare string passed instead of dict).
        """
        if not isinstance(d, dict):
            raise ValueError(
                f"PortablePath.from_json: expected dict, got "
                f"{type(d).__name__}={d!r}"
            )
        required = {"template", "cached_absolute", "kind"}
        missing = required - set(d.keys())
        if missing:
            raise ValueError(
                f"PortablePath.from_json: missing required keys "
                f"{sorted(missing)}; got {d!r}"
            )
        return cls(
            template=d["template"],
            cached_absolute=d["cached_absolute"],
            kind=d["kind"],
            env_var=d.get("env_var"),
        )


# ----------------------------------------------------------------------
# argv tokenization
# ----------------------------------------------------------------------


def tokenize_argv_paths(
    argv: list[str],
    sample_output_dir: Union[Path, str],
) -> list[str]:
    """Replace absolute paths in ``argv`` with their PortablePath templates.

    Used at sidecar-write time to build ``argv_template`` from raw ``argv``
    so that two runs that differ only in ``$L_SCRATCH`` resolution produce
    the same ``config_hash``.

    Tokens that don't start with ``/`` pass through unchanged. Tokens that
    look like absolute paths are tokenized; if the result is
    ``kind='absolute'`` (no env match, not under sample_output_dir), they
    pass through unchanged too (the absolute form IS the template). Other
    tokens (``sample_relative``, ``env_relative``) are replaced with their
    template.

    ``transient=True`` is passed to ``PortablePath.from_path`` for argv
    tokenization because we're recording the invocation, not declaring a
    persistent output — the ephemeral-env-var invariant should NOT fire here.

    Args:
        argv: the argv list (e.g. ``sys.argv``).
        sample_output_dir: the run's sample output directory.

    Returns:
        New list with absolute path tokens replaced by templates. Same
        length as input.
    """
    out: list[str] = []
    for arg in argv:
        if not arg.startswith("/"):
            out.append(arg)
            continue
        try:
            pp = PortablePath.from_path(
                arg,
                sample_output_dir=sample_output_dir,
                transient=True,
            )
        except Exception:
            # Pathological path — pass through.
            out.append(arg)
            continue
        if pp.kind == "absolute":
            out.append(arg)
        else:
            out.append(pp.template)
    return out
