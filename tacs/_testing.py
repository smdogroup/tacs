"""Internal test helpers shipped with TACS."""

from __future__ import annotations

import functools
from importlib.metadata import version

from packaging.version import Version


def fails_at_version(removalVersion: str):
    """Decorator that fails a unittest test method once TACS reaches ``removalVersion``.

    Use to self-decommission a test that covers a deprecated API. When the installed
    TACS version reaches the stated removal version, the test raises an
    ``AssertionError`` reminding the developer to delete both the deprecated code
    path and the test.

    Parameters
    ----------
    removalVersion : str
        Removal version, e.g. ``"3.14"``. Parsed with
        :class:`packaging.version.Version`.

    Examples
    --------
    >>> class TestDeprecated(unittest.TestCase):
    ...     @fails_at_version("3.14")
    ...     def test_old_api_still_warns(self):
    ...         ...
    """
    removal = Version(removalVersion)

    def decorator(method):
        @functools.wraps(method)
        def wrapper(self, *args, **kwargs):
            current = Version(version("tacs"))
            if current >= removal:
                raise AssertionError(
                    f"{type(self).__name__}.{method.__name__} is testing a "
                    f"deprecation that should be removed in TACS v{removal}, "
                    f"current version is v{current}. Please delete the "
                    f"deprecated API and this test."
                )
            return method(self, *args, **kwargs)

        return wrapper

    return decorator
