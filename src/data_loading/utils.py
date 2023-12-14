import functools
import pickle
from pathlib import Path
import inspect
import hashlib


CACHE_DIR = Path(__file__).resolve().parents[2] / "cache"
DATA_DIR = Path(__file__).resolve().parents[2] / "data"


def _get_default_args(func):
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }


def cache_with_pickle(func):
    """
    BUG: if you pass a keyword argument by position, this will crash
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs_passed):
        kwargs = _get_default_args(func)
        kwargs.update(kwargs_passed)
        cached_file = CACHE_DIR / (
            func.__name__
            + "_"
            + hashlib.md5(str(args).encode()).hexdigest()
            + "_"
            + hashlib.md5(str(kwargs).encode()).hexdigest()
            + ".pkl"
        )
        # TODO: this only works if the underlying function takes kwargs
        use_cache = kwargs.get("use_cache", True)
        if cached_file.exists() and use_cache:
            with open(cached_file, "rb") as f:
                print("reading from cache")
                ret = pickle.load(f)
            return ret
        else:
            ret = func(*args, **kwargs)
            with open(cached_file, "wb") as f:
                pickle.dump(ret, f)
            return ret

    return wrapper
