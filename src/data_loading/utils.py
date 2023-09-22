import functools
import pickle
from pathlib import Path


CACHE_DIR = Path(__file__).resolve().parents[2] / "cache"
DATA_DIR = Path(__file__).resolve().parents[2] / "data"


def cache_with_pickle(func):
    """NOTE: only works for functions without arguments"""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        cached_file = CACHE_DIR / (func.__name__ + ".pkl")
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
