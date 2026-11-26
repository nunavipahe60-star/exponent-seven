"""Recover the known field of discriminant -5923 with the divisor engine."""
from contextlib import redirect_stdout
from io import StringIO
from unittest.mock import patch

import pair_check


def main():
    stream = StringIO()
    with patch.dict(pair_check.BANDS, {1: (5923, 5923, 19, 29)}), redirect_stdout(stream):
        counts = pair_check.run(1)
    witness = "HIT 1 19 29 758 3413 13128 2243 5923"
    if counts[3] != 1 or stream.getvalue().splitlines() != [witness]:
        raise ValueError("The two-prime positive control for -5923 failed")
    print(witness)


if __name__ == "__main__":
    main()
