# LiCSBAS2 tests

pytest-based test suite. No network access and no real data are needed;
everything runs on small synthetic datasets (see `synth.py`) in a few
tens of seconds.

## Running

```bash
python -m pip install pytest   # once, into the LiCSBAS conda env
python -m pytest               # whole suite (~30 s)
python -m pytest -m "not smoke"  # fast unit tests only (~5 s)
```

`LiCSBAS_lib` is put on the import path by `pytest.ini`/`conftest.py`,
so sourcing `bashrc_LiCSBAS.sh` is not required.

## Layout

- `test_tools_lib.py`, `test_inv_lib.py`, `test_loop_lib.py`,
  `test_io_lib.py`, `test_plot_lib.py` — unit tests of `LiCSBAS_lib`.
- `test_bin_smoke.py` (marker `smoke`) — runs steps 11–16 and
  `LiCSBAS_cum2vel.py` as subprocesses on a synthetic 10x10-pixel,
  5-epoch dataset with a known velocity field, and checks the inverted
  velocity against the truth.
- `synth.py` — builder of the synthetic GEOCml-format dataset.

## Known bugs pinned with xfail

Tests marked `xfail(strict=True)` document existing latent bugs (e.g.
calls to the nonexistent `np.linalg.leastsq` in
`LiCSBAS_inv_lib.censored_lstsq*`). They are expected to fail until the
bugs are fixed in separate PRs; when a fix lands, the test starts
XPASS-ing and must be updated to assert the correct behavior.

## Manual end-to-end tests

`bin/test_LiCSBAS.sh` and `bin/test2_LiCSBAS.sh` (if present) download a
real LiCSAR frame and run the full batch pipeline. They are
network-heavy manual tests, intentionally outside pytest.
