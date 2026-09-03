# LiCSBAS2 tests

pytest-based test suite. No network access and no real data are needed;
everything runs on small synthetic datasets (see `synth.py`) in a few
tens of seconds.

## Running

```bash
conda install -c conda-forge pytest   # once, into the LiCSBAS conda env
python -m pytest                      # whole suite (~30 s)
python -m pytest -m "not smoke"       # fast unit tests only (~5 s)
```

`LiCSBAS_lib` is put on the import path by `pytest.ini`/`conftest.py`,
so sourcing `bashrc_LiCSBAS.sh` is not required.

CI builds `LiCSBAS.yml` (the same environment users get) plus pytest, on
every python version in the matrix in `.github/workflows/ci.yml`. Because
`LiCSBAS.yml` sets only lower bounds, CI installs current releases and so
also acts as an early warning when a new upstream release breaks LiCSBAS.

## Layout

- `test_tools_lib.py`, `test_inv_lib.py`, `test_loop_lib.py`,
  `test_io_lib.py`, `test_plot_lib.py` — unit tests of `LiCSBAS_lib`.
- `test_bin_smoke.py` (marker `smoke`) — runs steps 11–16 and
  `LiCSBAS_cum2vel.py` as subprocesses on a synthetic 10x10-pixel,
  5-epoch dataset with a known velocity field, and checks the inverted
  velocity against the truth.
- `synth.py` — builder of the synthetic GEOCml-format dataset.

## Known bugs pinned with xfail

Latent bugs found while writing tests are pinned with
`xfail(strict=True)` markers: such a test is expected to fail until the
bug is fixed in a separate PR; when a fix lands, the test starts
XPASS-ing and must be updated to assert the correct behavior. There are
currently no open xfail tests.
