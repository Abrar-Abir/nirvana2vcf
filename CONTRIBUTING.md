# Contributing to nirvana2vcf

Thanks for your interest in improving nirvana2vcf. Bug reports, feature requests,
and pull requests are all welcome.

## Reporting bugs

Please open an issue at https://github.com/Abrar-Abir/nirvana2vcf/issues with:

- A short description of the problem
- Steps to reproduce (a minimal Nirvana JSON snippet is ideal)
- Expected vs. actual VCF output
- nirvana2vcf version (`pip show nirvana2vcf`) and Python version

## Development setup

```bash
git clone https://github.com/Abrar-Abir/nirvana2vcf.git
cd nirvana2vcf
pip install -e ".[dev]"
pytest -v
```

## Pull requests

1. Fork the repo and create a feature branch off `main`.
2. Add tests for any new behavior or bug fix — `tests/` is laid out by module.
3. Run `pytest -v` locally and make sure everything passes.
4. Keep changes focused; unrelated refactors are easier to review as separate PRs.
5. Open a PR describing the change and the motivation.

## Code style

- Plain Python, no formatter is enforced — match the surrounding style.
- The runtime package depends only on `orjson`. New runtime dependencies should
  be discussed in an issue first.
- New CLI flags need both a unit test and a short note in the README.

## Validation

Concordance work against bcftools, VEP, and SnpEff lives under `validation/`.
See [validation/docs/strategy.md](validation/docs/strategy.md) for the
five-phase strategy if you want to extend it.
