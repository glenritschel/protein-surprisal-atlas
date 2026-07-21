.PHONY: env test download build_pilot score_pilot validate_scoring run_analysis

env:
	conda env create -f environment.yml

test:
	pytest tests/

download:
	python scripts/download_data.py --config config.yaml

build_pilot:
	python scripts/build_pilot.py --config config.yaml

score_pilot:
	python scripts/score_pilot.py --config config.yaml --method sampled_mask

validate_scoring:
	python scripts/score_pilot.py --config config.yaml --method exact --validation-only
	python scripts/validate_scoring.py --config config.yaml

run_analysis:
	python scripts/run_analysis.py --config config.yaml
