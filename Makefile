test:
	eython3 tests/run_test.py

coverage:
	python3 -m coverage run tests/run_test.py
	@echo "computing coverage.."
	python3 -m coverage report > "unittest"
	python3 tests/uni_report.py -i "unittest" -o "uni_report"
	rm "unittest"
	@echo "check uni_report.."

package:
	python3 setup.py sdist
	rm -rf ANNOgesic.egg-info
	ls dist/*

build:
	python3 setup.py bdist

package_to_pypi:
	python setup.py sdist register upload
	@echo "Go to https://pypi.python.org/pypi/ANNOgesic/"

html_doc:
	cd docs && make html && cd ..

new_release:
	new_release:
	@echo "* Create/checkout a release branch"
	@echo "  git branch release_v0.3.X"
	@echo "  git checkout release_v0.3.X"
	@echo "* Change bin/reademption"
	@echo "* Change setup.py"
	@echo "* Change docs/source/conf.py"
	@echo "* Change CHANGELOG.txt"
	@echo "* Create new docs"
	@echo "* Test package creation"
	@echo "* Test doc creation"
	@echo "* make package_to_pypi"
	@echo "* git add CHANGELOG.txt bin/reademption docs/source/conf.py setup.py"
	@echo "* Commit changes e.g. 'git commit -m \"Set version to 0.3.X\"'"
	@echo "* Tag the commit e.g. 'git tag -a v0.3.X -m \"version v0.3.X\"'"
	@echo "* Merge release into dev and master"
	@echo "* Push it to github: git push"
	@echo "* Generate a new release based on this tag at"
	@echo "  https://github.com/konrad/READemption/releases/new"
	@echo "* Upload new docs using 'make upload_doc'"
