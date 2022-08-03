"""
Tasks for maintaining the project.

Execute 'invoke --list' for guidance on using invoke
"""
import shutil
import pprint

from invoke import task
import webbrowser
from pathlib import Path
import utils

Path().expanduser()

ROOT_DIR = Path(__file__).parent
SETUP_FILE = ROOT_DIR / "setup.py"
TEST_DIR = ROOT_DIR / "tests"
DOCS_DIR = ROOT_DIR / "docs"
DOCS_BUILD_DIR = DOCS_DIR / "_build"
SOURCE_DIR = ROOT_DIR / "tdsm"
TOX_DIR = ROOT_DIR / ".tox"
COVERAGE_FILE = ROOT_DIR / ".coverage"
COVERAGE_DIR = ROOT_DIR / "htmlcov"
COVERAGE_REPORT = COVERAGE_DIR / "index.html"
PYTHON_DIRS = [str(d) for d in [SOURCE_DIR, TEST_DIR]]


@task(help={"check": "Checks if source is formatted without applying changes"})
def format(c, check=False):
    """Format code"""
    python_dirs_string = " ".join(PYTHON_DIRS)
    black_options = "--diff" if check else ""
    c.run("pipenv run black {} {}".format(black_options, python_dirs_string))
    isort_options = "{}".format("--check-only" if check else "")
    c.run("pipenv run isort {} {}".format(isort_options, python_dirs_string))


@task
def lint(c):
    """Lint code"""
    c.run("pipenv run flake8 {}".format(SOURCE_DIR))


@task
def test(c, min_coverage=None):
    """Run tests"""
    pytest_options = "--cov-fail-under={}".format(min_coverage) if min_coverage else ""
    c.run("pipenv run pytest --cov={} {}".format(SOURCE_DIR, pytest_options))


@task
def docs(c, target="html", serve=True):
    """Generate documentation and serve locally"""
    c.run(f"make -C {DOCS_DIR} {target}")
    if serve and target == "html":
        utils.serve_dir(DOCS_BUILD_DIR / "html")


@task
def type_check(c):
    """Check types"""
    c.run("pipenv run mypy")


@task
def install_hooks(c):
    """Install pre-commit hooks"""
    c.run("pipenv run pre-commit install -t pre-commit")
    c.run("pipenv run pre-commit install -t pre-push")


@task
def pre_commit(c):
    """Run all pre-commit checks"""
    c.run("pipenv run pre-commit run --all-files")


@task(
    pre=[test],
    help=dict(
        publish="Publish the result (default False)",
        provider="The provider to publish (default codecov)",
    ),
)
def coverage(c, publish=False, provider="codecov"):
    """Create coverage report"""
    if publish:
        # Publish the results via provider (e.g. codecov or coveralls)
        c.run("pipenv run {}".format(provider))
    else:
        # Build a local report
        c.run("pipenv run coverage html -d {}".format(COVERAGE_DIR))
        webbrowser.open_new_tab(COVERAGE_REPORT.as_uri())


@task
def clean_docs(c):
    """Clean up documentation build files"""
    utils.delete_dir(DOCS_BUILD_DIR)


@task
def clean_build(c):
    """Clean up files from package building"""
    c.run("rm -fr build/")
    c.run("rm -fr dist/")
    c.run("rm -fr .eggs/")
    c.run("find . -name '*.egg-info' -exec rm -fr {} +")
    c.run("find . -name '*.egg' -exec rm -f {} +")


@task
def clean_python(c):
    """Clean up python file artifacts"""
    c.run("find . -name '*.pyc' -exec rm -f {} +")
    c.run("find . -name '*.pyo' -exec rm -f {} +")
    c.run("find . -name '*~' -exec rm -f {} +")
    c.run("find . -name '__pycache__' -exec rm -fr {} +")


@task
def clean_tests(c):
    """Clean up files from testing"""
    utils.delete_file(COVERAGE_FILE)
    utils.delete_dir(TOX_DIR)
    utils.delete_dir(COVERAGE_DIR)


@task(pre=[clean_build, clean_python, clean_tests, clean_docs])
def clean(c):
    """Runs all clean sub-tasks"""
    pass
