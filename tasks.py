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

REPO_ROOT = Path(__file__).parent
SETUP_FILE = REPO_ROOT / "setup.py"
TEST_DIR = REPO_ROOT / "tests"
DOCS_DIR = REPO_ROOT / "docs"
DOCS_BUILD_DIR = DOCS_DIR / "_build"
SOURCE_DIR = REPO_ROOT / "tdsr"
TOX_DIR = REPO_ROOT / ".tox"
COVERAGE_FILE = REPO_ROOT / ".coverage"
COVERAGE_DIR = REPO_ROOT / "htmlcov"
COVERAGE_REPORT = COVERAGE_DIR / "index.html"
PYTHON_FILES = list(REPO_ROOT.glob("**/*.py"))


@task(help={"check": "Checks if source is formatted without applying changes"})
def format(c, check=False):
    """Format code"""
    python_files = " ".join([str(f) for f in PYTHON_FILES])
    black_options = "--diff" if check else ""
    c.run(f"pipenv run black {black_options} {python_files}")


@task
def lint(c):
    """Lint code"""
    c.run(f"pipenv run flake8 {SOURCE_DIR}")


@task
def test(c, min_coverage=None, parallel=False, fail_fast=True):
    """Run tests"""
    pytest_options = []
    if min_coverage:
        pytest_options.append(f"--cov-fail-under={min_coverage}")
    if parallel:
        pytest_options.append("-n auto")
    if fail_fast:
        pytest_options.append("-x")

    c.run(f"pipenv run pytest --force-sugar {' '.join(pytest_options)} {TEST_DIR}")


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
        c.run(f"pipenv run {provider}")
    else:
        # Build a local report
        c.run(f"pipenv run coverage html -d {COVERAGE_DIR}")
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
