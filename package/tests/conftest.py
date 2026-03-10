import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )
    parser.addoption(
        "--ccp4", action="store_true", default=False, help="run CCP4 only tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")
    config.addinivalue_line("markers", "ccp4: mark test as for CCP4 only")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--runslow"):
        skip_slow = pytest.mark.skip(reason="need --runslow option to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)

    if not config.getoption("--ccp4"):
        skip_ccp4 = pytest.mark.skip(reason="need --ccp4 option to run")
        for item in items:
            if "ccp4" in item.keywords:
                item.add_marker(skip_ccp4)